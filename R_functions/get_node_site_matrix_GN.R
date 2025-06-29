# Required packages.
require(ape)
require(phytools)
require(phangorn)

################################################################################
# getAncestorPath
################################################################################
# Get the sequence of ancestors of a node or tip of a phylogeny.
# `tree`: a phylogenetic tree in phylo format.
# `node`: a node/tip label or node/tip index.
# return: a vector with a sequence of heights of ancestor nodes of the focal
# node. The names of elements of this vector corresponds to the indexes of the
# ancestor nodes.
getAncestorPath <- function(tree, node){

  # Verify if the focal node informed is in an appropriate format.
  if(is.character(node)){
    node_index <- which(tree$tip.label==node)
  }
  else if(is.numeric(node)){
    node_index <- node
  }
  else{
    print("Error: node format unrecognized. Please, use node label or node index.")
  }

  # Get the indexes of the ancestor nodes of the focal node.
  anc <- nodepath(tree, from=node_index, to=(tree$Nnode+2))

  # Create a vector to get the heights of the ancestor nodes.
  node_heights <- rep(NA, length(anc))
  # Set the indexes of ancestor nodes as names of the vector elements.
  names(node_heights) <- as.character(anc)

  # Get the heights of the ancestor nodes and sum it with the height of the
  # root.
  for(i in 1:length(anc)){
    node_heights[i] <- tree$root.edge + phytools::nodeheight(tree, node=anc[i])
  }

  # Return the vector with the indexes of ancestor nodes an its respective
  # heights.
  return(node_heights)
}

################################################################################
# mapAncestorSpecies
################################################################################
# This function map the ancestor species index in their correspondent species
# names at that time.
# This function work based in the fact that species names are assigned in the
# same order that species emerge.
# `tree`: phylogenetic tree in phylo format.
# `ancestor_indexes`: a vector with a sequence of ancestor indexes.
# return: a vector with the ids (gen3sis species ids) of the ancestors.
mapAncestorSpecies <- function(tree, ancestor_indexes){

  # Number of ancestor nodes.
  n_nodes <- length(ancestor_indexes)
  # Coerce the ancestor indexes to a numeric format.
  ancestor_indexes <- as.numeric(ancestor_indexes)

  # Create a vector to get the ids (gen3sis species ids) associated to each
  # ancestor.
  ancestor_ids <- rep(NA, n_nodes)

  # Get the id for each ancestor of the sequence.
  for(i in 1:n_nodes){
    # Index (in the phylogenetic tree) of the ancestor (node) i.
    node <- ancestor_indexes[i]
    # Get the list of indexes of the descendants (tips) of the node i.
    descendant_indexes <- phangorn::Descendants(tree, node)[[1]]
    # Get the labels (tip names) of the descendants of node i.
    descendant_names <- tree$tip.label[descendant_indexes]
    # Get the gen3sis ids of the descendants from the tip names.
    descendant_ids <- as.numeric(gsub("species", "", descendant_names))
    # Use the pattern of nomenclature of species in the gen3sis package to get
    # the id of the ancestor node.
    # The lesser id among the descendants corresponds to the id of the founder
    # of the clade, i. e. the id of the ancestor node.
    ancestor_ids[i] <- min(descendant_ids)
  }

  return(ancestor_ids)
}

################################################################################
# listAncestorSpecies
################################################################################
# Function to create a temporal list of ancestor species from a species ancestor
# path (a data.frame within species ancestor id and the time step in which it
# speciate.).
# This list can be used to get information of a lineage along the phylogeny.
# This data can be used in plots, or to calculate metrics for this lineage,
# beyond tracking the evolution of its parameters.
# `species_ancestor_path`: a data.frame with the heights of the nodes of the
# focal lineage and the ids of the nodes (ancestors).
# `Ntimesteps`: number of time-steps of the gen3sis simulation.
# return: a vector with a temporal series of the gen3sis id of the ancestor of
# focal species at each time.
listAncestorSpecies <- function(species_ancestor_path, Ntimesteps){

  # Create a vector to get the id of the ancestors of focal species at each
  # time-step.
  ancestor_species_id <- rep(NA, Ntimesteps)

  # Number of different ids of ancestors.
  k <- length(species_ancestor_path$ancestor_height) + 1
  # Create a vector with the times when occur speciations in the lineage.
  time_threshold <- c(Ntimesteps, species_ancestor_path$ancestor_height)
  # Create a sequence with the ids of the ancestor of the focal species.
  ancestor_ids <- c(species_ancestor_path$ancestor_id[1], species_ancestor_path$ancestor_id)

  # For each time-step get the id of the ancestor of the focal species.
  for(t in 1:Ntimesteps){

    # Set the ancestor id at time-step `t`.
    ancestor_species_id[t] <- ancestor_ids[k]

    # Verify if the time exceed the next `time_threshold`.
    # If TRUE (exceed), then change the id to set in the next time-steps.
    if(t>=time_threshold[k]){
      k <- k-1
    }
  }

  return(ancestor_species_id)
}


# Function to rename tips and nodes of the phylogeny.
# `tree`: a phylogenetic tree in phylo format.
renamePhylo <- function(tree){

  # Number of nodes + tips.
  nTipNodes <- 2*tree$Nnode + 1
  # Create a data.frame to get information of node heights.
  nameIndexGuide <- data.frame(id=rep(NA, nTipNodes),
                               height=rep(NA, nTipNodes))

  # For each node/tip, it gets the id and the height of the node and use these
  # information to make a new name for the nodes.
  for(i in 1:nTipNodes){
    #i=1130
    #
    species_ancestors <- getAncestorPath(tree, node=i)
    ancestor_ids <- mapAncestorSpecies(tree, names(species_ancestors))
    ancestor_ids <- rev(ancestor_ids)

    nAncestors <- length(ancestor_ids)
    new_name <- paste("sp", ancestor_ids[1], sep="_")
    id <- ancestor_ids[1]

    if(nAncestors>1){
      for(id in ancestor_ids[2:nAncestors]){
        new_name <- paste(new_name, id, sep=".")
      }
    }

    rownames(nameIndexGuide)[i] <- new_name
    nameIndexGuide$id[i] <- id
    nameIndexGuide$height[i] <- species_ancestors[1]

  }

  # Rename the tips and the nodes of the phylogeny.
  tree$tip.label <- rownames(nameIndexGuide)[1:(tree$Nnode + 1)]
  tree$node.label <- rownames(nameIndexGuide)[(tree$Nnode + 2):nTipNodes]

  new_tree <- list(tree=tree, nameIndexGuide=nameIndexGuide)

  return(new_tree)
}


################################################################################
# Function to plot species abundance.
plotSpeciesAbundance <- function(metacommunity,
                                 landscape_basis,
                                 species_id=1,
                                 colorscale=c(0, 1, 3, 10, 30, 100, 300),
                                 colors=c("#ffffff", "#ffff80", "#ffff00", "#ff8000", "#ff4000", "#ff0000", "#800000")){


  # Creating a a basic data.frame to species abundance data.
  species_data <- data.frame(
    "long"=landscape_basis$long,
    "lat"=landscape_basis$lat,
    "abundance"=landscape_basis$hab)

  rownames(species_data) <- rownames(landscape_basis)

  # Site names of the sites where the species occur.
  select_site_names <- names(metacommunity[[species_id]]$abundance)
  species_data[select_site_names, "abundance"] <-
    metacommunity[[species_id]]$abundance[select_site_names]

  # Defining color scale and breaks position.
  colorbreaks <- truncatedLog(colorscale)
  colorscale_labels <- as.character(colorscale[-length(colorscale)])

  # Plot.
  species_abundance_plot <- ggplot(data=species_data) +
    geom_raster(mapping=aes(x=long, y=lat, fill=truncatedLog(abundance))) +
    scale_fill_gradientn(
      colors=colors,
      name='Abundance',
      space="Lab", # I don't know what this is for.
      limits=c(min(colorbreaks), max(colorbreaks)),
      values=scales::rescale(colorbreaks),
      breaks=colorbreaks[-length(colorbreaks)],
      labels=colorscale_labels) +
    coord_equal() +
    theme_minimal()

  return(species_abundance_plot)
}


#' Extract a compostion matrix of communities and nodes
#'
#' @param simulation_path A string indicating the path from the root directory
#'     to the simulation results from genesis. Each folder must be separeted
#'     by a forward slash "/"
#' @param save Logical. If false (default), the list with the results won't be
#'     saved in your machine
#' @param output_file A string indicating the path where the output of this
#'     function will be saved. If NULL (default), and if save is TRUE the
#'     output will be saved in the directory root
#'
#' @returns A list with four components:
#'     \itemize{
#'         \item{node_site_matrix}{A matrix containing the sites in the rows
#'             and the nodes of the tree in the columns}
#'         \item{tree}{A phylo object with the species produced in the simulation
#'             process. The nodes were renamed from the original phylogeny}
#'         \item{nameIndexGuide}{A data frame with node information}
#'         \item{coordinates}{A data frame with coordinate information for
#'             the sites in node_site_matrix}
#'     }
#' @export
#'
#' @examples
get_node_site_matrix <-
  function(simulation_path,
           save = FALSE,
           output_file = NULL){
    # Load the gen3sis simulation data.

    # Species data at present time.
    # metacommunity_at_present <- readRDS(file.path("species", "species_t_0.rds"))
    # metacommunity_at_present <- readRDS(here::here(simulation_path, "species", "species_t_0.rds"))
    simulation_path_partial <- here::here(simulation_path, "species")
    timestep0 <- dir(simulation_path_partial)[1]
    metacommunity_at_presence <- readRDS(here::here(simulation_path, "species", timestep0))

    # Simulation summary
    #gen3sis_simulation <- readRDS("sgen3sis.rds")
    gen3sis_simulation <- readRDS(here::here(simulation_path, "sgen3sis.rds"))

    # Get phylogeny.
    # tree <- read.nexus("phy.nex")
    tree <- read.nexus(here::here(simulation_path, "phy.nex"))

    # Basic parameters of the simulation.
    Ntimesteps <- length(gen3sis_simulation$summary$phylo_summary[,"total"]) - 1
    time_initial <- as.numeric(rownames(gen3sis_simulation$summary$phylo_summary)[2])
    time_final <- as.numeric(rownames(gen3sis_simulation$summary$phylo_summary)[Ntimesteps + 1])
    Nsites <- length(gen3sis_simulation$summary$`richness-final`[,1])
    nSpeciesTotal <- gen3sis_simulation$summary$phylo_summary[Ntimesteps+1, "total"]
    nSpeciesAlive <- gen3sis_simulation$summary$phylo_summary[Ntimesteps+1, "alive"]

    #########################################################################
    # Rename nodes and tips of the phylogeny.
    #########################################################################

    # Get the id and the height of the nodes and rename the phylogeny.
    renamedPhylo <- renamePhylo(tree)
    tree_n <- renamedPhylo$tree
    #View(renamedPhylo$nameIndexGuide)
    #plot.phylo(tree_n)

    # For each node/tip get the species ocurrences at specific time.
    # For tips the time is the last timestep and for nodes is the last timestep
    # before the cladogenesis
    # The results are agregated in a list of ocurrences.
    ocurrences <- list()

    renamedPhylo$nameIndexGuide$timestep <- NA

    for(node in 1:(2*nSpeciesTotal - 1)){
      id <- renamedPhylo$nameIndexGuide$id[node]
      if(node<=nSpeciesTotal){
        # Tips
        time <- time_initial + sign(time_final - time_initial)*renamedPhylo$nameIndexGuide$height[node]
      }
      else{
        # Nodes
        time <- time_initial + sign(time_final - time_initial)*(renamedPhylo$nameIndexGuide$height[node] - 1)
      }

      renamedPhylo$nameIndexGuide$timestep[node] <- time

      file_to_load <- paste(paste("species_t", time, sep="_"), "rds", sep=".")

      # metacommunity_at_time <- readRDS(file.path("species", file_to_load))
      metacommunity_at_time <- readRDS(here::here(simulation_path, "species", file_to_load))

      ocurrences_node <- list(as.numeric(names(metacommunity_at_time[[id]]$abundance)))
      names(ocurrences_node) <- rownames(renamedPhylo$nameIndexGuide)[node]

      ocurrences <- c(ocurrences, ocurrences_node)
    }

    # Verify the list with the ocurrences of the nodes.

    # Create a ocurrence matrix of the nodes.
    node_site_matrix <- matrix(0, nrow=Nsites, ncol=(2*nSpeciesTotal - 1))
    rownames(node_site_matrix) <- names(gen3sis_simulation$summary$`richness-final`[,1])
    colnames(node_site_matrix) <- names(ocurrences)

    for(node in names(ocurrences)){
      node_site_matrix[ocurrences[[node]], node] <- 1
    }
    #View(node_site_matrix)

    # Remove the extinct species from phylogeny and from the node_site_matrix .
    tree_n_alive <- drop.fossil(tree_n)
    node_names_alive <- c(tree_n_alive$tip.label, tree_n_alive$node.label)
    node_site_matrix_alive <- node_site_matrix[, node_names_alive]
    #View(node_site_matrix_alive)

    # Site coordinates.
    coordinates <- data.frame(long=gen3sis_simulation$summary$`richness-final`[,1], lat=gen3sis_simulation$summary$`richness-final`[,2])
    rownames(coordinates) <- rownames(node_site_matrix)

    # Save node-site data.
    node_site_data <- list(node_site_matrix = node_site_matrix_alive,
                           tree = tree_n_alive,
                           nameIndexGuide = renamedPhylo$nameIndexGuide,
                           coordinates = coordinates)
    if(save){
      if(is.null(output_file)){
        saveRDS(node_site_data, file="node_site_data.rds")
      }
      else{
        saveRDS(node_site_data, file= here::here(output_file))
      }
    }
    return(node_site_data)
  }

