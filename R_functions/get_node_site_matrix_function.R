# Packages.
library(tidyr)
library(ape)
library(gen3sis)
source("~/Dropbox/ada_simulations/mygen3sisfunctions/phylo_gen3sis_tools_1.R")

# Path names.
simulation_main <- "~/Dropbox/ada_simulations/simulation_output/"
sim <- "config_neutral_zero_sum" # Name of simulation dir.
simulation_path <- file.path(simulation_main, sim)

setwd(simulation_path)

# Set the output path. 
output_file <- file.path(simulation_main, "simulation_data", sim)

node_site_data <- get_node_site_matrix(simulation_path, sim)

# Check
node_site_data$nameIndexGuide
plot.phylo(node_site_data$tree, no.margin=T, show.node.label=T, root.edge=T)
raster::plot(rasterFromXYZ(data.frame(node_site_data$coordinates, z=node_site_data$node_site_matrix[,"sp_1.2"])))

################################################################################
# `get_node_site_data`:
# Function do built a community matrix (sites X species) for gen3sis simulation 
# data with the occurrences of species and nodes. 
# For tips/species, occurrences are get from their last timestep alive in simulation. 
# For nodes, the occurrences are get from their last timestep before the cladogenesis.
# Inputs:
# `simulation_path`: a string with the path where are the simulation results.
# `save`: a boolean value. If TRUE the results are saved in a file. Default is 
# FALSE.
# `output_file`: a string containing the file name/path to save the results if 
# save=T. If nothing, the results are saved as "node_site_data.rds" in the 
# current directory.
# And it has as output a list containing:
# `node_site_matrix`: a matrix with data of occurrence of each species/tip 
# and node in the sites of the landscape.
# `tree`: the phylogenetic tree of the simulation with tips and nodes renamed.
# `nameIndexGuide`: a data.frame with a index with the new names, the old labels
# (ids) of species and nodes (when species in the simulation), the height of the 
# tip/nodes in the phylogeny (distance to the root) and the time-step of which 
# the occurrence of tip/node was got.
# `coordinates`: a data.frame with the coordinates of the sites.
get_node_site_matrix <- function(simulation_path, save=F, output_file=NULL){
  # Load the gen3sis simulation data.
  # Species data at present time.
  metacommunity_at_present <- readRDS(file.path("species", "species_t_0.rds"))
  
  # Simulation resume.
  gen3sis_simulation <- readRDS("sgen3sis.rds")
  # Get phylogeny.
  tree <- read.nexus("phy.nex")
  
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
    
    metacommunity_at_time <- readRDS(file.path("species", file_to_load))
    
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
  node_site_data <- list(node_site_matrix=node_site_matrix_alive, tree=tree_n_alive, nameIndexGuide=renamedPhylo$nameIndexGuide, coordinates=coordinates)
  if(save){
    if(is.null(output_file)){
      saveRDS(node_site_data, file="node_site_data.rds")
    }
    else{
      saveRDS(node_site_data, file=output_file)
    }  
  }
  return(node_site_data)
}
