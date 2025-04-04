################################################################################
# workshop_2.R
# Date: December 14, 2023.
# Author: Marcelo Henrique Schwade
# Based on: workshop_1.R
#
# This script built a node-species community matrix using species presence data 
# at each particular timestep (final to tips and the first timestep before 
# speciation). The script also save data of species ocurrences for `rase` and 
# `BioGeoBEARS` (FINISHED).
# This script is an automated and very improved version of the workshop_1.R.

################################################################################

# Packages.
library(tidyr)
library(ape)
source("~/Desktop/simple_simulation_20231130/phylo_gen3sis_tools_1.R")

# Path names.
simulation_main <- "/Users/marceloschwade/Desktop/simple_simulation_20231130/simulation_outputs5"
simulation_path <- file.path(simulation_main, "default_config")

simulation_list <- read.table(file.path(simulation_main, "list_simulation_outputs.txt"), header=F)[,1]

################################################################################
# For each sim in simulation_list:
# 1) Get species ocurrences at present time.
# 2) Get phylogeny (final/present).
# 3) Create a list with the extinct species.
# 4) Remove the extinct species of the phylogeny.
# 5) Get ocurrences of nodes (ancestor species) in the past time.
# 6) Built a community matrix with the species ocurrences of tips (present time)
#    and nodes (different past times).
# 7) Classify sites in biorregions (large square site agregations).
# 8) Classify species ocurrences in biorregions, save data in phylip format to 
#    `BioGeoBEARS`.
# 9) Save the coordinates of species ocurrences to `rase`.

for(sim in simulation_list){

  ##############################################################################
  # Built Community Matrix at present.
  cat("Start processing of a new simulation\n")

  # Set the path with the simulation results.
  setwd(file.path(simulation_path, sim))

  # Set the output path. 
  output_path <- file.path(simulation_main, "simulation_data", sim)
  dir.create(output_path, recursive=T)

  # Load the gen3sis simulation data.
  # Species data at present time.
  metacommunity_at_present <- readRDS(file.path("species", "species_t_0.rds"))
  # Simulation resume.
  gen3sis_simulation <- readRDS("sgen3sis.rds")
  # Get phylogeny.
  tree <- read.nexus("phy.nex")
  
  # Basic parameters of the simulation.
  Ntimesteps <- length(gen3sis_simulation$summary$phylo_summary[,"total"]) - 2
  Nsites <- length(gen3sis_simulation$summary$`richness-final`[,1])
  nSpeciesTotal <- gen3sis_simulation$summary$phylo_summary[Ntimesteps+2, "total"]
  nSpeciesAlive <- gen3sis_simulation$summary$phylo_summary[Ntimesteps+2, "alive"]

  cat("Basic parameters imported.\n")

  # This suits only to verify the richness compatibility.
  ##############################################################################

  species_presence <- list()
  
  for(species_id in 1:nSpeciesTotal){
    if(length(metacommunity_at_present[[species_id]]$abundance)>0){
      species_id_presence <- list(id=species_id, 
                                  sites_with_presence=as.numeric(names(metacommunity_at_present[[species_id]]$abundance)))
      species_presence <- c(species_presence, list(species_id_presence))
      names(species_presence)[length(species_presence)] <- as.character(species_id)
    }
  }
  
  # Verify the names of the list (names or the ids of species)
  # names(species_presence)
  cat("Species presences obtained.\n")
  
  # Create the community matrix.
  community_matrix <- matrix(0, nrow=Nsites, ncol=nSpeciesAlive)
  colnames(community_matrix) <- names(species_presence)

  cat("Community matrix created\n")
  
  for(species_id in names(species_presence)){
    community_matrix[species_presence[[species_id]]$sites_with_presence, species_id] <- 1
  }
  
  cat("Community matrix built.\n")

  # Verify if the richness in each site matches with the final richness of the simulation.
  richness_final <- replace_na(gen3sis_simulation$summary$`richness-final`[,3], 0)
  richness <- rowSums(community_matrix)
  if(all(richness==richness_final)){
    print("Richness verification is OK!")
  }
  else{
    print("An incompatibility was found in the richness verification. The process was aborted.")
    break
  }

  # Viewing phylogeny.
  ##############################################################################

  #par(mfrow=c(1,1))
  #plot.phylo(tree, show.tip.label=F, root.edge=T)

  # Rename nodes and tips of the phylogeny.
  ##############################################################################

  # Get the id and the height of the nodes and rename the phylogeny.
  renamedPhylo <- renamePhylo(tree)
  tree_n <- renamedPhylo$tree
  cat("Phylogeny renamed\n")

  # Get the ocurrences of ancestor nodes.
  ##############################################################################

  # For each node/tip get the species ocurrences at specific time.
  # The results are agregated in a list of ocurrences.
  ocurrences <- list()

  for(node in 1:(2*nSpeciesTotal - 1)){
  
    id <- renamedPhylo$nameIndexGuide$id[node]
    time <- Ntimesteps - renamedPhylo$nameIndexGuide$height[node]
  
    file_to_load <- paste(paste("species_t", time, sep="_"), "rds", sep=".")
  
    metacommunity_at_time <- readRDS(file.path("species", file_to_load))
  
    ocurrences_node <- list(as.numeric(names(metacommunity_at_time[[id]]$abundance)))
    names(ocurrences_node) <- rownames(renamedPhylo$nameIndexGuide)[node]
  
    ocurrences <- c(ocurrences, ocurrences_node)
  
  }

  cat("Node ocurrences obtained\n")

  # Verify the list with the ocurrences of the nodes.

  # Create a ocurrence matrix of the nodes.
  node_site_matrix <- matrix(0, nrow=Nsites, ncol=(2*nSpeciesTotal - 1))
  rownames(node_site_matrix) <- names(gen3sis_simulation$summary$`richness-final`[,1])
  colnames(node_site_matrix) <- names(ocurrences)

  for(node in names(ocurrences)){
    node_site_matrix[ocurrences[[node]], node] <- 1
  }

  # Remove the extinct species from phylogeny and from the node_site_matrix .
  tree_n_alive <- drop.fossil(tree_n)
  node_names_alive <- c(tree_n_alive$tip.label, tree_n_alive$node.label)
  node_site_matrix_alive <- node_site_matrix[, node_names_alive]

  # Site coordinates.
  coordinates <- data.frame(long=gen3sis_simulation$summary$`richness-final`[,1], lat=gen3sis_simulation$summary$`richness-final`[,2])
  rownames(coordinates) <- rownames(node_site_matrix)

  cat("Node - site data obtained.\n")

  # Save node-site data.
  node_site_data <- list(node_site_matrix=node_site_matrix_alive, nameIndexGuide=renamedPhylo$nameIndexGuide, coordinates=coordinates)
  file_to_save <- paste("node_site_data", "rds", sep=".")
  saveRDS(node_site_data, file=file.path(output_path, file_to_save))
  
  cat("Node - site data saved.\n")

  # Save phylogenies.
  write.nexus(tree_n, file=file.path(output_path, "phy_total.nex"))
  write.nexus(tree_n_alive, file=file.path(output_path, "phy_alive.nex"))

  cat("Phylogenies saved.\n")

  # Create data to BioGeoBEARS and rase.
  ##############################################################################
  
  # Number of bioregions.
  Nbioregions <- 9
  longmin <- min(coordinates$long)
  longmax <- max(coordinates$long)
  latmin <- min(coordinates$lat)
  latmax <- max(coordinates$lat)
  
  # Thresholds.
  bioregion_limits <- data.frame(lat=rep(NA, 4), long=rep(NA, 4))
  bioregion_limits$lat <- seq(latmin, latmax, by=(latmax-latmin)/3)
  bioregion_limits$long <- seq(longmin, longmax, by=(longmax-longmin)/3)
  
  # Classify the sites between bioregions.
  classification <- matrix(1, nrow=Nsites, ncol=3)
  colnames(classification) <- c("long", "lat", "bioregion")
  rownames(classification) <- rownames(coordinates)
  for(i in 1:Nsites){
    # Longitudinal classification.
    if(coordinates[i, "long"]>bioregion_limits$long[2]){
      classification[i, "long"] <- 2
      #print("2")
    }
    if(coordinates[i, "long"]>bioregion_limits$long[3]){
      classification[i, "long"] <- 3
      #print("3")
    }
    
    # Latitudinal classification.
    if(coordinates[i, "lat"]<bioregion_limits$lat[3]){
      classification[i, "lat"] <- 2
      #print("2")
    }
    if(coordinates[i, "lat"]<bioregion_limits$lat[2]){
      classification[i, "lat"] <- 3
      #print("3")
    }
    
    classification[i, "bioregion"] <- 3*(classification[i, "lat"] - 1) + classification[i, "long"]
  }

  #View(data.frame(coordinates, classification))
  
  # Create a matrix site x bioregion
  site_bioregion_matrix <- matrix(0, nrow=Nsites, ncol=Nbioregions)
  rownames(site_bioregion_matrix) <- rownames(classification)
  colnames(site_bioregion_matrix) <- as.character(1:Nbioregions)
  for(site in 1:Nsites){
    site_bioregion_matrix[site, classification[site, "bioregion"]] <- 1
  }
  
  #View(site_bioregion_matrix)
  #View(community_matrix) # Site X Species
  #dim(community_matrix)
  #dim(site_bioregion_matrix)
  
  species_bioregion_matrix <- t(community_matrix)%*%site_bioregion_matrix
  #dim(species_bioregion_matrix)
  #View(species_bioregion_matrix)
  
  species_incidence_bioregion_matrix <- ifelse(species_bioregion_matrix>0, yes=1, no=0)
  
  
  # Sites of coordinates of species 
  species_coordinates <- list()
  for(i in 1:nSpeciesAlive){
    species_i <- species_presence[[i]]
    species_coordinates_i <- coordinates[species_i$sites_with_presence,]
    species_coordinates <- c(species_coordinates, list(species_coordinates_i))
  }
  
  #
  data_for_rase <- species_coordinates
  
  #species_incidence_bioregion_matrix
  text_line <- paste(rownames(species_incidence_bioregion_matrix), " ", sep="")
  
  for(i in 1:Nbioregions){
    text_line <- paste(text_line, species_incidence_bioregion_matrix[,i], sep="")  
  }
  
  #text_line
  text_head <- paste(paste(nSpeciesAlive, Nbioregions, sep=" "), "(A B C D E F G H I)", sep=" ")
  
  # Write a biogeobears file in phylip format.
  biogeobears_datafile <- file.path(output_path, "data_biogeobears.txt")
  
  sink(biogeobears_datafile)
  cat(text_head, "\n")
  for(i in 1:nSpeciesAlive){
    cat(text_line[i], "\n")
  }
  sink()
  
  # Saving the data of species coordinates for rase in .rds files.
  rase_datafile <- file.path(output_path, "data_rase.rds")
  saveRDS(data_for_rase, file=rase_datafile)
  
}

