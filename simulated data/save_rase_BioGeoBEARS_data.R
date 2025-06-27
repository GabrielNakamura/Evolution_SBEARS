#' Format data for BioGeoBEARS and rase analyses
#'
#' @description
#' This function divide the landscape in Nx*Ny rectangular bioregions, with Nx
#' horizontal partitions and Ny vertical partitions. Then, the sites of the
#' landscape are classified in the bioregions. Then, it compute the presence of
#' the species in the bioregions. This data are saved in phylip format to use
#' with the BioGeoBEARs package ("data_biogeobears.txt").
#' This function also creates a list with the coordinates of the sites where each
#' species is present for use with rase package ("data_rase.rds"). This object is
#' also returned when the function is used.
#'
#' @param node_site_matrix matrix/array. a matrix or data.frame with the tips(+nodes)
#'      presences and absences. Tips(+nodes) in columns, sites in lines.
#' @param coordinates data.frame. an object with the coordinates (lat and long) of
#'      all the sites.
#' @param Nx integer. The number of horizontal partitions in the landscape to create
#'      the bioregions. Default: 2.
#' @param Ny integer. The number of vertical partitions in the landscape to create
#'      the bioregions. Default: 2.
#' @param output_path string. The path to save the data. Default: "./" (local path).
#' @param filename string. A character vector specifying the name of the files of
#     BioGeoBEARS and rase data. Default: "data".
#'
#' @returns
#' @export
#'
#' @examples
save_rase_BioGeoBEARS_data <-
  function(node_site_matrix, coordinates, Nx=2, Ny=2, output_path="./", filename="data"){

  # Get basic values.
  Nsites <- nrow(coordinates)
  Ncol <- sqrt(Nsites)
  nSpeciesAlive <- (ncol(node_site_matrix) + 1)/2

  # Number of bioregions.
  Nbioregions <- Nx*Ny
  longmin <- min(coordinates$long)
  longmax <- max(coordinates$long)
  latmin <- min(coordinates$lat)
  latmax <- max(coordinates$lat)

  # Thresholds.
  #dx <- (longmax-longmin)/Nx
  #dy <- (latmax-latmin)/Ny
  #bioregion_limits <- list(lat=rep(NA, Ny+1), long=rep(NA, Nx+1))
  #bioregion_limits$lat <- seq(latmin, latmax, by=dy)
  #bioregion_limits$long <- seq(longmin, longmax, by=dx)

  # Classify the sites in bioregions.
  classification <- matrix(1, nrow=Nsites, ncol=3)
  colnames(classification) <- c("x", "y", "bioregion")
  rownames(classification) <- rownames(coordinates)
  for(i in 1:Nsites){
    # In this classification, I used the index's sequence of site coordinates to
    # get the horizontal and the vertical positions of sites in terms of
    # bioregions and then obtain the bioregion of each site.
    # I tried making the classification using the coordinates directly, but
    # until now I wasn't successful.
    classification[i, "x"] <- ((i-1)%%Ncol)%/%(Ncol/Nx)
    classification[i, "y"] <- ((i-1)%/%Ncol)%/%(Ncol/Ny)

    classification[i, "bioregion"] <- Nx*(classification[i, "y"]) + classification[i, "x"] + 1
  }
  # For test and view bioregions.
  #bioras <- raster::rasterFromXYZ(data.frame(coordinates, classification[,3]))
  #raster::plot(bioras)
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

  node_bioregion_matrix <- t(node_site_matrix)%*%site_bioregion_matrix
  #dim(species_bioregion_matrix)
  #View(species_bioregion_matrix)

  node_bioregion_matrix <- ifelse(node_bioregion_matrix>0, yes=1, no=0)

  # Sites of coordinates of species
  species_coordinates <- list()
  for(i in 1:nSpeciesAlive){
    sites_species_i <- names(node_site_matrix[node_site_matrix[, i]==1, i])#species_presence[[i]]
    species_coordinates[[colnames(node_site_matrix)[i]]] <- coordinates[sites_species_i,]
  }

  # Data for create polygons for rase.
  data_for_rase <- species_coordinates

  # Write text to save the data in phylip format for BioGeoBEARS.
  text_line <- paste(rownames(node_bioregion_matrix), " ", sep="")

  for(i in 1:Nbioregions){
    text_line <- paste(text_line, node_bioregion_matrix[,i], sep="")
  }

  #text_line
  text_head <- paste(paste(nSpeciesAlive, Nbioregions, sep=" "), paste0("(", paste(LETTERS[1:Nbioregions], collapse=" "), ")"), sep=" ")

  # Write a biogeobears file in phylip format.
  biogeobears_datafile <- file.path(output_path, paste0(filename, "_biogeobears.txt"))

  sink(biogeobears_datafile)
  cat(text_head, "\n")
  for(i in 1:nSpeciesAlive){
    cat(text_line[i], "\n")
  }
  sink()

  # Saving the data of species coordinates for rase in .rds files.
  rase_datafile <- file.path(output_path, paste0(filename, "_rase.rds"))
  saveRDS(data_for_rase, file=rase_datafile)

  return(data_for_rase)
}
