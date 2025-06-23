#' Run RASE analysis on simulated biogeographic data
#'
#' This function performs RASE (Range Ancestral State Estimation) analysis on simulated biogeographic data.
#' It processes the input data, runs the RASE algorithm, and compares the estimated ancestral ranges
#' with the "true" ancestral ranges from the simulation.
#'
#' @param sim_file Path to the simulation RDS file (e.g., "output_dispersionX.rds")
#' @param output_dir Directory to save results (default = "rase_results")
#' @param niter Number of MCMC iterations for RASE (default = 120)
#' @param sigma2_scale Scale parameter for RASE (default = 0.01)
#' @param params0 Initial parameters for RASE (default = NA)
#' @param nGQ Number of quadrature points for RASE (default = 20)
#' @param burnin Percentage of iterations to discard as burn-in (default = 20)
#' @param Nx Number of grid cells in x-direction for range maps (default = 2)
#' @param Ny Number of grid cells in y-direction for range maps (default = 2)
#' @param min_coords Minimum number of coordinates required for a species to be included (default = 5)
#' @param contour_level Contour level for KDE estimation (default = 0.8)
#' @param save_intermediate Should intermediate files be saved? (default = TRUE)
#'
#' @return A list containing RASE results and comparison data
#' @export
#' 
#' # Load required packages (order is important)
library(ape)
library(sf)
library(dplyr)
library(rase)
library(polyCub)
library(spatstat)
library(spatstat.geom)
library(raster)
library(sp)
library(prettyGraphs)
library(ks)
library(letsR)
library(terra)

run_rase_analysis <- function(sim_file, 
                              output_dir = "rase_results", 
                              niter = 120, 
                              sigma2_scale = 0.01, 
                              params0 = NA, 
                              nGQ = 20, 
                              burnin = 20, 
                              Nx = 2, 
                              Ny = 2, 
                              min_coords = 5, 
                              contour_level = 0.8,
                              save_intermediate = TRUE) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  tryCatch({
    ##### Load simulation ####
    message(paste("Processing:", sim_file))
    res_sim_comm <- readRDS(sim_file)
    occ_nodes_spp <- res_sim_comm$node_site_matrix
    phy <- res_sim_comm$tree
    phy$edge.length <- phy$edge.length + 0.001
    coords <- res_sim_comm$coordinates
    comm <- occ_nodes_spp[, match(phy$tip.label, colnames(occ_nodes_spp))]
    real_node.anc.area <- occ_nodes_spp[, -match(phy$tip.label, colnames(occ_nodes_spp))]
    
    ##### Transform simulation into range maps to run rase ####
    source("save_rase_BioGeoBEARS_data.R")
    data_rase <- save_rase_BioGeoBEARS_data(res_sim_comm$node_site_matrix, 
                                            res_sim_comm$coordinates, 
                                            Nx = Nx, Ny = Ny, 
                                            filename = file.path(output_dir, "dados"))
    
    if (save_intermediate) {
      saveRDS(data_rase, file.path(output_dir, "dados_rase.rds"))
    }
    
    # Create polygons for each live species
    spcoords <- data_rase[sapply(data_rase, nrow) >= min_coords]
    
    # Iterate transformation for all
    polys <- vector(mode = "list", length(spcoords))
    for (i in 1:length(spcoords)) {
      xys <- st_as_sf(spcoords[[i]], coords = c("long", "lat"))
      polys[[i]] <- xys %>% 
        dplyr::summarise() %>%
        st_cast("POLYGON") %>%
        st_convex_hull()
    }
    
    ##### Setup for rase ####
    # Convert sf polygons to SpatialPolygonsDataFrame to use rase
    polys_spat <- vector(mode = "list", length(polys))
    for (i in 1:length(polys)) {
      polys_spat[[i]] <- as_Spatial(polys[[i]])
    }
    
    # Convert list of SpatialPolygons to single one
    z <- lapply(polys_spat, function(i) {
      SpatialPolygonsDataFrame(i, data.frame(id = 1:length(i)), match.ID = FALSE)
    })
    sp <- bind(z)
    spdf <- SpatialPolygonsDataFrame(sp, data.frame(id = 1:length(sp))) 
    
    # Transform into owin shapes
    owin.shapes <- shape.to.rase(spdf)
    y <- as(spdf, "SpatialPolygons")
    p <- slot(y, "polygons")
    v <- lapply(p, function(z) { SpatialPolygons(list(z)) })
    winlist <- lapply(v, as.owin)
    
    # Remove species with insufficient coordinates
    tree <- keep.tip(phy, which(sapply(spcoords, nrow) >= min_coords))
    
    # Attribute names to polys
    winlist <- name.poly(owin.shapes, tree, poly.names = tree$tip.label)
    
    ##### Run rase ####
    area <- spatstat.geom::area
    is.empty <- spatstat.geom::is.empty
    rase_results <- rase(tree, winlist, niter = niter, sigma2_scale = sigma2_scale,
                         params0 = params0, nGQ = nGQ)
    
    ##### Process rase results ####
    burnin_samples <- floor(niter * (burnin/100))
    rase_results_subset <- rase_results[(burnin_samples + 1):nrow(rase_results), ]
    
    # Create KDE contours for each node
    hpi_all <- vector(mode = "list", (ncol(rase_results_subset)/2 - 1))
    for (i in 1:(ncol(rase_results_subset)/2 - 1)) {
      df <- data.frame(rase_results_subset[, i], 
                       rase_results_subset[, i + (ncol(rase_results)/2 - 1)])
      hh <- Hpi(df, binned = TRUE) * 1
      dd <- kde(df, H = hh)
      
      # Generate contour lines
      num_levels_requested <- 20
      cc <- contourLines(x = dd$eval.points[[1]], y = dd$eval.points[[2]],
                         z = dd$estimate, levels = seq(0, max(dd$estimate), 
                                                       length = num_levels_requested))
      
      if (length(cc) > 0) {
        target_index <- round(contour_level * length(cc))
        target_index <- max(1, min(target_index, length(cc)))
        hpi_all[[i]] <- cbind(cc[[target_index]]$x, cc[[target_index]]$y)
      }
    }
    
    ##### Create node polygons ####
    nodepolys <- vector(mode = "list", length(hpi_all))
    for (i in 1:length(hpi_all)) {
      if (!is.null(hpi_all[[i]])) {
        colnames(hpi_all[[i]]) <- c("long", "lat")
        nodexys <- st_as_sf(as.data.frame(hpi_all[[i]]), coords = c("long", "lat"))
        nodepolys[[i]] <- nodexys %>% 
          dplyr::summarise() %>%
          st_cast("POLYGON") %>%
          st_convex_hull()
      }
    }
    
    # Convert to SpatialPolygonsDataFrame
    nodepolys_spat <- vector(mode = "list", length(nodepolys))
    for (i in 1:length(nodepolys)) {
      if (!is.null(nodepolys[[i]])) {
        nodepolys_spat[[i]] <- as_Spatial(nodepolys[[i]])
      }
    }
    
    z <- lapply(nodepolys_spat, function(i) {
      if (!is.null(i)) {
        SpatialPolygonsDataFrame(i, data.frame(id = 1:length(i)), match.ID = FALSE)
      }
    })
    z <- Filter(Negate(is.null), z)
    sp <- bind(z)
    nodespdf <- SpatialPolygonsDataFrame(sp, data.frame(id = 1:length(sp))) 
    names(nodespdf) <- "sciname"
    
    ##### Process true ancestral areas ####
    nodes_to_keep <- tree$node.label
    real_node.anc.area <- real_node.anc.area[, colnames(real_node.anc.area) %in% nodes_to_keep]
    
    ##### Fit node ranges into existing grid cells ####
    st_grid <- st_as_sf(coords, coords = c("long", "lat"))
    st_grid3 <- terra::vect(st_grid)
    st_grid3$ID <- 1:length(st_grid3)
    
    teste <- lets.presab.grid(nodespdf, st_grid3, sample.unit = "ID", remove.sp = FALSE)
    rasePAMnodes <- teste$PAM[, -1] # remove ID column
    colnames(rasePAMnodes) <- colnames(real_node.anc.area)
    rase_node.anc.area <- rasePAMnodes
    
    ##### Save and return results ####
    result <- list(
      rase_results = rase_results,
      tree = tree,
      real_node.anc.area = real_node.anc.area,
      rase_node.anc.area = rase_node.anc.area
    )
    
    base_name <- tools::file_path_sans_ext(basename(sim_file))
    save_path <- file.path(output_dir, paste0("rase_", base_name, "_results.rds"))
    saveRDS(result, save_path)
    
    return(list(
      real_node.anc.area = real_node.anc.area,
      rase_node.anc.area = rase_node.anc.area,
      tree = tree
    ))
    
  }, error = function(e) {
    message(paste("Error processing file", sim_file, ":", e$message))
    return(NULL)
  })
}
