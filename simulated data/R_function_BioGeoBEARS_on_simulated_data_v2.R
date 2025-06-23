#' Run BioGeoBEARS analysis on multiple simulation outputs.
#'
#' This function takes a directory containing simulation output files (RDS format),
#' runs a BioGeoBEARS analysis on each, and returns key results.
#'
#' @param input_dir The directory containing the RDS files (e.g.,
#'   "gen3sis_output_dispersion_all/").
#' @param output_dir The directory where the results (RDS files) will be saved.
#' @param Nx The number of horizontal partitions in the landscape. Default: 2.
#' @param Ny The number of vertical partitions in the landscape. Default: 2.
#' @param num_cores The number of cores to use for parallel processing in
#'   BioGeoBEARS. Default: 1.
#'
#' @return A list containing the results of the BioGeoBEARS analysis for each
#'   input file. Each element of the list is itself a list containing:
#'   \item{restable_AIC_rellike}{A table summarizing the AIC and relative
#'     likelihood of each model.}
#'   \item{best_model_name}{The name of the model with the lowest AIC.}
#'   \item{node_site_presence_absence}{A matrix with nodes as rows and sites as
#'     columns, indicating presence/absence of each node for each site.}
#'    \item{biogeobears_results}{A list containing the BioGeoBEARS results for each model.}
#'
#' @details The function assumes that the RDS files in `input_dir` contain
#'   simulation outputs with the following structure (as in the original script):
#'   \itemize{
#'     \item{`node_site_matrix`}{A matrix of tip/node occurrences at sites.}
#'     \item{`tree`}{A phylogenetic tree (phylo object).}
#'     \item{`coordinates`}{A data frame of site coordinates.}
#'   }
#'   The function performs a BioGeoBEARS analysis using the DEC, DEC+J,
#'   DIVALIKE, DIVALIKE+J, BAYAREALIKE, and BAYAREALIKE+J models.  It then
#'   selects the best model based on AIC and extracts the ancestral state
#'   reconstructions.
#'
#' @export
run_biogeobears_analysis <- function(input_dir, output_dir, Nx = 2, Ny = 2, num_cores = 1) {
  # Get list of RDS files in the input directory
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  
  # Check if there are any RDS files
  if (length(rds_files) == 0) {
    stop(paste("No RDS files found in the specified input directory:", input_dir))
  }
  
  # Load the save_rase_BioGeoBEARS_data.R script
  source("save_rase_BioGeoBEARS_data.R") # Make sure this file is in the same directory, or specify the correct path
  
  # Initialize a list to store the results from each file
  all_results <- list()
  
  # Loop through each RDS file
  for (rds_file in rds_files) {
    # Get the filename without the extension for naming outputs
    filename_base <- tools::file_path_sans_ext(basename(rds_file))
    cat("\nProcessing file:", rds_file, "\n")
    
    # Load the simulation output from the RDS file
    res_sim_comm <- readRDS(rds_file)
    occ_nodes_spp <- res_sim_comm$node_site_matrix
    phy <- res_sim_comm$tree
    phy$edge.length <- phy$edge.length + 0.001 # Add a small value to edge lengths
    coords <- res_sim_comm$coordinates
    comm <- occ_nodes_spp[, match(phy$tip.label, colnames(occ_nodes_spp))]
    real_node.anc.area <- occ_nodes_spp[, -match(phy$tip.label, colnames(occ_nodes_spp)), drop = FALSE]
    
    # Put simulation data in phylip format
    data_BGB <- save_rase_BioGeoBEARS_data(res_sim_comm$node_site_matrix, res_sim_comm$coordinates, Nx = Nx, Ny = Ny, filename = filename_base)
    
    # Load BioGeoBEARS and dependencies
    library(ape)
    library(optimx)
    library(GenSA)
    library(rexpokit)
    library(cladoRcpp)
    library(parallel) # Use parallel, should work on most systems
    library(BioGeoBEARS)
    
    # Setup BioGeoBEARS run
    geogfn <- np(paste0(filename_base, "_biogeobears.txt"))
    tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
    max_range_size <- 4
    trfn <- write.tree(phy, file = paste0("tree_", filename_base, ".nwk"))
    trfn <- np(paste0("tree_", filename_base, ".nwk"))
    
    # Function to run a BioGeoBEARS model
    run_biogeobears_model <- function(model_name, res_object_name, include_null = TRUE, ...) {
      BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
      BioGeoBEARS_run_object$trfn <- trfn
      BioGeoBEARS_run_object$geogfn <- geogfn
      BioGeoBEARS_run_object$max_range_size <- max_range_size
      BioGeoBEARS_run_object$min_branchlength <- 0.000001
      BioGeoBEARS_run_object$include_null_range <- include_null
      BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
      BioGeoBEARS_run_object$return_condlikes_table <- TRUE
      BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
      BioGeoBEARS_run_object$calc_ancprobs <- TRUE
      BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
      check_BioGeoBEARS_run(BioGeoBEARS_run_object)
      BioGeoBEARS_run_object$num_cores_to_use <- num_cores # Use user-specified number of cores
      
      # Set model-specific parameters
      if (model_name == "DEC+J") {
        dstart <- if (exists("resDEC")) get("resDEC")$outputs@params_table["d", "est"] else 0.1
        estart <- if (exists("resDEC")) get("resDEC")$outputs@params_table["e", "est"] else 0.1
        jstart <- 0.0001
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] <- dstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] <- dstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] <- estart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] <- estart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart
      } else if (model_name %in% c("DIVALIKE", "DIVALIKE+J")) {
        BioGeoBEARS_run_object$include_null_range <- TRUE
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y", "type"] <- "ysv*1/2"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "init"] <- 0.5
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "est"] <- 0.5
        if (model_name == "DIVALIKE+J") {
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
          jstart <- if (exists("resDIVALIKE")) get("resDIVALIKE")$outputs@params_table["j", "est"] else 0.0001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "min"] <- 0.00001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "max"] <- 1.99999
        }
      } else if (model_name %in% c("BAYAREALIKE", "BAYAREALIKE+J")) {
        BioGeoBEARS_run_object$include_null_range <- TRUE
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "est"] <- 0.0
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/1"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y", "type"] <- "1-j"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "init"] <- 0.9999
        BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "est"] <- 0.9999
        if (model_name == "BAYAREALIKE+J") {
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
          jstart <- if (exists("resBAYAREALIKE")) get("resBAYAREALIKE")$outputs@params_table["j", "est"] else 0.0001
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart
          BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "max"] <- 0.99999
        }
      }
      BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
      check_BioGeoBEARS_run(BioGeoBEARS_run_object)
      resfn <- paste0("res_", model_name, "_", filename_base, ".Rdata") # Unique filename
      res <- bears_optim_run(BioGeoBEARS_run_object)
      save(res, file = resfn)
      assign(res_object_name, res, envir = .GlobalEnv) #assign the result to the global environment
      return(res)
    }
    
    # Run BioGeoBEARS models
    resDEC <- run_biogeobears_model("DEC", "resDEC", include_null = FALSE)
    resDECj <- run_biogeobears_model("DEC+J", "resDECj", include_null = FALSE)
    resDIVALIKE <- run_biogeobears_model("DIVALIKE", "resDIVALIKE", include_null = TRUE)
    resDIVALIKEj <- run_biogeobears_model("DIVALIKE+J", "resDIVALIKEj", include_null = TRUE)
    resBAYAREALIKE <- run_biogeobears_model("BAYAREALIKE", "resBAYAREALIKE", include_null = TRUE)
    resBAYAREALIKEj <- run_biogeobears_model("BAYAREALIKE+J", "resBAYAREALIKEj", include_null = TRUE)
    
    # CALCULATE SUMMARY STATISTICS TO COMPARE MODELS
    restable <- NULL
    teststable <- NULL
    
    # Function to perform model comparison and extract results
    compare_models <- function(res1, res2, model_name1, model_name2) {
      LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(res2)
      LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(res1)
      numparams1 <- length(res1$outputs@params_table[, "est"][res1$outputs@params_table[, "type"] == "free"])
      numparams2 <- length(res2$outputs@params_table[, "est"][res2$outputs@params_table[, "type"] == "free"])
      stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
      res2_params <- extract_params_from_BioGeoBEARS_results_object(results_object = res2, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
      res1_params <- extract_params_from_BioGeoBEARS_results_object(results_object = res1, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
      tmp_tests <- conditional_format_table(stats)
      return(list(restable = rbind(res2_params, res1_params), teststable = tmp_tests, aic_stats = stats))
    }
    
    # Compare models
    comparison_results <- list(
      DEC_vs_DECj = compare_models(resDECj, resDEC, "DEC+J", "DEC"),
      DIVALIKE_vs_DIVALIKEj = compare_models(resDIVALIKEj, resDIVALIKE, "DIVALIKE+J", "DIVALIKE"),
      BAYAREALIKE_vs_BAYAREALIKEj = compare_models(resBAYAREALIKEj, resBAYAREALIKE, "BAYAREALIKE+J", "BAYAREALIKE")
    )
    # Assemble results tables
    for (comparison in comparison_results) {
      restable <- rbind(restable, comparison$restable)
      teststable <- rbind(teststable, comparison$teststable)
    }
    
    teststable$alt <- c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
    teststable$null <- c("DEC", "DIVALIKE", "BAYAREALIKE")
    row.names(restable) <- c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
    restable <- put_jcol_after_ecol(restable)
    
    # Calculate AICs
    AICtable <- calc_AIC_column(LnL_vals = restable$LnL, nparam_vals = restable$numparams)
    restable_AIC_rellike <- AkaikeWeights_on_summary_table(restable = cbind(restable, AICtable), colname_to_use = "AIC")
    restable_AIC_rellike <- put_jcol_after_ecol(restable_AIC_rellike)
    
    # Select the best model
    best_model_name <- rownames(restable_AIC_rellike)[which.min(restable_AIC_rellike$AIC)]
    cat("The best model based on AIC is:", best_model_name, "\n")
    
    # Create a list to hold results objects
    results_objects <- list(
      DEC = resDEC,
      `DEC+J` = resDECj,
      DIVALIKE = resDIVALIKE,
      `DIVALIKE+J` = resDIVALIKEj,
      BAYAREALIKE = resBAYAREALIKE,
      `BAYAREALIKE+J` = resBAYAREALIKEj
    )
    best_model_output <- results_objects[[best_model_name]]
    
    # Extract probabilities of states/ranges at each node
    areas <- c("A", "B", "C", "D")
    marginal_probs_at_top <- best_model_output$ML_marginal_prob_each_state_at_branch_top_AT_node
    numeric_states_list <- best_model_output$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
    
    # 1. Identify Unique Regions
    unique_regions <- areas
    num_regions <- length(unique_regions)
    
    # 2. Initialize the presence/absence matrix
    num_nodes <- nrow(marginal_probs_at_top)
    presence_absence_matrix <- matrix(0, nrow = num_nodes, ncol = num_regions,
                                      dimnames = list(rownames(marginal_probs_at_top),
                                                      unique_regions))
    
    # 3. Iterate through each internal node
    for (i in 1:num_nodes) {
      node_probs <- marginal_probs_at_top[i, ]
      
      # 4. Find the maximum probability for the current node
      max_prob <- max(node_probs)
      
      # 5. Identify all states with the maximum probability
      most_probable_state_indices <- which(node_probs == max_prob)
      
      # 6. Iterate through the most probable states and mark the constituent regions as present
      for (state_index in most_probable_state_indices) {
        numeric_state <- numeric_states_list[[state_index]]
        regions_in_state <- areas[numeric_state + 1] # Assuming 0-based indexing
        
        # Mark each region in the most probable state(s) as present (1)
        for (region in regions_in_state) {
          if (region %in% colnames(presence_absence_matrix)) {
            presence_absence_matrix[i, region] <- 1
          }
        }
      }
    }
    
    # Create presence/absence for sites instead of areas
    node_site_matrix <- res_sim_comm$node_site_matrix
    coordinates <- res_sim_comm$coordinates
    Nsites <- nrow(coordinates)
    Ncol <- sqrt(Nsites)
    Nbioregions <- Nx * Ny
    classification <- matrix(1, nrow = Nsites, ncol = 3)
    colnames(classification) <- c("x", "y", "bioregion")
    rownames(classification) <- rownames(coordinates)
    for (i in 1:Nsites) {
      classification[i, "x"] <- ((i - 1) %% Ncol) %/% (Ncol / Nx)
      classification[i, "y"] <- ((i - 1) %/% Ncol) %/% (Ncol / Ny)
      classification[i, "bioregion"] <- Nx * (classification[i, "y"]) + classification[i, "x"] + 1
    }
    site_bioregion_matrix <- matrix(0, nrow = Nsites, ncol = Nbioregions)
    rownames(site_bioregion_matrix) <- rownames(classification)
    colnames(site_bioregion_matrix) <- as.character(1:Nbioregions)
    for (site in 1:Nsites) {
      site_bioregion_matrix[site, classification[site, "bioregion"]] <- 1
    }
    
    # Create a mapping from numeric bioregion to letter
    area_mapping <- LETTERS[1:4]
    names(area_mapping) <- as.character(1:4) # Name the mapping for easy lookup
    
    num_nodes <- nrow(presence_absence_matrix)
    num_sites <- nrow(site_bioregion_matrix)
    node_site_presence_absence <- matrix(0, nrow = num_nodes, ncol = num_sites,
                                         dimnames = list(rownames(presence_absence_matrix),
                                                         rownames(site_bioregion_matrix)))
    
    # Iterate through each internal node (row in presence_absence_matrix)
    for (i in 1:num_nodes) {
      # Get the presence/absence for each area for the current node
      area_presence <- presence_absence_matrix[i, ]
      
      # Iterate through each site (row in site_bioregion_matrix)
      for (j in 1:num_sites) {
        # Determine the bioregion (area number) for the current site
        site_bioregion <- which(site_bioregion_matrix[j, ] == 1)
        
        # Map the numeric bioregion to a letter
        site_area_letter <- area_mapping[as.character(site_bioregion)]
        
        # Find the presence/absence value for that area for the current node
        if (site_area_letter %in% colnames(presence_absence_matrix)) {
          node_site_presence_absence[i, j] <- area_presence[site_area_letter]
        } else {
          # Handle potential issues
          warning(paste("Area", site_area_letter, "not found in BioGeoBEARS results."))
        }
      }
    }
    # Store the three outputs.
    all_results[[filename_base]] <- list(
      restable_AIC_rellike = restable_AIC_rellike,
      best_model_name = best_model_name,
      node_site_presence_absence = node_site_presence_absence,
      biogeobears_results = list(resDEC = resDEC, resDECj = get("resDECj"), resDIVALIKE = resDIVALIKE, resDIVALIKEj = resDIVALIKEj, resBAYAREALIKE = resBAYAREALIKE, resBAYAREALIKEj = resBAYAREALIKEj) # Store all results
    )
  }
  return(all_results)
}
