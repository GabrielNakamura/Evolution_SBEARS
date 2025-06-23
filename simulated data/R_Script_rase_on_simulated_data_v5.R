devtools::install_github("ignacioq/rase")
# Load the functions
source("R_function_rase_on_simulated_data_v5.R")
source("save_rase_BioGeoBEARS_data.R")

getwd()
# Input directory containing your simulation files
input_dir <- "/Users/leandroduarte/Desktop/ADA/teste_metodos/Matrizes simuladas/gen3sis_output_lowdiv_100/data"

# Output directory for results (will be created if doesn't exist)
output_dir <- "rase_results"

# File pattern to identify simulation files
file_pattern <- "output_dispersion.*\\.rds"

# RASE parameters (adjust as needed)
rase_params <- list(
  niter = 12000,
  sigma2_scale = 0.05,
  params0 = NA,
  nGQ = 20,
  burnin = 20,
  min_coords = 5,
  contour_level = 0.8,
  save_intermediate = TRUE
)

# Run the analysis
# Get list of simulation files
sim_files <- list.files(input_dir,
                        pattern = file_pattern,
                        full.names = TRUE)

# Check if any files were found
if (length(sim_files) == 0) {
  stop("No simulation files found matching the pattern: ", file_pattern)
}

message("Found ", length(sim_files), " simulation files to process")

# Process all files (using do.call to pass parameters as list)
all_results <- lapply(sim_files, function(f) {
  message("\nProcessing file: ", basename(f))

  args <- c(list(sim_file = f, output_dir = output_dir), rase_params)

  tryCatch({
    do.call(run_rase_analysis, args)
  }, error = function(e) {
    message("Failed to process ", basename(f), ": ", e$message)
    return(NULL)
  })
})

# Remove NULL results from failed runs
all_results <- Filter(Negate(is.null), all_results)

# Save and report results
if (length(all_results) > 0) {
  # Save combined results
  saveRDS(all_results, file.path(output_dir, "all_rase_comparison_results.rds"))

  # Summary report
  message("\nAnalysis complete!")
  message("Successfully processed ", length(all_results), " out of ", length(sim_files), " files")
  message("Results saved in: ", normalizePath(output_dir))

  # Print summary of first result for verification
  message("\nStructure of first result:")
  str(all_results[[1]], max.level = 2)
} else {
  warning("No files were successfully processed")
}

#### all_results contains real and rase estimated ancestral locations
#### in a grid, site-by-site basis
#### each element of the list contains real and estimated anc for the same simulation
all_results
