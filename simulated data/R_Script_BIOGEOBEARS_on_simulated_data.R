# Run biogeobears on simulated data:
dir= #define directory path to data
input_directory <- dir
output_directory <- "./BioGeoBEARS_Results/"
# Create the output directory if it doesn't exist
if (!file.exists(output_directory)) {
  dir.create(output_directory)
}

source("save_rase_BioGeoBEARS_data.R")# dependance
source("R_function_BioGeoBEARS_on_simulated_data_v2.R")


results <- run_biogeobears_analysis(input_dir = input_directory, output_dir = output_directory, Nx = 2, Ny = 2, num_cores = 8)

# Save the results to an RDS file in the output directory
saveRDS(results, file = file.path(output_directory, "all_biogeobears_results.rds"))

# The most important outputs are here:
results$output_dispersion1$node_site_presence_absence

# for each simulation, $output_dispersionX$node_site_presence_absence
# where X is the simulation number
# and node_site_presence_absence have nodes as rows and sites as columns
# where each site has the presence/absence value of its corresponding area

# Taking only the lineages_by_site matrices:
all_presence_absence_matrices <- list()
# Iterate through the 'results' list
for (name in names(results)) {
  # Check if the element starts with 'output_dispersion' and has the target matrix
  if (grepl("^output_dispersion\\d+$", name) &&
      !is.null(results[[name]]$node_site_presence_absence) &&
      is.matrix(results[[name]]$node_site_presence_absence)) {

    # Add the matrix to our list
    all_presence_absence_matrices[[name]] <- results[[name]]$node_site_presence_absence
  }
}

# Print the resulting list of lineages_by_site matrices
print(all_presence_absence_matrices)

# NOTE: The final results have live species + Nodes
# The nodes start at the end of the live species
