<!-- README.md is generated from README.Rmd. Please edit that file -->

# SBEARS: Data & Code Repository

This repository contains all data and R scripts required to replicate
the analyses in the manuscript:

**SBEARS – A site-based method to estimate ancestral ranges of species**
(in review).

------------------------------------------------------------------------

## Repository Structure

To view the main repository structure, run:

``` r
library(fs)
tracked_files <- system("git ls-files", intern = TRUE)
top_level_dirs <- unique(dirname(tracked_files))
top_level_dirs <- top_level_dirs[top_level_dirs != "."]
all_top_dirs <- dir_ls(path = ".", type = "directory", recurse = 0)
tracked_top_dirs <- all_top_dirs[basename(all_top_dirs) %in% top_level_dirs]
purrr::walk(tracked_top_dirs, ~ dir_tree(.x, recurse = 0))
```

------------------------------------------------------------------------

## Main Folders

- **R_functions**  
  Auxiliary functions used in scripts and RMarkdown files for analyses.

- **R_scripts**  
  All RMarkdown files needed to run the simulations described in the
  paper.

- **gen3sis_simulations**  
  Files and folders for simulation parameter settings:

  - `config`:
    - `config_presence_high_diversity`: Parameters for high species
      richness scenario (`config_presence.R`).
    - `config_presence_less_coex`: Parameters for low species richness
      scenario (`config_presence_less_coex.R`).
  - `landscape_plain`: Files for simulations in a homogeneous
    environment (100 grids).
  - `landscape_plain_900`: Files for simulations in a 900-grid
    landscape.

- **gen3sis_simulations_grid100_high_species**  
  100 replicates for high species richness in a 100-grid landscape.

  - `gen3sis_output_dispersion_res_100_highdiv`: Processed outputs
    (.rds) with site composition and phylogenetic trees.

- **gen3sis_simulations_high_diversity_grid900**  
  100 replicates for high species richness in a 900-grid landscape.

  - `gen3sis_output_high_diversity_grid900`: Processed outputs (.rds)
    with site composition and phylogenetic trees.

- **gen3sis_simulations_low_diversity_grid100**  
  100 replicates for low species richness in a 100-grid landscape.

  - `gen3sis_output_lowdiv_100`: Processed outputs (.rds) with site
    composition and phylogenetic trees.

- **gen3sis_simulations_low_diversity_grid900**  
  100 replicates for low species richness in a 900-grid landscape.

  - `gen3sis_output_low_diversity_grid900`: Processed outputs (.rds)
    with site composition and phylogenetic trees.

- **empirical data**  
  Data and code for empirical analyses of ancestral area reconstruction
  and method comparison (SBEARS, rase, BioGeoBEARS) at different spatial
  grains.

  - `Script_correlation_empirical_methods.R`: Code for analyses (see
    Figure 4 in the manuscript).
  - `find_threshold_v2.R`: Auxiliary function for threshold definition
    in dispersal.
  - `function_SBEARS_spat_v5.R`: Computes SBEARS results.
  - Data folders by spatial scale:
    - `Escala 12 celulas`
    - `Escala 20 celulas`
    - `Escala 46 celulas`
    - `Escala 98 celulas`
    - `Escala 512 celulas`
    - `Escala 996 celulas`
    - `Escala 2178 celulas`

- **simulated data**  
  Data and code for simulated ancestral area reconstruction and method
  comparison.

  - `R_Script_BIOGEOBEARS_on_simulated_data.R`: Extracts species
    presence/absence using BioGeoBEARS.
  - `R_Script_BIOGEOBEARS_on_simulated_data.R`: Extracts species
    presence/absence using rase.
  - `save_rase_BioGeoBEARS_data.R`: Formats data for BioGeoBEARS and
    rase.
  - `sim_metacomm_featrues.R`: Reads gen3sis simulation output.
  - `R_function_BioGeoBEARS_on_simulated_data_v2.R`: Runs BioGeoBEARS on
    simulated data.
  - `R_function_rase_on_simulated_data_v5.R`: Runs rase on simulated
    data.
  - `function_SBEARS_on_simulated_data.R`: Compares simulated data and
    SBEARS reconstruction.
  - `function_SBEARS_spat_v5.R`: Computes SBEARS results.

------------------------------------------------------------------------

For questions or issues, please open an issue
