
######################################
###            METADATA            ###
######################################
# gen3sis configuration
#
# Version: 1.0
#
# Author: Marcelo Henrique Schwade
#
# Date: 04.02.2025 (MM.DD.YYYY)
#
# Landscape: landscape
#
# Publications: these data will are used in a publication in colaboration with
# Leandro Duarte, Gabriel Nakamura, Renan Maestri, Fabricio Villalobos, Carolina
# Prauchner and Juliene costa.
#
# Description: this is a first simulation to generate data to test models of
# species ancestral range reconstruction.
#
#
######################################


######################################
###         General settings       ###
######################################

# set the random seed for the simulation.
random_seed = rnorm(3)

# set the starting time step or leave NA to use the earliest/highest time-step.
start_time = 900#8000

# set the end time step or leave as NA to use the latest/lowest time-step (0).
end_time = 0

# maximum total number of species in the simulation before it is aborted.
max_number_of_species = 1000

# maximum number of species within one cell before the simulation is aborted.
max_number_of_coexisting_species = 200

# a list of traits to include with each species
# a "dispersal" trait is implicitly added in any case
trait_names = c("temp")

# ranges to scale the input environments with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list() #"temp"=NA

# Parameters to change.
dispersal_rate <- 0.15

# Extinction probability (local).
extinction_prob <- 0

# Initial distribution of the ancestor species.
range_initial <- c(-5, 0, -5, 0) # longmin, longmax, latmin, latmax.


######################################
###            Observer            ###
######################################

# a place to inspect the internal state of the simulation and collect additional information if desired.
end_of_timestep_observer = function(data, vars, config){
  # the list of all species can be found in data$all_species
  # the current landscape can be found in data$landscape

  # Save data of the species without make a copy of landscape data.
  config <- dynGet("config")
  data <- dynGet("data")
  vars <-  dynGet("vars")

  dir.create(file.path(config$directories$output, "species"), showWarnings=FALSE, recursive=TRUE)
  species <- data$all_species
  saveRDS(object = species,
          file = file.path(config$directories$output, "species", paste0("species_t_", vars$ti, ".rds")))

  # Plot data at each time step.
  #par(mfrow=c(2,1)) # To plot a mosaic of two graphs.

  # Plot richness.
  plot_richness(data$all_species, data$landscape)
  # Plot species 1 distribution.
  #plot_species_presence(data$all_species[[1]], data$landscape)

  mtext("STATUS", 1)
}


######################################
###         Initialization         ###
######################################

# the initial abundance of a newly colonized cell, both during setup and later when
# colonizing a cell during the dispersal.
initial_abundance = 1

# place species in the landscape:
create_ancestor_species <- function(landscape, config) {
  # stop("create the initial species here")

  # Define the distribution range of the ancestral species.
  # longmin, longmax, latmin, latmax.
  range <- config$user$range_initial
  # Using an initial range fully distributed.
  # The speciation starts only after the first ancestor colonizing all the sites.
  # Then, it's more practical and simple to start with the ancestor fully dispersed.

  # Select the sites within distribution range.
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] & co[, "x"] <= range[2] &
    co[, "y"] >= range[3] & co[, "y"] <= range[4]

  initial_cells <- rownames(co)[selection]

  # Create the initial species in this sites.
  new_species <- create_species(initial_cells, config)

  # Set the species traits.
  new_species$traits[, "temp"] <- 0 #landscape$environment[initial_cells, "temp"]
  # Here, we set the initial trait of each population species as the
  # environmental condition on that local.

  # Return the list with species information.
  return(list(new_species))

}


######################################
###             Dispersal          ###
######################################

# the maximum range to consider when calculating the distances from local distance inputs.
max_dispersal <- 10*dispersal_rate # Cutting radius.
# Distance where the accumulated probability of dispersal beyond is lesser than 0.001.
# This is computed to an exponential dispersal kernel.

# returns n dispersal values.
get_dispersal_values <- function(n, species, landscape, config) {

  # Generate values of dispersal.
  # Here, dispersal values are ever 1 or 2, than the dispersal happen only
  # between first neighbors.
  values <- 1 + rbinom(n, 1, prob=config$user$dispersal_rate)

  # Other option (equivalent to above one using rbinom)
  # prob control the dispersal rate
  # prob(dist = 1): 0.6 (don't disperse);
  # prob(dist = 2): 0.4 (disperse)
  #sample(c(1,2), n, prob=c(1-config$user$dispersal_rate, config$user$dispersal_rate), replace=T)

  # Probabilistic dispersal using a dispersal kernel exponential.
  #values <- rexp(n, rate=(1/config$user$dispersal_rate)) # is this a proper rate?

  return(values)
}


######################################
###          Speciation            ###
######################################

# threshold for genetic distance after which a speciation event takes place.
divergence_threshold = 10#500

# Define the factor by which the divergence is increased between geographically isolated population.
# can also be a matrix between the different population clusters.
# `species`:
# `cluster_indices`:
# `landscape`:
# `config`:
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  # stop("calculate divergence factor here")
  # Return 1 for each time step at which is accumulated genetic incompatibilities (?).
  return(1)
}


######################################
###            Evolution           ###
######################################

# mutate the traits of a species and return the new traits matrix.
apply_evolution <- function(species, cluster_indices, landscape, config) {

  # Get the species traits.
  traits <- species[["traits"]]

  # Get the cells.
  cells <- rownames(traits)

  # Return the modified traits.
  return(traits)
}


######################################
###             Ecology            ###
######################################

# called for every cell with all occurring species, this function calculates abundances and/or
# who survives for each sites.
# returns a vector of abundances.
# set the abundance to 0 for every species supposed to die.
# `abundance`:
# `traits`:
# `environment`:
# `config`:
apply_ecology <- function(abundance, traits, environment, config) {

  # Raffle if species population go extinct or survive.
  abundance <- abundance*rbinom(length(abundance), 1, prob=(1 - config$user$extinction_prob))

  # Return the final abundance.
  return(abundance)
}


