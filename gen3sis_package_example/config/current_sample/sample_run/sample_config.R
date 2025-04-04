

########################
### General settings ###
########################
# set the random seed for the simulation
random_seed = 666

# set the starting time step or leave NA to use the earliest/highest timestep
start_time = NA

# set the end time step or leave as NA to use the lates/lowest timestep (0)
end_time = NA

# maximum total number of species in the simulation before it is aborted
max_number_of_species = 50000

# maximum number of species within one cell before the simulation is aborted
max_number_of_coexisting_species = 1e5

# a list of traits to include with each species
# a "dispersion" trait is implictly added in any case
#trait_names = c("t_min", "a_min", "competition", "dispersion")
trait_names = c("t_min", "p_min", "dispersal")

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list("temp" = c(-40, 60))

# a place to inspect the internal state of the simulation and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  # a list of all species can be found in data$all_species
  # the current landscape can be found in data$landscape
  save_landscape()
  save_species()
  save_traits()
}

######################
### Initialization ###
######################
# the intial abundace of a newly colonized cell, both during setup and later when colonizing a cell during the dispersal
initial_abundance = 1

# place species within rectangle:
create_ancestor_species <- function(landscape, config) {
  #range <- c(-180, 180, -90, 90)
  range <- c(-180, 180, -90, 90)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
               co[, "x"] <= range[2] &
               co[, "y"] >= range[3] &
               co[, "y"] <= range[4]
  initial_cells <- rownames(co)[selection]
  new_species <- create_species(NA, initial_cells, config)
  #t_min <- min(landscape$environment[initial_cells, "temp", drop = F])
  ### Here, you give the t_min and a_min (species trait)
  new_species$traits[ , "t_min"] <- 0.6
  new_species$traits[ , "p_min"] <- 0.6
  new_species$traits[ , "dispersal"] <- 1

  return(list(new_species))
}



#################
### Dispersal ###
#################
# returns n dispersal values
get_dispersal_values <- function(n, species, landscape, config) {
  dispersal_range = c(0, 666)
  disp_factor <- 1
  scale <- ((dispersal_range[2] - dispersal_range[1]) * disp_factor ) + dispersal_range[1]

  values <- rweibull(n, shape = 2.5, scale = scale)
  return(values)
}


##################
### Speciation ###
##################
# threshold for genetic distance after which a speciation event takes place
# speciation after every timestep : 0.9
divergence_threshold = 24

# ifactor by which the genetic distance is increased between geographicaly isolated population of a species
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}


################
### Mutation ###
################
# mutate the traits of a species and return the new traits matrix
source("../modules/mutation_mode_adaptive.R", local = T)
apply_evolution <- mutation_mode_adaptive(evolving_traits = c("t_min", "p_min"),
                                     trait_rules = list(relation=NULL, tr_planes=c("t_min", "p_min")),
                                     mutation_factor = 0.1,
                                     trait_evolutionary_power = 0.1)


###############
### Ecology ###
###############
# called for every cell with all occuring species, this functin calculates who survives in the current cells
# returns a vector of abundances
# set the abundance to 0 for every species supposed to die
apply_ecology <- function(abundance, traits, local_environment, config) {
  # browser()
  rt= 0.1
  rt2=rt/2
  #temp
  t_min <- traits[,"t_min"]
  t_max <- t_min+rt
  t_env <- local_environment[,"temp"] #env temp
  #equation---------
  #growth
  g <- (1/(rt2^2)*(t_env-(t_min+rt2))^2)
  #logical if inside range
  at <- as.numeric((t_min<=t_env & t_max>=t_env))
  #final biomasses
  abd <- g*(at*100)
  #end equation-------
  return(abd)
}


