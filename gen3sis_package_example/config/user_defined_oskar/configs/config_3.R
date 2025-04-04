### general settings ###
random_seed = 3269
start_time = NA
end_time = NA

max_number_of_species = 50000
max_number_of_coexisting_species = 1e5

trait_names = c("dispersal") # would be added by the package anyway
environmental_ranges = list() #"temp" = c(-45, 25)


end_of_timestep_observer = function(...){}


### initialization ###
initial_abundance = 1

create_initial_species <- function(landscape, config) {
  spacing <- 3220
  world <- raster::rasterFromXYZ(cbind(landscape$coordinates, landscape$environment[, 1, drop = F]))
  resolution <- raster::res(world)
  xs <- sort(unique(landscape$coordinates[, "x", drop = F]))
  ys <- sort(unique(landscape$coordinates[, "y", drop = F]))
  step_x <- spacing / (resolution[1] * 111)
  step_y <- spacing / (resolution[2] * 111)

  all_species <- list()

  for( cell in c("488", "496", "503", "510", "783") ) {
    new_species <- create_species(NA, cell, config)
    new_species$traits[ , "dispersal"] <- 1
    all_species <- append(all_species, list(new_species))
  }
  return(all_species)
}


### dispersal ###
max_dispersal <- 900
get_dispersal_values <- function(n, species, landscape, config) {
  dispersal_range = c(0, 900)
  disp_factor <- species[["traits"]][1, "dispersal"]
  scale <- ((dispersal_range[2] - dispersal_range[1]) * disp_factor ) + dispersal_range[1]

  values <- rnorm(n, mean = scale, sd=0)
  return(values)
}


### speciation ###
divergence_threshold = 6.9

a.mut <- function(t) {
  a_dist <- pnorm(t, mean=15, sd=15)*2
  return(a_dist)
}

get_divergence_factor = function(species, clu_geo_spi_ti, landscape, config) {
  temp_idi <- names(species[["abundance"]]) # get occurences
  t_env <- mean(landscape[["environment"]][temp_idi, "temp"]) #mean temp value of occurence
  ifactor <- config$user$a.mut(t_env) #get ifactor from a.mut function
 }


### mutation ###
source("../modules/mutation_mode_none.R", local = T)
apply_evolution <- mutation_mode_none


### ecology ###
source("../modules/ecology_none.R", local = T)
apply_ecology <- ecology_function_none

