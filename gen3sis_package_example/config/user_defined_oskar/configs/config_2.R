### general settings ###
random_seed = 326962
start_time = NA
end_time = 100

max_number_of_species = 50000
max_number_of_coexisting_species = 1e5

trait_names = c("t_min", "competition", "dispersal")
environmental_ranges = list("temp" = c(-45, 25) )


end_of_timestep_observer = function(...){}


### initialization ###
initial_abundance = 1

create_initial_species <- function(landscape, config) {
  spacing <- 2220
  world <- raster::rasterFromXYZ(cbind(landscape$coordinates, landscape$environment[, 1, drop = F]))
  resolution <- raster::res(world)
  xs <- sort(unique(landscape$coordinates[, "x", drop = F]))
  ys <- sort(unique(landscape$coordinates[, "y", drop = F]))
  step_x <- spacing / (resolution[1] * 111)
  step_y <- spacing / (resolution[2] * 111)

  all_species <- list()


  for( y in rev(seq(length(ys), 1, -step_y))){
    force(y)
      for( x in seq(1, length(xs), step_x)){
        force(x)
      if( !is.na(world[y, x])) {
        cell <- as.character(raster::cellFromRowCol(world, y, x))
        new_species <- create_species(NA, as.character(cell), config)
        t_min <- min(landscape$environment[cell, "temp"])
        new_species$traits[ , "dispersal"] <- 1
        new_species$traits[ , "t_min"] <- 1 - t_min + 0.1
        new_species$traits[ , "competition"] <- 0.8 - ( 1 - t_min + 0.1 )
        all_species <- append(all_species, list(new_species))
      }
    }
  }
  return(all_species)
}

### dispersal ###
get_dispersal_values <- function(n, species, landscape, config) {
  dispersal_range = c(0, 900)
  weighted_dispersal <- sum(species[["traits"]][, "dispersal"] * species[["abundance"]])
  disp_factor <- weighted_dispersal/sum(species[["abundance"]])
  scale <- ((dispersal_range[2] - dispersal_range[1]) * disp_factor ) + dispersal_range[1]

  values <- rweibull(n, shape = 2.5, scale = scale)
  return(values)
}


### speciation ###
divergence_threshold = 10.9

get_divergence_factor = function(...) { return(1) }


### mutation ###
source("../modules/mutation_mode_abundance.R", local = T)
apply_evolution <- mutation_mode_abundance(evolving_traits = c("t_min","competition"),
                                          trait_rules = list(relation=NULL, tr_planes=c("t_min", "competition")),
                                          mutation_factor = 0.2,
                                          trait_evolutionary_power = 0.05)


### ecology ###
source("../modules/ecology_vl_cpp_vectorized.R", local = T)
apply_ecology <- ecology_function_fast_cpp_solver(abundance_scale = 100,
                                                         abundance_threshold = 1)

