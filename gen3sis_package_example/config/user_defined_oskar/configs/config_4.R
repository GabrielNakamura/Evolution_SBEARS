### general settings ###
random_seed = 326962
start_time = NA
end_time = 161

max_number_of_species = 50000
max_number_of_coexisting_species = 1e5

trait_names = c("t_min", "competition", "dispersal")
environmental_ranges = list("temp" = c(-45, 25) )


end_of_timestep_observer = function(...){}


### initialization ###
initial_abundance = 1

create_initial_species <- function(landscape, config) {
  range <- c(-180, 180, -90, 90)
  co <- landscape$coordinates
  selection <- co[, "x"] >= range[1] &
    co[, "x"] <= range[2] &
    co[, "y"] >= range[3] &
    co[, "y"] <= range[4]
  initial_cells <- rownames(co)[selection]
  new_species <- create_species(NA, initial_cells, config)
  t_min <- min(landscape$environment[initial_cells, "temp", drop = F])
  new_species$traits[ , "dispersal"] <- 1
  new_species$traits[ , "t_min"] <- 1 - t_min + 0.1
  new_species$traits[ , "competition"] <- 0.8 - ( 1 - t_min + 0.1 )
  return(list(new_species))
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
divergence_threshold = 2.9

get_divergence_factor = function(...) { 1 }


### mutation ###
source("../modules/mutation_mode_cell.R", local = T)
apply_evolution <- mutation_mode_cell(evolving_traits = c("t_min","competition"),
                                     trait_rules = list(relation=NULL, tr_planes=c("t_min", "competition")),
                                     mutation_factor = 0.2,
                                     trait_evolutionary_power = 0.05)


### ecology ###
source("../modules/ecology_vl_cpp_vectorized.R", local = T)
apply_ecology <- ecology_function_fast_cpp_solver(abundance_scale = 100,
                                                         abundance_threshold = 1)

