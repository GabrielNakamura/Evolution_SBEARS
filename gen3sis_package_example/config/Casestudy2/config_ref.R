########################
### General settings ###
########################
# set the random seed for the simulation
random_seed = NA

# set the starting time step or leave NA to use the earliest/highest timestep
start_time = 390

# set the end time step or leave as NA to use the lates/lowest timestep (0)
end_time = NA

# maximum total number of species in the simulation before it is aborted
max_number_of_species = 10000

# maximum number of species within one cell before the simulation is aborted
max_number_of_coexisting_species = 5000

# a list of traits to include with each species
# a "dispersion" trait is implictly added in any case
#trait_names = c("t_min", "a_min", "competition", "dispersion")
trait_names = c("temp",  "dispersal") # "prec",

# ranges to scale the input environemts with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# lsited with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list("temp" = c(-45, 55))#, "prec" = NA )

# a place to inspect the internal state of the simulation and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  par(mfrow=c(2,2))
  plot_raster_single(data$landscape$environment[,"temp"], data$landscape, "temp", NA)
 # plot_raster_single(data$landscape$environment[,"prec"], data$landscape, "prec", NA)
  plot_richness(data$all_species, data$landscape)
  plot_species_presence(data$all_species[[1]], data$landscape)
  #plot_landscape(data$landscape)
  # a list of all species can be found in data$all_species
  # the current landscape can be found in data$landscape
}



######################
### Initialization ###
######################
# the intial abundace of a newly colonized cell, both during setup and later when colonizing a cell during the dispersal
initial_abundance =1

create_initial_species <- function(landscape, config) {
  spacing <- 2220
  world <- raster::rasterFromXYZ(cbind(landscape$coordinates, landscape$environment[, 1, drop = F])) 
  world <- extend(world, extent(-180, 180, -90, 90)) 
  
  #FOR ALL WORLD - START
  #     resolution <- raster::res(world) 
  #     xs <- sort(unique(landscape$coordinates[, "x", drop = F]))
  #     ys <- sort(unique(landscape$coordinates[, "y", drop = F]))
  #     step_x <- spacing / (resolution[1] * 111)
  #     step_y <- spacing / (resolution[2] * 111)
  #     
  #     all_species <- list()
  #     
  #     for( x in seq(1, length(xs), step_x)){
  #       force(x)
  #       for( y in seq(1, length(ys), step_y)){
  #         force(y)
  #         if( !is.na(world[y, x])) {
  #           cell <- as.character(raster::cellFromRowCol(world, y, x))
  #           new_species <- create_species(NA, as.character(cell), config)
  #           new_species$traits[ , "dispersal"] <- 1
  #           new_species$traits[ , "temp"] <- landscape$environment[cell, "temp"]
  #        #   new_species$traits[ , "prec"] <- landscape$environment[cell, "prec"]
  #           all_species <- append(all_species, list(new_species))
  #         }
  #       }
  #     }
  #FOR ALL WORLD - END
  
  #   #FOR TROPICS - START
       t <- raster('input/world_scotese/data_driven/1_koeppen_zone_bands/land_temp/1d/land_temp_0065.00.asc')
       p <- raster('input/world_scotese/data_driven/1_aridity/1_4_aridity_continental_rasters/065.00Ma.asc')
       #extent(t) <- c(-180, 180. -90, 90)
       #extent()
       #cat(res(t))
       tropics <- t
       tropics[!is.na(tropics)] <-0
       
       #tropics[t>25 & p<0.5] <- 1
       tropics[t>10 & p>0.5] <- 1
       tropics <- aggregate(tropics, fact=4)
       #plot(tropics)
       #all_species <- list()
       #get ID's of tropical cells

     #all_species <- list()
     #cell <- which(values(tropics==1)) #tropic IDs as integer/num
     #cell <- sample(cell, 50)
     #index <- as.character(cell) %in% rownames(landscape$coordinates)
     #cell <- cell[index]
     #new_species <- create_species(NA, as.character(cell), config) #as.character(cell) durch tropische zellen-liste ersetzen
     #new_species$traits[, "dispersal"] <- 1
     #new_species$traits[, "temp"] <- landscape$environment[as.character(cell), "temp"]
     ## new_species$traits[ , “prec”] <- landscape$environment[as.character(cell), “prec”]
     #all_species <- append(all_species, list(new_species))
       
       
       
       
       cell <- which(values(tropics==1)) #tropic IDs as integer/num
       cell <- sample(cell, 20)
       index <- as.character(cell) %in% rownames(landscape$coordinates)
       cell <- cell[index]
       all_species <- list()
       for(i in 1:length(cell)){
         current_cell <- cell[i]
         new_species <- create_species(NA, as.character(current_cell), config)
         new_species$traits[, "dispersal"] <- 1
         new_species$traits[, "temp"] <- landscape$environment[as.character(current_cell), "temp"]
         all_species <- append(all_species, list(new_species))
       }
       
       
      
  #FOR TROPICS - END 
  
  return(all_species)
}


#################
### Dispersal ###
#################
# returns n dispersal values
get_dispersal_values <- function(n, species, landscape, config) {
  dispersal_range = c(0,350)
  weighted_dispersal <- sum(species[["traits"]][, "dispersal"] * species[["abundance"]])
  disp_factor <- weighted_dispersal/sum(species[["abundance"]])
  scale <- ((dispersal_range[2] - dispersal_range[1]) * disp_factor ) + dispersal_range[1]
  
  values <- rweibull(n, shape = 5, scale = scale)
  return(values)
}


##################
### Speciation ###
##################
# threshold for genetic distance after which a speciation event takes place
divergence_threshold =72

# factor by which the divergence is increased between geographicaly isolated population
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  return(1)
}


################
### Mutation ###
################
# mutate the traits of a species and return the new traits matrix
apply_evolution <- function(species, cluster_indices, landscape, config) {
  traits <- species[["traits"]]
  traits[, "temp"] <- traits[, "temp"] + rnorm(length(traits[, "temp"]), mean = 0, sd =0.001) #evt. erhöhen
  #traits[, "prec"] <- traits[, "prec"] + rnorm(length(traits[, "prec"]), mean = 0, sd =0.001) #evt. erhöhen
  return(traits)
}


###############
### Ecology ###
###############
# called for every cell with all occuring species, this functin calculates the who survives in the current cells
# returns a vector of abundances
# set the abundance to 0 for every species supposed to die
apply_ecology <- function(abundance, traits, landscape, config) {
  abundance <- ( abs( traits[, "temp"] - landscape[, "temp"]) <0.02 ) #& #evt. erhöhen
 #   ( abs( traits[, "prec"] - landscape[, "prec"]) < 0.3 ) #evt. erhöhen
  return(abundance)
}

