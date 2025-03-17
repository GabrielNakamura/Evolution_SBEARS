run_simulation<-function (config = NA, landscape = NA, output_directory = NA, 
                          timestep_restart = NA, save_state = NA, call_observer = "all", 
                          enable_gc = FALSE, verbose = 1) 
{require(gen3sis)
  system_time_start <- Sys.time()
  directories <- prepare_directories(config_file = config, 
                                     input_directory = landscape, output_directory = output_directory)
  if (is.na(config)[1]) {
    stop("please provide either a config file or a config object")
  }else if (is(config, "gen3sis_config")) {
    config[["directories"]] <- directories
  }else if (is(config, "character")) {
    file.copy(config, directories$output)
    config <- create_input_config(config_file = config)
    config[["directories"]] <- directories
  }else {
    stop("this is not a known config, please provide either a config file or a config object")
  }
  if (!verify_config(config)) {
    stop("config verification failed")
  }
  val <- list(data = list(), vars = list(), config = config)
  val$config <- gen3sis:::complete_config(val$config)
  val$config$gen3sis$general$verbose <- verbose
  val <- gen3sis:::setup_inputs(val$config, val$data, val$vars) 
  val <- gen3sis:::setup_variables(val$config, val$data, val$vars)
  val <- gen3sis:::setup_landscape(val$config, val$data, val$vars)
  val$data$landscape$id <- val$data$landscape$id + 1
  val <- gen3sis:::init_attribute_ancestor_distribution(val$config, 
                                                        val$data, val$vars)
  val <- gen3sis:::init_simulation(val$config, val$data, val$vars)
  val <- gen3sis:::init_summary_statistics(val$data, val$vars, val$config)
  if (is.na(call_observer)) {
    save_steps <- c(val$config$gen3sis$general$start_time, 
                    val$config$gen3sis$general$end_time)
  }else if (call_observer == "all") {
    save_steps <- val$config$gen3sis$general$start_time:val$config$gen3sis$general$end_time
  }else {
    steps <- as.integer(call_observer) + 2
    save_steps <- ceiling(seq(val$config$gen3sis$general$start_time, 
                              val$config$gen3sis$general$end_time, length.out = steps))
  }
  val$vars$save_steps <- save_steps
  val$vars$steps <- val$config$gen3sis$general$start_time:val$config$gen3sis$general$end_time
  if (!is.na(timestep_restart)) {
    val <- gen3sis:::restore_state(val, timestep_restart)
  }
  presence<-list()
  for (ti in val$vars$steps) {
    val$vars$n_new_sp_ti <- 0
    val$vars$n_ext_sp_ti <- 0
    val$vars$n_sp_added_ti <- 0
    val$vars$ti <- ti
    if (verbose >= 2) {
      cat("loop setup \n")
    }
    val <- gen3sis:::setup_landscape(val$config, val$data, val$vars)
    val <- gen3sis:::restrict_species(val$config, val$data, val$vars)
    val <- gen3sis:::setup_distance_matrix(val$config, val$data, val$vars)
    if (verbose >= 2) {
      cat("speciation \n")
    }
    val <- gen3sis:::loop_speciation(val$config, val$data, val$vars)
    habitable<-val$data$landscape$coordinates
    habitable<-data.frame(habitable,
                          sapply(val$data$all_species,function(x) 
                            as.numeric(rownames(habitable)%in%names(which(x$abundance>0)))))
    colnames(habitable)[-c(1,2)]<-paste("species",sapply(val$data$all_species,"[[",1),sep="")
    presence[[which(val$vars$steps==ti)]]<-habitable
    val <- gen3sis:::update1.n_sp.all_geo_sp_ti(val$config, val$data, 
                                                val$vars)
    val <- gen3sis:::update2.n_sp_alive.geo_sp_ti(val$config, val$data, 
                                                  val$vars)
    if (verbose >= 2) {
      cat("dispersal \n")
    }
    val <- gen3sis:::loop_dispersal(val$config, val$data, val$vars)
    if (verbose >= 2) {
      cat("evolution \n")
    }
    val <- gen3sis:::loop_evolution(val$config, val$data, val$vars)
    if (verbose >= 2) {
      cat("ecology \n")
    }
    val <- gen3sis:::loop_ecology(val$config, val$data, val$vars)
    if (verbose >= 2) {
      cat("end of loop updates \n")
    }
    val$vars$n_sp_alive <- sum(sapply(val$data$all_species, 
                                      function(sp) {
                                        ifelse(length(sp[["abundance"]]), 1, 0) 
                                      }))
    val$vars$n_sp <- length(val$data$all_species)
    val <- gen3sis:::update_extinction_times(val$config, val$data, 
                                             val$vars)
    if (verbose >= 1) {
      cat("step =", ti, ", species alive =", val$vars$n_sp_alive, 
          ", species total =", val$vars$n_sp, "\n")
    }
    if (val$vars$ti %in% val$vars$save_steps) {
      gen3sis:::call_main_observer(val$data, val$vars, val$config)
    }
    val <- gen3sis:::update_summary_statistics(val$data, val$vars, 
                                               val$config)
    gen3sis:::save_val(val, save_state)
    if (val$vars$n_sp_alive >= val$config$gen3sis$general$max_number_of_species) {
      val$vars$flag <- "max_number_species"
      print("max number of species reached, breaking loop")
      break
    }
    if (val$vars$flag == "max_number_coexisting_species") {
      print("max number of coexisting species reached, breaking loop")
      break
    }
  }
  names(presence)<-val$vars$steps
  if (verbose >= 0 & val$vars$flag == "OK") {
    cat("Simulation finished. All OK \n")
  }else if (verbose >= 0 & val$vars$flag == "max_number_species") {
    cat("Simulation finished. Early abort due to exceeding max number of species")
  }else if (verbose >= 0 & val$vars$flag == "max_number_coexisting_species") {
    cat("Simulation finished. Early abort due to exceeding max number of co-occuring species")
  }
  val <- gen3sis:::update.phylo(val$config, val$data, val$vars)
  write.table(val$data$phy, file = file.path(val$config$directories$output, 
                                             "phy.txt"), sep = "\t")
  gen3sis:::write_nex(phy = val$data$phy, label = "species", file.path(output_location = val$config$directories$output, 
                                                                       "phy.nex"))
  system_time_stop <- Sys.time()
  total_runtime <- difftime(system_time_stop, system_time_start, 
                            units = "hours")[[1]]
  gen3sis:::write_runtime_statisitics(val$data, val$vars, val$config, 
                                      total_runtime)
  sgen3sis <- gen3sis:::make_summary(val$config, val$data, val$vars, 
                                     total_runtime, save_file = TRUE)
  if (verbose >= 1) {
    cat("Simulation runtime:", total_runtime, "hours\n")
  }
  return(list(sgen3sis=sgen3sis,presence=presence))
}
