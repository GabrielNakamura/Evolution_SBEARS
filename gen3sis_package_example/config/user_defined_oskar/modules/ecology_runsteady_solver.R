ecology_function_runsteady_solver <- function(species, landscape, config){
  params <- cbind(species, landscape)
  init <- rep(1, nrow(params))
  abundance <- runsteady(y=init,#tr_sp_ti_idi$"abd",
                         #fun=config$exp$eco$eco.equilibria,
                         fun = config$user$d_abundance,
                         parms=params[, -1, drop=F],
                         stol=1e-8,
                         mf = 10,
                         maxsteps = 1000,
                         abd_scale)$y
  return(abundance)
}