
ecology_function_fast_cpp_solver <- function(abundance_scale, abundance_threshold) {
  force(abundance_scale)
  force(abundance_threshold)
  return(function(abundance, traits, landscape, config) {
    col_names <- colnames(traits)
    params <- cbind(abundance, traits, landscape[1, "temp"])
    colnames(params) <- c("abd", col_names, "temp")
    init <- rep(1, nrow(params))
    config$user$solver_cpp_vlr10k(init, params[,-1, drop = F], abundance_scale, abundance_threshold)
  })
}


solver_cpp_vlr10k <- function(init, traits, abundance_scale, abundance_threshold) {
  par_c=0.2
  init <- init
  N <- length(init)
  g <- (1-traits[, "t_min"])*(traits[, "temp"] -1+traits[, "t_min"])
  c <- rep((1-par_c)/abundance_scale, N)
  l <- (1-traits[, "competition"])/abundance_scale

  t0 <- 0
  tf <- 10000
  dt <- 0.1 #probably useless for adaptive methods

  out_B <- integrateModelVectorized(init, g, c, l, t0, tf, dt)
  y <- unname(unlist((out_B[length(out_B$Time),-1])))
  y[y < abundance_threshold] <- 0
  return(y)
}
