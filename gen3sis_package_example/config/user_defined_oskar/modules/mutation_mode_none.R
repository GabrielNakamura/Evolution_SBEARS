#' mutation_mode_none
#'
#' @param species the current species
#' @param cluster_indices indices to assign cells to geographic clusters
#' @param landscape the current landscape
#' @param config the general config
#'
#' @export
mutation_mode_none <- function(species, cluster_indices, landscape, config){
  return(invisible(species[["traits"]]))
}
