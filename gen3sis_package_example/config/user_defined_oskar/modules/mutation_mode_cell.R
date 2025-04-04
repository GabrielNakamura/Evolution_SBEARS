
#' mutation_mode_cell
#'
#' @param evolving_traits a vector contianing the names of all traits to be evolved
#' @param trait_rules a list of rules to tradeoff, invert, and limit species traits
#' @param mutation_factor the chance for a mutation to occur
#' @param trait_evolutionary_power the amount of change per mutaiton
#'
#' @export
mutation_mode_cell <- function(evolving_traits,
                               trait_rules,
                               mutation_factor,
                               trait_evolutionary_power) {
  force(evolving_traits)
  force(trait_rules)
  force(mutation_factor)
  force(trait_evolutionary_power)
  mutation_mode_cell_internal <- function(species, cluster_indices, landscape, config) {
    length_cluster <- length(cluster_indices)
    traits_clusters <- split(as.data.frame(species[["traits"]]), cluster_indices)
    abundance_clusters <- split(species[["abundance"]], cluster_indices)
    traits_clusters <- mapply(config$user$tr.homogenize,
                              traits_clusters,
                              abundance_clusters,
                              MoreArgs = list(
                                tr_names = config$gen3sis$general$trait_names,
                                tr_rules = trait_rules,
                                config = config),
                              SIMPLIFY = FALSE
    )
    # getting probabilities of mutation for each cluster
    mutations <- rbinom(length_cluster, 1, mutation_factor)

    mutations <- as.logical(mutations)
    #bind clu_tr_spi #attention to the new order (wow, that sounds strange...)
    NEW_clu_tr_spi <- do.call(rbind, unname(traits_clusters))

    # converting traits from data.frame into matrix (alternative would be bypass split dependency on data.frame)
    NEW_clu_tr_spi <- as.matrix(NEW_clu_tr_spi)

    if (any(mutations)) { # if mutations hapened, then...
      # mutation mechanism ONE where all traits mutate
      num_mutation_values <- sum(mutations)*length(evolving_traits)
      mutation_deltas <- sample(c(-trait_evolutionary_power, trait_evolutionary_power),
                                num_mutation_values, replace=T)

      NEW_clu_tr_spi[mutations, evolving_traits] <- NEW_clu_tr_spi[mutations, evolving_traits, drop=F] + mutation_deltas

      # apply tradeoff.traits
      NEW_clu_tr_spi[mutations,] <- config$user$tradeoff.traits(traits=NEW_clu_tr_spi[mutations, ,drop=F],
                                                                rules = trait_rules$relation)

      #avoid traits outside the range
      NEW_clu_tr_spi[mutations,] <- config$user$limit.traits.range(traits = NEW_clu_tr_spi[mutations, ,drop=F],
                                                                   trait_names = config$gen3sis$general$trait_names,
                                                                   range = c(0,1))

      if (length(trait_rules$tr_plane)!=0) { # check if there is any trait off surface rules to be applied
        # calculate new tr sum
        sum_trs <- rowSums(NEW_clu_tr_spi[ ,trait_rules$tr_plane, drop=F])

        # to limit traits to the surface where sum of g_max,c and l = 1

        mut_and_above <- sum_trs > 1
        if (any(mut_and_above)){
          NEW_clu_tr_spi[mut_and_above,] <-
            config$user$bring.traits.to.surface(traits_sp= NEW_clu_tr_spi[mut_and_above,,drop=F],
                                                traits_plane=trait_rules$tr_plane)
        }

      } # end check if there is any trait off surface rules to be applied

      # apply tradeoff.traits
      NEW_clu_tr_spi[mutations,] <- config$user$tradeoff.traits(traits=NEW_clu_tr_spi[mutations, , drop=F],
                                                                rules = trait_rules$relation)

    }
    species_traits <- NEW_clu_tr_spi[order(as.numeric(rownames(NEW_clu_tr_spi))), ,drop=F]
    return(invisible(species_traits))
  }
  return(invisible(mutation_mode_cell_internal))
}



tradeoff.traits <- function(traits, rules){
  #enforce declared trait tradeoff relations
  n_tr_rules_relation <- length(rules)
  if (n_tr_rules_relation>=1){
    for (x in 1:n_tr_rules_relation){
      tr_dom=rules[[x]][1]
      tr_sub=rules[[x]][2]
      if (names(rules)[x]=="inverse"){
        traits[,tr_sub] <- 1- traits[,tr_dom]
      }else if (names(rules)[x]=="equal"){
        traits[,tr_sub] <- traits[,tr_dom]
      }else{
        stop("name of traid should be either 'inverse' or 'equal'")
      }
    } # end loop over trait rules relation...
  }
  return(traits)
}


limit.traits.range <- function(traits, trait_names, range){
  # we force each trait listed by traits to be fix to a range!
  traits[,trait_names][traits[,trait_names]<range[1]] <- range[1]
  traits[,trait_names][traits[,trait_names]>range[2]] <- range[2]
  return(traits)
}


bring.traits.to.surface <- function(traits_sp, traits_plane){
  #balance out traits in case they are above the traitoff surface
  # print("traits were above the trait tradeoff surface... bringing them back to planet earth ;)")
  sums <- rowSums(traits_sp[,traits_plane, drop=F])
  id <- sums >= 1
  for (traits_i in traits_plane){
    traits_sp[id,traits_i] <-  traits_sp[id, traits_i, drop=F]/sums[id]
  }
  return(traits_sp)
}

tr.homogenize <- function(clu_tr_spi_x, abundance, tr_names, tr_rules, config){
  # updating eco (traits and abundance.... here I homegenize the traits per cluster!)
  # clu_tr_spi_x is one of the lists of the geo_cluster of the species
  nrow_clu_tr_spi_x <- nrow(clu_tr_spi_x)
  clu_tr_spi_x <- as.matrix(clu_tr_spi_x)
  #if (nrow_clu_tr_spi_x>1) { #check if more than one line
  mean_abd <- mean(abundance)

  weight_abd <- abundance/mean_abd


  hom_clu_tr_spi <- colMeans(clu_tr_spi_x[,tr_names, drop=F]*weight_abd)



  if (length(tr_rules$tr_plane)!=0) { # check if there is any trait off surface rules to be applied
    # calculate new tr sum
    sum_trs <- sum(hom_clu_tr_spi[tr_rules$tr_plane, drop=F])
    if ( sum_trs > 1 ) {
      # to limit traits to the surface where sum of g_max,c and l = 1
      hom_clu_tr_spi <- config$user$bring.traits.to.surface(traits_sp=t(as.matrix(hom_clu_tr_spi)), traits_plane=tr_rules$tr_plane)
    }
  } # end check if there is any trait off surface rules to be applied

  # apply trait traits
  hom_clu_tr_spi <- tradeoff.traits(traits=hom_clu_tr_spi, rules = tr_rules$relation)

  #clu_tr_spi_x[,c(tr_names, "sum_trs")] <- rep(as.matrix(hom_clu_tr_spi), each = nrow_clu_tr_spi_x)
  clu_tr_spi_x[ , tr_names] <- (rep(as.matrix(hom_clu_tr_spi), each = nrow_clu_tr_spi_x)+clu_tr_spi_x[, tr_names])/2

  clu_tr_spi_x <- config$user$bring.traits.to.surface(traits_sp=clu_tr_spi_x, traits_plane=tr_rules$tr_plane)
  clu_tr_spi_x <- config$user$tradeoff.traits(traits=clu_tr_spi_x, rules = tr_rules$relation)

  #}# end check if more than one line
  return(as.data.frame(clu_tr_spi_x))
}
