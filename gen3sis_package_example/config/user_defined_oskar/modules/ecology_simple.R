calculate_equilibria <- function(abundance, traits, landscape, config) {
  browser()
  rt= 0.1  # temperature tolerance
  ra= 1 # aridity/precipitation tolerance
  rt2=rt/2
  ra2=ra/2
  #temp
  tmin <- traits[,"t_min"]
  tmax <- tmin+rt
  tenv <- landscape[,"t_env"] #env temp
  #aridity
  amin <- traits[,"a_min"]
  amax <- amin+ra
  aenv <- landscape[,"a_env"] #env aridity
  #equation---------
  #growth
  g <- pmin(1-(1/(rt2^2)*(tenv-(tmin+rt2))^2),1-(1/(ra2^2)*(aenv-(amin+ra2))^2))
  #logical if inside range
  at <- as.numeric((tmin<=tenv & tmax>=tenv)&(amin<=aenv & amax>=aenv))
  #final biomasses
  abd <- g*at*100
  #end equation-------
  return(abd)
}
