#' Auxiliary function for reading gen3sis simulation output
#'
#' @param data A string containing the path to the simulation output
#' @param n_sim A scalar indicating the number of simulations
#'
#' @returns
#' @export
#'
#' @examples
sim_met<-function(data, n_sim){
    res_table<-matrix(NA,n_sim,4,dimnames=list(1:n_sim,c("N_spp","Mean_richness","eveness_spp_dist", "n_doubletons")))
      for (i in 1:n_sim){
        res_sim_comm<-readRDS(paste(data,i,".rds",sep = ""))
        occ_nodes_spp <- res_sim_comm$node_site_matrix
        phy<- res_sim_comm$tree
        comm <- occ_nodes_spp[, match(phy$tip.label, colnames(occ_nodes_spp))]
        spp_dist<-colSums(comm)
        res_table[i,1]<-ncol(comm)
        res_table[i,2]<-mean(rowSums(comm))
        res_table[i,3]<-vegan::renyi(spp_dist,12)/vegan::renyi(spp_dist,0)
        res_table[i,4]<-length(which(colSums(comm)<=2))
        }
    return(res_table)
}
