library(Herodotools)
source("function_SBEARS_spat_v5.R")

  sim_test<-function(data="/Users/leandroduarte/Desktop/ADA/teste_metodos/Matrizes simuladas/gen3sis_output_lowdiv_100/data/output_dispersion_low_spp_100_",
                     n_sim,method=c("single_site","disp_assembly"),
                     w_slope=2){
    res_table_node<-res_table_site<-matrix(NA,n_sim,3,dimnames=list(1:n_sim,c("t=0.5","t=0.7","t=0.9")))
    res_test<-list()
      for (i in 1:n_sim){
        res_sim_comm<-readRDS(paste(data,i,".rds",sep = ""))
        occ_nodes_spp <- res_sim_comm$node_site_matrix
        phy<- res_sim_comm$tree
        phy$edge.length<-phy$edge.length+0.001
        coords <- res_sim_comm$coordinates
        comm <- occ_nodes_spp[, match(phy$tip.label, colnames(occ_nodes_spp))]
        real_node.anc.area<- occ_nodes_spp[, -match(phy$tip.label, colnames(occ_nodes_spp))]
        test<-sbears(x=comm,phy=phy,coords=coords,method = method,w_slope=w_slope,compute.node.by.sites = F)
        x<-t(test$reconstruction)
        res_test[[i]]<-test

        proc_site<-purrr::possibly(function(x,y){
          test<-vegan::protest(x,y,symmetric = T,permutations = 1)
          res<-sqrt(1-test$ss)
        },otherwise=NA)
        proc_node<-purrr::possibly(function(x,y){
          test<-vegan::protest(x,y,symmetric = T,permutations = 1)
          res<-sqrt(1-test$ss)
        },otherwise=NA)

        x_0.5<-ifelse(x>=0.5,1,0)
          if(max(x_0.5)==0){
            res_table_site[i,1]<-NA
            res_table_node[i,1]<-NA
          } else { x_0.5 = x_0.5
            prot_site<-proc_site(real_node.anc.area,x_0.5)
            res_table_site[i,1]<-prot_site
            prot_node<-proc_node(t(real_node.anc.area),t(x_0.5))
            res_table_node[i,1]<-prot_node
            }

        x_0.7<-ifelse(x>=0.7,1,0)
          if(max(x_0.7)==0){
            res_table_site[i,2]<-NA
            res_table_node[i,2]<-NA
          } else { x_0.7 = x_0.7
            prot_site<-proc_site(real_node.anc.area,x_0.7)
            res_table_site[i,2]<-prot_site
            prot_node<-proc_node(t(real_node.anc.area),t(x_0.7))
            res_table_node[i,2]<-prot_node
            }

        x_0.9<-ifelse(x>=0.9,1,0)
          if(max(x_0.9)==0){
            res_table_site[i,3]<-NA
            res_table_node[i,3]<-NA
          } else { x_0.9 = x_0.9
                   proc<-proc_site(real_node.anc.area,x_0.9)
                   res_table_site[i,3]<-prot_site
                   prot_node<-proc_node(t(real_node.anc.area),t(x_0.9))
                   res_table_node[i,3]<-prot_node
            }
        print(i)
      }
    return(list(Correlations_per_site=res_table_site,Correlations_per_node=res_table_node,Reconstruction=res_test))
  }
