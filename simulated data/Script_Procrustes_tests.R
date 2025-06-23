source("function_SBEARS_on_simulated_data.R")
source("sim_metacomm_features.R")


#SBEARS

#Basic features of the metacommunity
data=#define path
metacomm_features<-sim_met(data=data,n_sim=100)
hist(metacomm_features[,1])
hist(metacomm_features[,2])
plot(metacomm_features[,4],metacomm_features[,3])

# Run Procrustes tests
test_single<-sim_test(data,n_sim=100,method="single_site")
test_single$Correlations_per_site
test_single$Correlations_per_node
test_single$Reconstruction

w_slope=10 # define w_slope
test_disp_assembly<-sim_test(data,n_sim=100,method="disp_assembly",min_disp_prob=0.8,w_slope=w_slope)
test_disp_assembly$Correlations_per_site
test_disp_assembly$Correlations_per_node
test_disp_assembly$Reconstruction

#BioGeoBEARS

# Enter and organize data
n_sim=100
res_biogeo<-readRDS("all_biogeobears_results.rds")
order_names<-c("output_dispersion_low_spp_100_1",
               "output_dispersion_low_spp_100_2",
               "output_dispersion_low_spp_100_3",
               "output_dispersion_low_spp_100_4",
               "output_dispersion_low_spp_100_5",
               "output_dispersion_low_spp_100_6",
               "output_dispersion_low_spp_100_7",
               "output_dispersion_low_spp_100_8",
               "output_dispersion_low_spp_100_9",
               "output_dispersion_low_spp_100_10",
               "output_dispersion_low_spp_100_11",
               "output_dispersion_low_spp_100_12",
               "output_dispersion_low_spp_100_13",
               "output_dispersion_low_spp_100_14",
               "output_dispersion_low_spp_100_15",
               "output_dispersion_low_spp_100_16",
               "output_dispersion_low_spp_100_17",
               "output_dispersion_low_spp_100_18",
               "output_dispersion_low_spp_100_19",
               "output_dispersion_low_spp_100_20",
               "output_dispersion_low_spp_100_21",
               "output_dispersion_low_spp_100_22",
               "output_dispersion_low_spp_100_23",
               "output_dispersion_low_spp_100_24",
               "output_dispersion_low_spp_100_25",
               "output_dispersion_low_spp_100_26",
               "output_dispersion_low_spp_100_27",
               "output_dispersion_low_spp_100_28",
               "output_dispersion_low_spp_100_29",
               "output_dispersion_low_spp_100_30",
               "output_dispersion_low_spp_100_31",
               "output_dispersion_low_spp_100_32",
               "output_dispersion_low_spp_100_33",
               "output_dispersion_low_spp_100_34",
               "output_dispersion_low_spp_100_35",
               "output_dispersion_low_spp_100_36",
               "output_dispersion_low_spp_100_37",
               "output_dispersion_low_spp_100_38",
               "output_dispersion_low_spp_100_39",
               "output_dispersion_low_spp_100_40",
               "output_dispersion_low_spp_100_41",
               "output_dispersion_low_spp_100_42",
               "output_dispersion_low_spp_100_43",
               "output_dispersion_low_spp_100_44",
               "output_dispersion_low_spp_100_45",
               "output_dispersion_low_spp_100_46",
               "output_dispersion_low_spp_100_47",
               "output_dispersion_low_spp_100_48",
               "output_dispersion_low_spp_100_49",
               "output_dispersion_low_spp_100_50",
               "output_dispersion_low_spp_100_51",
               "output_dispersion_low_spp_100_52",
               "output_dispersion_low_spp_100_53",
               "output_dispersion_low_spp_100_54",
               "output_dispersion_low_spp_100_55",
               "output_dispersion_low_spp_100_56",
               "output_dispersion_low_spp_100_57",
               "output_dispersion_low_spp_100_58",
               "output_dispersion_low_spp_100_59",
               "output_dispersion_low_spp_100_60",
               "output_dispersion_low_spp_100_61",
               "output_dispersion_low_spp_100_62",
               "output_dispersion_low_spp_100_63",
               "output_dispersion_low_spp_100_64",
               "output_dispersion_low_spp_100_65",
               "output_dispersion_low_spp_100_66",
               "output_dispersion_low_spp_100_67",
               "output_dispersion_low_spp_100_68",
               "output_dispersion_low_spp_100_69",
               "output_dispersion_low_spp_100_70",
               "output_dispersion_low_spp_100_71",
               "output_dispersion_low_spp_100_72",
               "output_dispersion_low_spp_100_73",
               "output_dispersion_low_spp_100_74",
               "output_dispersion_low_spp_100_75",
               "output_dispersion_low_spp_100_76",
               "output_dispersion_low_spp_100_77",
               "output_dispersion_low_spp_100_78",
               "output_dispersion_low_spp_100_79",
               "output_dispersion_low_spp_100_80",
               "output_dispersion_low_spp_100_81",
               "output_dispersion_low_spp_100_82",
               "output_dispersion_low_spp_100_83",
               "output_dispersion_low_spp_100_84",
               "output_dispersion_low_spp_100_85",
               "output_dispersion_low_spp_100_86",
               "output_dispersion_low_spp_100_87",
               "output_dispersion_low_spp_100_88",
               "output_dispersion_low_spp_100_89",
               "output_dispersion_low_spp_100_90",
               "output_dispersion_low_spp_100_91",
               "output_dispersion_low_spp_100_92",
               "output_dispersion_low_spp_100_93",
               "output_dispersion_low_spp_100_94",
               "output_dispersion_low_spp_100_95",
               "output_dispersion_low_spp_100_96",
               "output_dispersion_low_spp_100_97",
               "output_dispersion_low_spp_100_98",
               "output_dispersion_low_spp_100_99",
               "output_dispersion_low_spp_100_100")
res_biogeo_sort<-res_biogeo[order_names]

# Run Procrustes
res_procrustes_biogeo<-matrix(NA,n_sim,2,dimnames=list(1:n_sim,c("Correlations_per_site","Correlations_per_node")))
  for (i in 1:n_sim){
    res_sim_comm<-readRDS(paste("output_dispersion_low_spp_100_",i,".rds",sep = ""))
    occ_nodes_spp <- res_sim_comm$node_site_matrix
    phy<- res_sim_comm$tree
    phy$edge.length<-phy$edge.length+0.001
    real_node.anc.area<- occ_nodes_spp[, -match(phy$tip.label, colnames(occ_nodes_spp))]
    dim(real_node.anc.area)
    biogeo<-t(res_biogeo_sort[[i]]$node_site_presence_absence[-(1:length(phy$tip.label)),])
    prot_site<-vegan::protest(real_node.anc.area,biogeo,symmetric = T,permutations = 1)
    res_procrustes_biogeo[i,1]<-sqrt(1-prot_site$ss)
    prot_node<-vegan::protest(t(real_node.anc.area),t(biogeo),symmetric = T,permutations = 1)
    res_procrustes_biogeo[i,2]<-sqrt(1-prot_node$ss)
    print(i)
  }

# Extract the best model
mod_sel_name<-vector()
for (i in 1:n_sim){
  mod_sel_name[i]<-res_biogeo_sort[[i]]$best_model_name
}

# rase

# Enter and organize data
data= #enter data path to rase results ("all_rase_comparison_results.rds")
res_rase<-readRDS(data)
n_sim<-length(res_rase)
res_procrustes_rase<-matrix(NA,n_sim,2,dimnames=list(1:n_sim,c("Correlations_per_site","Correlations_per_node")))

# Run Procrustes
proc_site<-purrr::possibly(function(x,y){
  test<-vegan::protest(x,y,symmetric = T,permutations = 1)
  res<-sqrt(1-test$ss)
},otherwise=NA)
proc_node<-purrr::possibly(function(x,y){
  test<-vegan::protest(x,y,symmetric = T,permutations = 1)
  res<-sqrt(1-test$ss)
},otherwise=NA)

for (i in 1:n_sim){
  rownames(res_rase[[i]]$rase_node.anc.area)<-rownames(res_rase[[i]]$real_node.anc.area)
  res_procrustes_rase[i,1]<-proc_site(res_rase[[i]]$real_node.anc.area,res_rase[[i]]$rase_node.anc.area)
  res_procrustes_rase[i,2]<-proc_node(t(res_rase[[i]]$real_node.anc.area),t(res_rase[[i]]$rase_node.anc.area))
  print(i)
}





