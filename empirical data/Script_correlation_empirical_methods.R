library(Herodotools)
source("function_SBEARS_spat_v5.R")
source("find_threshold_v2.R")

# Enter data (Input example using 12-sized cell grid)
comm<-read.table("comm_atual_12cells.txt",h=T)
phy<-ape::read.tree("tree260sp.tre")
coords<-read.table("comm_coords_12cells.txt",h=T)

# biogeobears reconstruction
biogeo<-read.table("comm_nodes_biogeobears_12cells.txt",h=T)
y=biogeo
rownames(y)<-c(1:12)
colnames(y)<-c(1:259)
dim(y)

# rase reconstruction
rase<-read.table("comm_nodes_rase_12cells.txt",h=T)
z<-rase
rownames(z)<-c(1:12)
colnames(z)<-c(1:259)


teste_single<-sbears(x=comm,
                     phy=phy,
                     coords=coords,
                     method = "single_site",
                     compute.node.by.sites = FALSE)

x<-t(teste_single$reconstruction)
threshold_single<-find_threshold(x)
sel_threshold<-threshold_single$Select_threshold
# Define discrete ada
threshold=0.65
x<-ifelse(x>=threshold,1,0)
dim(x)
rownames(x)<-c(1:12)
colnames(x)<-c(1:259)

#Procrustes test
res_tab_single<-matrix(NA,2,3,dimnames=list(c("per_site","per_node"),c("ada_biogeo","ada_rase","biogeo_rase")))
#per site
proc_ada_rase<-vegan::protest(x,z,symmetric = T,permutations = 1)
res_tab_single[1,2]<-sqrt(1-proc_ada_rase$ss)
proc_ada_biogeo<-vegan::protest(x,y,symmetric = T,permutations = 1)
res_tab_single[1,1]<-sqrt(1-proc_ada_biogeo$ss)
proc_rase_biogeo<-vegan::protest(z,y,symmetric = T,permutations = 1)
res_tab_single[1,3]<-sqrt(1-proc_rase_biogeo$ss)

#per node
proc_ada_rase_trans<-vegan::protest(t(x),t(z),symmetric = T,permutations = 1)
res_tab_single[2,2]<-sqrt(1-proc_ada_rase_trans$ss)
proc_ada_biogeo_trans<-vegan::protest(t(x),t(y),symmetric = T,permutations = 1)
res_tab_single[2,1]<-sqrt(1-proc_ada_biogeo_trans$ss)
proc_rase_biogeo_trans<-vegan::protest(t(z),t(y),symmetric = T,permutations = 1)
res_tab_single[2,3]<-sqrt(1-proc_rase_biogeo_trans$ss)

write.table(res_tab_single,"Res_Single12.txt",sep=" ")

# Test using w_slope = 10

teste_w10<-sbears(x=comm,
                            phy=phy,
                            coords=coords,
                            w_slope=10,
                            method="disp_assembly",
                            compute.node.by.sites = FALSE)# gave error.

x<-t(teste_w10$reconstruction)
threshold_w10<-find_threshold(x)
sel_threshold<-threshold_w10$Select_threshold
# Define discrete ada
threshold=0.7
x<-ifelse(x>=threshold,1,0)
dim(x)
rownames(x)<-c(1:12)
colnames(x)<-c(1:259)

#Procrustes test
res_tab_w10<-matrix(NA,2,3,dimnames=list(c("per_site","per_node"),c("ada_biogeo","ada_rase","biogeo_rase")))
#per site
proc_ada_rase<-vegan::protest(x,z,symmetric = T,permutations = 1)
res_tab_w10[1,2]<-sqrt(1-proc_ada_rase$ss)
proc_ada_biogeo<-vegan::protest(x,y,symmetric = T,permutations = 1)
res_tab_w10[1,1]<-sqrt(1-proc_ada_biogeo$ss)
proc_rase_biogeo<-vegan::protest(z,y,symmetric = T,permutations = 1)
res_tab_w10[1,3]<-sqrt(1-proc_rase_biogeo$ss)

#per node
proc_ada_rase_trans<-vegan::protest(t(x),t(z),symmetric = T,permutations = 1)
res_tab_w10[2,2]<-sqrt(1-proc_ada_rase_trans$ss)
proc_ada_biogeo_trans<-vegan::protest(t(x),t(y),symmetric = T,permutations = 1)
res_tab_w10[2,1]<-sqrt(1-proc_ada_biogeo_trans$ss)
proc_rase_biogeo_trans<-vegan::protest(t(z),t(y),symmetric = T,permutations = 1)
res_tab_w10[2,3]<-sqrt(1-proc_rase_biogeo_trans$ss)

