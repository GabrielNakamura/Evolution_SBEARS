required<-c("ape","sp","scales","rgeos","adehabitatHS","raster","sf","terra","parallel",
            "spatialEco","doParallel","ggplot2","cowplot","doSNOW","gdistance",
            "leastcostpath","ks","devtools","gen3sis","pbapply")
if(any(!required%in%installed.packages()[,1]))
  install.packages(required[which(!required%in%installed.packages()[,1])])
packs <- c('Rphylopars', 'RRphylo', 'dismo', 'ecospat', 'PresenceAbsence', 'biomod2')
install.packages(packs)

install.packages("RRgeo_0.0.1.tar.gz", repos = NULL, type = "source")



# genesis configuration ---------------------------------------------------

gen.conf <- "30_seed2"
sp.type <- "branching"

#datapath <- file.path(getwd(),"gen3sis/SouthAmerica/")
datapath <- here::here("gen3sis", "SouthAmerica")
file.show(here::here("gen3sis", "SouthAmerica", "config", "config_southamerica30_seed2.R"))


# Running gen3sis simulations ---------------------------------------------

require(ape)
require(gen3sis)
#source("R_functions/run_simulation_custom.R")
source(here::here("R_functions", "run_simulation_custom.R"))
#output.dir<-file.path(getwd(),"gen3sis/gen3sis outputs/")
output.dir <- here::here("gen3sis", "gen3sis_outputs")

sim <-
  run_simulation(config = file.path(datapath,"config",
                                    paste("config_southamerica", gen.conf, ".R",sep="")),
                 landscape = file.path(datapath, "landscape"),
                 output_directory = output.dir, verbose = 3)
save(sim,file=file.path(getwd(),"gen3sis",
                        paste("gen3sis output start",gen.conf,".Rda",sep="")))


# extracting gen3sis output -----------------------------------------------

gen.out.dir <- sim$sgen3sis$parameters$directories$output
read.table(paste(gen.out.dir, "phy.txt", sep="/"))[-1, , drop=FALSE] -> treetime
treetime[, 1:2] <- t(apply(treetime[, 1:2], 1, function(x) paste("species", x, sep="")))
read.nexus(paste(gen.out.dir, "phy.nex", sep = "/")) -> gentree
genocc <- lapply(gentree$tip.label, function(k){
  do.call(rbind,lapply(1:length(sim$presence),function(a){
    sim$presence[[a]]->aa
    if(k%in%colnames(aa)&&any(aa[,k]!=0))
      data.frame(species=k,aa[which(aa[,k]!=0),c("x","y")],
                 age=names(sim$presence)[a])
  }))
})
names(genocc)<-gentree$tip.label
gentree<-keep.tip(gentree,names(which(!sapply(genocc,is.null))))
genocc<-genocc[which(!sapply(genocc,is.null))]
sapply(genocc,nrow)->totspocc
if(any(totspocc<10)){
  drop.tip(gentree,names(totspocc)[which(totspocc<10)])->gentree
  genocc[which(!names(genocc)%in%names(genocc)[which(totspocc<10)])]->genocc
}

