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
  do.call(rbind, lapply(1:length(sim$presence), function(a){
    sim$presence[[a]] -> aa
    if(k%in%colnames(aa)&&any(aa[,k]!=0))
      data.frame(species = k, aa[which(aa[, k] != 0), c("x", "y")],
                 age = names(sim$presence)[a])
  }))
})
names(genocc) <- gentree$tip.label
gentree <- keep.tip(gentree, names(which(!sapply(genocc, is.null))))
genocc <- genocc[which(!sapply(genocc,is.null))]
sapply(genocc,nrow) -> totspocc
if(any(totspocc<10)){
  drop.tip(gentree,names(totspocc)[which(totspocc<10)])->gentree
  genocc[which(!names(genocc)%in%names(genocc)[which(totspocc<10)])]->genocc
}

# speciation mode - using only branching ----------------------------------

genocc -> startocc
gentree -> tree
treetime -> spectime
if(any(c("branching","anagenesis")%in%sp.type)){
  require(RRphylo)
  require(phytools)
  tree$edge[which(tree$edge[,2]<=Ntip(tree)),]->tipanc
  tipanc[which(tipanc[,1]%in%names(which(table(tipanc[,1])>1))),] -> tipanc
  lapply(split(tipanc,tipanc[,1]),function(x) tree$tip.label[x[3:4]]) -> tippair
  spcladogen <- data.frame(ancestor=do.call(rbind,tippair)[,1],
                         sptime=sapply(tippair,function(j){
                           spectime[which(apply(spectime[,1:2],1,function(k)
                             all(k%in%j))),3]->abc
                           if(length(abc)>0) abc else NA
                         }))
  tippair[which(!is.na(spcladogen[,2]))]->tippair
  tippair[which(spcladogen[,2]<3)]->tippair
  spcladogen[which(!is.na(spcladogen[,2])),,]->spcladogen
  spcladogen[which(spcladogen[,2]<3),,]->spcladogen
  if(sp.type=="branching"){#### Simulating branching cladogenesis ####
    sapply(1:nrow(spcladogen),function(j){
      which(names(startocc)==spcladogen[j,1])->ind
      startocc[[ind]][which(startocc[[ind]][,4]<=as.numeric(spcladogen[j,2])),1]<<-
        paste(spcladogen[j,1],"a",sep="")
      which(as.numeric(names(sim$presence))<=as.numeric(spcladogen[j,2]))->simpre
      sapply(simpre,function(w){
        sim$presence[[w]]->aa
        if(spcladogen[j,1]%in%colnames(aa))
          colnames(sim$presence[[w]])[which(colnames(aa)==spcladogen[j,1])]<<-
            paste(spcladogen[j,1],"a",sep="")
      })
      spectime[which(spectime[,1]==spcladogen[j,1]),1]<<-
        paste(spcladogen[j,1],"a",sep="")
      tree$tip.label[which(tree$tip.label==spcladogen[j,1])]<<-
        paste(spcladogen[j,1],"a",sep="")
      getMommy(tree,paste(spcladogen[j,1],"a",sep=""))[1]->mom
      tree$edge.length[which(tree$edge[,2]==mom)]/2->len
      tree<<-bind.tip(tree,spcladogen[j,1],where=mom,
                      edge.length = len, position = len)
    })
    split(do.call(rbind,startocc),do.call(rbind,startocc)[,1])->startocc
    startocc<-startocc[match(tree$tip.label,names(startocc))]
  }
}
