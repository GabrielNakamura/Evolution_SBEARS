
# external installation of Rtools may be required on windows!
# https://cran.r-project.org/bin/windows/Rtools/index.html

# load libraries
libs = c("BH",
         "raster",
         "sp",
         "stringr",
         "tools",
         "optparse")
sapply(libs, require, character.only = T)

#install.packages(libs)

#install.packages("rgasm_0.5.25.tar.gz", repos = NULL, type = "source")

library("rgasm")

#running rGams on the simple island landscape:
rgasm::run_simulation(config_file = paste(getwd(),"/config/my_experiment/simple_config.R",sep=""),
                      input_directory = paste(getwd(),"/input/my_experiment",sep=""),
                      output_directory = paste(getwd(),"/output/my_experiment",sep=""),
                      timestep_restart = NA,
                      save_intermediate_results = "all")



#Experiment 1: dynamic island with initial colonizers
Replicates<- 30
Timesteps<-140
Species_richness_TimeSeries<-matrix(0,nrow=Timesteps, ncol=Replicates)
for(i in 2:Replicates){
  rgasm::run_simulation(config_file = paste(getwd(),"/config/my_experiment/simple_config.R",sep=""),
                        input_directory = paste(getwd(),"/input/my_experiment",sep=""),
                        output_directory = paste(getwd(),"/output/my_experiment",sep=""),
                        timestep_restart = NA,
                        save_intermediate_results = "all")
  load(paste(getwd(),"/output/my_experiment/simple_config/sgasm.RData",sep=""))
  ##plot(sgasm$turnover[,3])
  Species_richness_TimeSeries[,i]<-sgasm$turnover[,3]  
}

Species_richness_TimeSeries2<-as.data.frame(Species_richness_TimeSeries)
stringtimestepsnames<-as.character(1:Replicates)
colnames(Species_richness_TimeSeries2)<-stringtimestepsnames
Species_richness_TimeSeries3<-Species_richness_TimeSeries2[,-which(Species_richness_TimeSeries2[90,]==0)]
summary(Species_richness_TimeSeries3)
Species_richness_TimeSeries3$Means<- rowMeans(Species_richness_TimeSeries3)
#plot(Species_richness_TimeSeries3$Means)
#Number of replicate with non-zero final communities:
Numberofreplicates<-ncol(Species_richness_TimeSeries3)-1
Numberofreplicates #14 for high phylo constraint (=environmental variability) #7 less phylo constraint (>environmental variability)

#dataframe:
Species_richness_TimeSeries_Dataframe<-matrix(0,nrow=Timesteps*(ncol(Species_richness_TimeSeries3)-1), ncol=3)
counting<-1
for(time in 1:nrow(Species_richness_TimeSeries3)){
  for(rep in 1:(ncol(Species_richness_TimeSeries3)-1)){
    Species_richness_TimeSeries_Dataframe[counting,1]<-rep
    Species_richness_TimeSeries_Dataframe[counting,2]<-time
    Species_richness_TimeSeries_Dataframe[counting,3]<-Species_richness_TimeSeries3[time, rep]
    
    counting<-counting+1
  }
}
Species_richness_TimeSeries_Dataframe<-as.data.frame(Species_richness_TimeSeries_Dataframe)
colnames(Species_richness_TimeSeries_Dataframe)<-c("Replicate", "Timestep", "Richness")
Species_richness_TimeSeries3$Means<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, mean)
Species_richness_TimeSeries3$SDs<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, sd)
par(mfrow=c(1,1))
#plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,12000),xlab="Time step", ylab="Number of species", bty="l")
#plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,600),xlab="Time step", ylab="Number of species", bty="l")
#plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,60),xlab="Time step", ylab="Number of species", bty="l")
plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,1300),xlab="Time step", ylab="Number of species", bty="l")
lines(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs, lwd=1)
lines(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs, lwd=1)


Mainline<-log(Species_richness_TimeSeries3$Means)
Mainline[which(Mainline<0)]<-0
plot(Mainline, las=1, type="l",lwd=4,ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
Higherline<-log(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs)
Higherline[which(Higherline<0)]<-0
lines(Higherline, lwd=1)

saveRDS(Species_richness_TimeSeries3, paste(getwd(),"/output/my_experiment/Richness_Timeseries_Experiment_1_newDivergThres_HigherPhyloConstraint.rds", sep=""))


#Experiment 2: dynamic island with initial colonizers, with dynamics including only area
Replicates<- 30
Timesteps<-140
Species_richness_TimeSeries<-matrix(0,nrow=Timesteps, ncol=Replicates)
for(i in 1:Replicates){
  rgasm::run_simulation(config_file = paste(getwd(),"/config/my_experiment/simple_config.R",sep=""),
                        input_directory = paste(getwd(),"/input/my_experiment_2",sep=""),
                        output_directory = paste(getwd(),"/output/my_experiment_2",sep=""),
                        timestep_restart = NA,
                        save_intermediate_results = "all")
  load(paste(getwd(),"/output/my_experiment_2/simple_config/sgasm.RData",sep=""))
  ##plot(sgasm$turnover[,3])
  Species_richness_TimeSeries[,i]<-sgasm$turnover[,3]  
}

Species_richness_TimeSeries2<-as.data.frame(Species_richness_TimeSeries)
stringtimestepsnames<-as.character(1:Replicates)
colnames(Species_richness_TimeSeries2)<-stringtimestepsnames
Species_richness_TimeSeries3<-Species_richness_TimeSeries2[,-which(Species_richness_TimeSeries2[80,]==0)]
Species_richness_TimeSeries3<-Species_richness_TimeSeries3[,-which(is.na(Species_richness_TimeSeries3[140,]))]
summary(Species_richness_TimeSeries3)
Species_richness_TimeSeries3$Means<- rowMeans(Species_richness_TimeSeries3)
#plot(Species_richness_TimeSeries3$Means)
#Number of replicate with non-zero final communities:
Numberofreplicates<-ncol(Species_richness_TimeSeries3)-1
Numberofreplicates  #10 for less Phylo constraint

#dataframe:
Species_richness_TimeSeries_Dataframe<-matrix(0,nrow=Timesteps*(ncol(Species_richness_TimeSeries3)-1), ncol=3)
counting<-1
for(time in 1:nrow(Species_richness_TimeSeries3)){
  for(rep in 1:(ncol(Species_richness_TimeSeries3)-1)){
    Species_richness_TimeSeries_Dataframe[counting,1]<-rep
    Species_richness_TimeSeries_Dataframe[counting,2]<-time
    Species_richness_TimeSeries_Dataframe[counting,3]<-Species_richness_TimeSeries3[time, rep]
    
    counting<-counting+1
  }
}
Species_richness_TimeSeries_Dataframe<-as.data.frame(Species_richness_TimeSeries_Dataframe)
colnames(Species_richness_TimeSeries_Dataframe)<-c("Replicate", "Timestep", "Richness")
Species_richness_TimeSeries3$Means<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, mean)
Species_richness_TimeSeries3$SDs<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, sd)
par(mfrow=c(1,1))
#plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,25000),xlab="Time step", ylab="Number of species", bty="l")
#plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,2600),xlab="Time step", ylab="Number of species", bty="l")
plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,37000),xlab="Time step", ylab="Number of species", bty="l")
lines(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs, lwd=1)
lines(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs, lwd=1)

plot(log(Species_richness_TimeSeries3$Means), las=1, type="l",lwd=4,ylim=c(0,10.1),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
lines(log(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs), lwd=1)

saveRDS(Species_richness_TimeSeries3, paste(getwd(),"/output/my_experiment_2/Richness_Timeseries_Experiment_2_HigherPhyloConstraint.rds", sep=""))




#Plotting experiments together:
Richness_Timeseries_Experiment_1 <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/output/my_experiment/Richness_Timeseries_Experiment_1_newDivergThres_HighPhyloConstraint.rds")
plot(log(Richness_Timeseries_Experiment_1$Means), las=1, type="l",lwd=4,ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Richness_Timeseries_Experiment_1$Means-Richness_Timeseries_Experiment_1$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
lines(log(Richness_Timeseries_Experiment_1$Means+Richness_Timeseries_Experiment_1$SDs), lwd=1)

Richness_Timeseries_Experiment_2 <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/output/my_experiment_2/Richness_Timeseries_Experiment_2.rds")
lines(log(Richness_Timeseries_Experiment_2$Means), lwd=4, lty=2)
Lowerline<-log(Richness_Timeseries_Experiment_2$Means-Richness_Timeseries_Experiment_2$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1, lty=2)
lines(log(Richness_Timeseries_Experiment_2$Means+Richness_Timeseries_Experiment_2$SDs), lwd=1, lty=2)

Ltys<-c(1,2)
legends<-c("Area + heterogeneity dynamics", "Only area dynamics")
legend(x=-5,y=10, lty=Ltys, legend=legends, bty="n")

#Plotting experiments together: less phylo constraint
Richness_Timeseries_Experiment_1 <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/output/my_experiment/Richness_Timeseries_Experiment_1_newDivergThres_LessPhyloConstraint.rds")
Mainline<-log(Richness_Timeseries_Experiment_1$Means)
Mainline[which(Mainline<0)]<-0
plot(Mainline, las=1, type="l",lwd=4,ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Richness_Timeseries_Experiment_1$Means-Richness_Timeseries_Experiment_1$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
Higherline<-log(Richness_Timeseries_Experiment_1$Means+Richness_Timeseries_Experiment_1$SDs)
Higherline[which(Higherline<0)]<-0
lines(Higherline, lwd=1)

Richness_Timeseries_Experiment_2 <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/output/my_experiment_2/Richness_Timeseries_Experiment_2_LessPhyloConstraint.rds")
lines(log(Richness_Timeseries_Experiment_2$Means), lwd=4, lty=2)
Lowerline<-log(Richness_Timeseries_Experiment_2$Means-Richness_Timeseries_Experiment_2$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1, lty=2)
lines(log(Richness_Timeseries_Experiment_2$Means+Richness_Timeseries_Experiment_2$SDs), lwd=1, lty=2)

#Plotting experiments together: heiher phylo constraint
Richness_Timeseries_Experiment_1 <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/output/my_experiment/Richness_Timeseries_Experiment_1_newDivergThres_HigherPhyloConstraint.rds")
Mainline<-log(Richness_Timeseries_Experiment_1$Means)
Mainline[which(Mainline<0)]<-0
plot(Mainline, las=1, type="l",lwd=4,ylim=c(0,10.1),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Richness_Timeseries_Experiment_1$Means-Richness_Timeseries_Experiment_1$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
Higherline<-log(Richness_Timeseries_Experiment_1$Means+Richness_Timeseries_Experiment_1$SDs)
Higherline[which(Higherline<0)]<-0
lines(Higherline, lwd=1)

Richness_Timeseries_Experiment_2 <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/output/my_experiment_2/Richness_Timeseries_Experiment_2_HigherPhyloConstraint.rds")
lines(log(Richness_Timeseries_Experiment_2$Means), lwd=4, lty=2)
Lowerline<-log(Richness_Timeseries_Experiment_2$Means-Richness_Timeseries_Experiment_2$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1, lty=2)
lines(log(Richness_Timeseries_Experiment_2$Means+Richness_Timeseries_Experiment_2$SDs), lwd=1, lty=2)


#Checking landscape
all_geo_hab <- readRDS("~/workshops/sELDIG/Genesis/rgasm_beta/input/my_experiment/all_geo_hab.rds")
is.list(all_geo_hab)
is.data.frame(all_geo_hab)
ls(all_geo_hab)
all_geo_hab$temp[1:10,1:4]
dim(all_geo_hab$temp)
which(!is.na(all_geo_hab$temp[,3]))
summary(all_geo_hab$temp [,1:4])
length(levels(as.factor(as.numeric(all_geo_hab$temp [,1]))))
length(levels(as.factor(as.numeric(all_geo_hab$temp [,2]))))





