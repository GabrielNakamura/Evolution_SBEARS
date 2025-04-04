
# external installation of Rtools may be required on windows!
# https://cran.r-project.org/bin/windows/Rtools/index.html

# load libraries
libs = c("BH",
         "raster",
         "sp",
         "stringr",
         "tools",
         "optparse",
         "gdistance")
sapply(libs, require, character.only = T)

#install.packages(libs)

#install.packages("gen3sis_0.9.2.tar.gz", repos = NULL, type = "source")

#library("rgasm")
library("gen3sis")

#Experiment: dynamic island with initial colonizers, with dynamics including only area
Replicates<- 30
Timesteps<-140
#Scenario 1: high mutation rate
Species_richness_TimeSeries<-matrix(0,nrow=Timesteps, ncol=Replicates)
for(i in 1:Replicates){
  gen3sis::run_simulation(config_file = paste(getwd(),"/config/CaseStudy1/high_trait_evolution/config_hte_randomseed.R",sep=""),
                        input_directory = paste(getwd(),"/landscape/CaseStudy1",sep=""),
                        output_directory = paste(getwd(),"/output/CaseStudy1/high_trait_evolution",sep=""),
                        timestep_restart = NA,
                        save_intermediate_results = "all")
  load(paste(getwd(),"/output/CaseStudy1/high_trait_evolution/config_hte_randomseed/sgen3sis.RData",sep=""))
  Species_richness_TimeSeries[,i]<-sgen3sis$turnover[,3] 
  print(i)
}

Species_richness_TimeSeries2<-as.data.frame(Species_richness_TimeSeries)
stringtimestepsnames<-as.character(1:Replicates)
colnames(Species_richness_TimeSeries2)<-stringtimestepsnames
summary(Species_richness_TimeSeries2)
Species_richness_TimeSeries2$Means<- rowMeans(Species_richness_TimeSeries2)
#plot(Species_richness_TimeSeries2$Means)
#Number of replicate with non-zero final communities:
Species_richness_TimeSeries3<-Species_richness_TimeSeries2[,-which(Species_richness_TimeSeries2[80,]==0)]
Numberofreplicates<-ncol(Species_richness_TimeSeries3)-1
Numberofreplicates  #14 for higher mutation (less Phylo constraint)

#dataframe:
Species_richness_TimeSeries_Dataframe<-matrix(0,nrow=Timesteps*(ncol(Species_richness_TimeSeries3)-1), ncol=4)
counting<-1
for(time in 1:nrow(Species_richness_TimeSeries3)){
  for(rep in 1:(ncol(Species_richness_TimeSeries3)-1)){
    Species_richness_TimeSeries_Dataframe[counting,1]<-"HighTraitEvolution"
    Species_richness_TimeSeries_Dataframe[counting,2]<-rep
    Species_richness_TimeSeries_Dataframe[counting,3]<-time
    Species_richness_TimeSeries_Dataframe[counting,4]<-Species_richness_TimeSeries3[time, rep]
    
    counting<-counting+1
  }
}
Species_richness_TimeSeries_Dataframe<-as.data.frame(Species_richness_TimeSeries_Dataframe)
colnames(Species_richness_TimeSeries_Dataframe)<-c("Scenario","Replicate", "Timestep", "Richness")
Species_richness_TimeSeries_Dataframe$Replicate<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Replicate))
Species_richness_TimeSeries_Dataframe$Timestep<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Timestep))
Species_richness_TimeSeries_Dataframe$Richness<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Richness))
Species_richness_TimeSeries3$Means<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, mean)
Species_richness_TimeSeries3$SDs<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, sd)
par(mfrow=c(1,1))
plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,2200),xlab="Time step", ylab="Number of species", bty="l")
lines(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs, lwd=1)
lines(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs, lwd=1)

plot(log(Species_richness_TimeSeries3$Means), las=1, type="l",lwd=4,ylim=c(0,8),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
lines(log(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs), lwd=1)

saveRDS(Species_richness_TimeSeries3, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_HighTraitEvo_OnlySurvivingReps.rds", sep=""))
saveRDS(Species_richness_TimeSeries2, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_HighTraitEvo_AllReps.rds", sep=""))

Species_richness_TimeSeries4<-Species_richness_TimeSeries2 #to bind next results


#Scenario 2: reference mutation rate
Species_richness_TimeSeries<-matrix(0,nrow=Timesteps, ncol=Replicates)
for(i in 1:Replicates){
  gen3sis::run_simulation(config_file = paste(getwd(),"/config/CaseStudy1/ref_trait_evolution/config_rte_randomseed.R",sep=""),
                          input_directory = paste(getwd(),"/landscape/CaseStudy1",sep=""),
                          output_directory = paste(getwd(),"/output/CaseStudy1/ref_trait_evolution",sep=""),
                          timestep_restart = NA,
                          save_intermediate_results = "all")
  load(paste(getwd(),"/output/CaseStudy1/ref_trait_evolution/config_rte_randomseed/sgen3sis.RData",sep=""))
  Species_richness_TimeSeries[,i]<-sgen3sis$turnover[,3] 
  print(i)
}

Species_richness_TimeSeries2<-as.data.frame(Species_richness_TimeSeries)
stringtimestepsnames<-as.character(1:Replicates)
colnames(Species_richness_TimeSeries2)<-stringtimestepsnames
summary(Species_richness_TimeSeries2)
Species_richness_TimeSeries2$Means<- rowMeans(Species_richness_TimeSeries2)
#plot(Species_richness_TimeSeries2$Means)
#Number of replicate with non-zero final communities:
Species_richness_TimeSeries_Dataframe2<-Species_richness_TimeSeries_Dataframe[,-which(Species_richness_TimeSeries_Dataframe[,4]==0)]
Species_richness_TimeSeries3<-Species_richness_TimeSeries2[,-which(Species_richness_TimeSeries2[80,]==0)]
Numberofreplicates<-ncol(Species_richness_TimeSeries3)-1
Numberofreplicates  #27 for reference mutation (reference Phylo constraint)

#dataframe:
Species_richness_TimeSeries_Dataframe<-matrix(0,nrow=Timesteps*(ncol(Species_richness_TimeSeries3)-1), ncol=4)
counting<-1
for(time in 1:nrow(Species_richness_TimeSeries3)){
  for(rep in 1:(ncol(Species_richness_TimeSeries3)-1)){
    Species_richness_TimeSeries_Dataframe[counting,1]<-"ReferenceTraitEvolution"
    Species_richness_TimeSeries_Dataframe[counting,2]<-rep
    Species_richness_TimeSeries_Dataframe[counting,3]<-time
    Species_richness_TimeSeries_Dataframe[counting,4]<-Species_richness_TimeSeries3[time, rep]
    
    counting<-counting+1
  }
}
Species_richness_TimeSeries_Dataframe<-as.data.frame(Species_richness_TimeSeries_Dataframe)
colnames(Species_richness_TimeSeries_Dataframe)<-c("Scenario","Replicate", "Timestep", "Richness")
Species_richness_TimeSeries_Dataframe$Replicate<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Replicate))
Species_richness_TimeSeries_Dataframe$Timestep<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Timestep))
Species_richness_TimeSeries_Dataframe$Richness<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Richness))

Species_richness_TimeSeries3$Means<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, mean)
Species_richness_TimeSeries3$SDs<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, sd)
par(mfrow=c(1,1))
plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,15000),xlab="Time step", ylab="Number of species", bty="l")
lines(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs, lwd=1)
lines(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs, lwd=1)

plot(log(Species_richness_TimeSeries3$Means), las=1, type="l",lwd=4,ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
lines(log(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs), lwd=1)

saveRDS(Species_richness_TimeSeries3, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_RefTraitEvo_OnlySurvivingReps.rds", sep=""))
saveRDS(Species_richness_TimeSeries2, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_RefTraitEvo_AllReps.rds", sep=""))

Species_richness_TimeSeries4<-rbind(Species_richness_TimeSeries4,Species_richness_TimeSeries2) #to bind next results
summary(Species_richness_TimeSeries4)

#Scenario 3: low mutation rate
Species_richness_TimeSeries<-matrix(0,nrow=Timesteps, ncol=Replicates)
for(i in 1:Replicates){
  gen3sis::run_simulation(config_file = paste(getwd(),"/config/CaseStudy1/low_trait_evolution/config_lte_randomseed.R",sep=""),
                          input_directory = paste(getwd(),"/landscape/CaseStudy1",sep=""),
                          output_directory = paste(getwd(),"/output/CaseStudy1/low_trait_evolution",sep=""),
                          timestep_restart = NA,
                          save_intermediate_results = "all")
  load(paste(getwd(),"/output/CaseStudy1/low_trait_evolution/config_lte_randomseed/sgen3sis.RData",sep=""))
  Species_richness_TimeSeries[,i]<-sgen3sis$turnover[,3] 
  print(i)
}

Species_richness_TimeSeries2<-as.data.frame(Species_richness_TimeSeries)
stringtimestepsnames<-as.character(1:Replicates)
colnames(Species_richness_TimeSeries2)<-stringtimestepsnames
summary(Species_richness_TimeSeries2)
Species_richness_TimeSeries2$Means<- rowMeans(Species_richness_TimeSeries2)
#plot(Species_richness_TimeSeries2$Means)
#Number of replicate with non-zero final communities:
Species_richness_TimeSeries3<-Species_richness_TimeSeries2[,-which(Species_richness_TimeSeries2[80,]==0)]
Numberofreplicates<-ncol(Species_richness_TimeSeries3)-1
Numberofreplicates  #30 for reference mutation (low Phylo constraint)
Species_richness_TimeSeries3<-Species_richness_TimeSeries2

#dataframe:
Species_richness_TimeSeries_Dataframe<-matrix(0,nrow=Timesteps*(ncol(Species_richness_TimeSeries3)-1), ncol=4)
counting<-1
for(time in 1:nrow(Species_richness_TimeSeries3)){
  for(rep in 1:(ncol(Species_richness_TimeSeries3)-1)){
    Species_richness_TimeSeries_Dataframe[counting,1]<-"LowTraitEvolution"
    Species_richness_TimeSeries_Dataframe[counting,2]<-rep
    Species_richness_TimeSeries_Dataframe[counting,3]<-time
    Species_richness_TimeSeries_Dataframe[counting,4]<-Species_richness_TimeSeries3[time, rep]
    
    counting<-counting+1
  }
}
Species_richness_TimeSeries_Dataframe<-as.data.frame(Species_richness_TimeSeries_Dataframe)
colnames(Species_richness_TimeSeries_Dataframe)<-c("Scenario","Replicate", "Timestep", "Richness")
Species_richness_TimeSeries_Dataframe$Replicate<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Replicate))
Species_richness_TimeSeries_Dataframe$Timestep<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Timestep))
Species_richness_TimeSeries_Dataframe$Richness<-as.numeric(as.character(Species_richness_TimeSeries_Dataframe$Richness))
Species_richness_TimeSeries3$Means<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, mean)
Species_richness_TimeSeries3$SDs<- tapply(Species_richness_TimeSeries_Dataframe$Richness, Species_richness_TimeSeries_Dataframe$Timestep, sd)
par(mfrow=c(1,1))
plot(Species_richness_TimeSeries3$Means, las=1, type="l",lwd=4,ylim=c(0,19000),xlab="Time step", ylab="Number of species", bty="l")
lines(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs, lwd=1)
lines(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs, lwd=1)

plot(log(Species_richness_TimeSeries3$Means), las=1, type="l",lwd=4,ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Species_richness_TimeSeries3$Means-Species_richness_TimeSeries3$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
lines(log(Species_richness_TimeSeries3$Means+Species_richness_TimeSeries3$SDs), lwd=1)

saveRDS(Species_richness_TimeSeries3, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_LowTraitEvo_OnlySurvivingReps.rds", sep=""))
saveRDS(Species_richness_TimeSeries2, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_LowTraitEvo_AllReps.rds", sep=""))

Species_richness_TimeSeries4<-rbind(Species_richness_TimeSeries4,Species_richness_TimeSeries2) #to bind next results

saveRDS(Species_richness_TimeSeries4, paste(getwd(),"/output/CaseStudy1/Richness_Timeseries_AllScenarios_AllReps.rds", sep=""))








###############################Plotting experiments together:
Richness_Timeseries_Experiment <- readRDS("~/workshops/sELDIG/Genesis/gen3sis_beta/output/CaseStudy1/Richness_Timeseries_LowTraitEvo_OnlySurvivingReps.rds")
plot(log(Richness_Timeseries_Experiment$Means), las=1, type="l",lwd=4,ylim=c(0,10),xlab="Time step", ylab="Number of species (log)", bty="l")
Lowerline<-log(Richness_Timeseries_Experiment$Means-Richness_Timeseries_Experiment$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1)
lines(log(Richness_Timeseries_Experiment$Means+Richness_Timeseries_Experiment$SDs), lwd=1)

Richness_Timeseries_Experiment_2 <- readRDS("~/workshops/sELDIG/Genesis/gen3sis_beta/output/CaseStudy1/Richness_Timeseries_RefTraitEvo_OnlySurvivingReps.rds")
ncol(Richness_Timeseries_Experiment_2)-2 #number of simulations: 27
lines(log(Richness_Timeseries_Experiment_2$Means), lwd=4, lty=2,las=1, type="l",ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Richness_Timeseries_Experiment_2$Means-Richness_Timeseries_Experiment_2$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1, lty=2)
lines(log(Richness_Timeseries_Experiment_2$Means+Richness_Timeseries_Experiment_2$SDs), lwd=1, lty=2)

Richness_Timeseries_Experiment_3 <- readRDS("~/workshops/sELDIG/Genesis/gen3sis_beta/output/CaseStudy1/Richness_Timeseries_HighTraitEvo_OnlySurvivingReps.rds")
ncol(Richness_Timeseries_Experiment_3)-2 #number of simulations: 14
lines(log(Richness_Timeseries_Experiment_3$Means), lwd=4, lty=3,las=1, type="l",ylim=c(0,10),xlab="Time step", ylab="Number of species", bty="l")
Lowerline<-log(Richness_Timeseries_Experiment_3$Means-Richness_Timeseries_Experiment_3$SDs)
Lowerline[which(Lowerline<0)]<-0
lines(Lowerline, lwd=1, lty=3)
lines(log(Richness_Timeseries_Experiment_3$Means+Richness_Timeseries_Experiment_3$SDs), lwd=1, lty=3)

Ltys<-c(1,2,3)
legends<-c("Low trait evolution", "Reference trait evolution", "High trait evolution")
legend(x=-7,y=10.5, lty=Ltys, legend=legends, bty="n")


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





