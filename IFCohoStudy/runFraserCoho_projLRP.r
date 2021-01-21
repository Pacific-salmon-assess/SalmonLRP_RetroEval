library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(tmbstan)
library(here)
library(zoo)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("ProjLRP_Functions.r")
  source("plotFunctions.r")
  source("helperFunctions.r")
}
sourceAll()

# Load TMB models

compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))


# ======================================================================
# Read-in Coho data:  
# =====================================================================
setwd(cohoDir)

CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")


# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)

# Roll up escapements, and get Gen Mean of that
genYrs <- 3
AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, genYrs, gm_mean, fill = NA, align="right"))


# ======================================================================
# Specify initial parameters / priors for TMB estimation:  
# =====================================================================

# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                      logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean_HM<-c(10.957092, 5.565526, 11.467815, 21.104274, 14.803877)

TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                               logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_HM, cap_sig=sqrt(2))

# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean_IM<-c(11.153583,  5.714955, 11.535779, 21.379558, 14.889006)

TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_IM, cap_sig=sqrt(2))


# ===================================================================
# Run Projections
# ==================================================================
 setwd(codeDir)
 devtools::install_github("Pacific-salmon-assess/samSim", ref="LRP")
 
BroodYrLag <- 2
pList <- seq(0.2,1,by=0.2)
year <- 2018
TMB_Inputs <- TMB_Inputs_IM
BMmodel <- "SR_IndivRicker_Surv"
 
 
# Only use SR data for brood years that have recruited by specified year
# (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 2 years) from current year)
SRDat <- CoSRDat %>%  filter(BroodYear <= year-BroodYrLag)
EscpDat.yy <- CoEscpDat %>% filter(yr <= year)

SRDat$yr_num <- group_by(SRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
SRDat$CU_ID <- group_by(SRDat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(SRDat[,c("CU_ID", "CU_Name")]))



# Create samSim input files for current scenario
scenarioName <- "IM.sigER.057"


projSpawners <-run_ScenarioProj(SRDat = SRDat, EscDat = EscDat, BMmodel = BMmodel, scenarioName=scenarioName,
                           useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                          nMCMC=10000, nProj=500, sigER = 0.057, recCorScalar=1)




 scenarioName <- "IM.sigER.057.corr0.75"

 projSpawners <-run_ScenarioProj(SRDat = SRDat, EscDat = EscDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                 useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                 nMCMC=10000, nProj=500, sigER = 0.057, recCorScalar=0.75)


  scenarioName <- "IM.sigER.057.corr0.5"
 
 
  projSpawners <-run_ScenarioProj(SRDat = SRDat, EscDat = EscDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                  useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                  nMCMC=10000, nProj=500, sigER = 0.057, recCorScalar=0.5)
 
 
  
  scenarioName <- "IM.sigER.07125"
  
  
  projSpawners <-run_ScenarioProj(SRDat = SRDat, EscDat = EscDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                  useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                  nMCMC=10000, nProj=500, sigER = 0.07125, recCorScalar=1)
  
 
  
   scenarioName <- "IM.sigER.0855"


 projSpawners <-run_ScenarioProj(SRDat = SRDat, EscDat = EscDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                 useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                 nMCMC=10000, nProj=500, sigER = 0.0855, recCorScalar=1)

 
 
# ===================================================================
# Estimate LRPs
# ==================================================================

# Read in projection outputs to create input lists for logistic regression

OMsToInclude<-c("IM.sigER.057", "IM.sigER.07125" ,"IM.sigER.0855", "IM.sigER.057.corr0.75","IM.sigER.057.corr0.5")



for (i in 1:length(OMsToInclude)) {
  filename<-paste("projLRPDat_",OMsToInclude[i],".csv",sep="")
  dat.i<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
  dat.i<-dat.i %>% filter(year > max(SRDat$yr_num)+4)
  dat.i$OM.Name<-OMsToInclude[i]
  if (i == 1) projLRPDat<-dat.i
  if (i > 1) projLRPDat<-rbind(projLRPDat,dat.i)
  
  filename<-paste( "projSpwnDat_",OMsToInclude[i],".csv",sep="")
  spDat.i<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))
  spDat.i$OM.Name<-OMsToInclude[i]
  if (i == 1) projCUSpDat<-spDat.i
  if (i > 1) projCUSpDat<-rbind(projCUSpDat,spDat.i)
}

 
 
# Calculate LRP and associated probability interval based on distribution of sAg 
  # -- Calculate and save LRPs by OM 

LRPs_byOM<-projLRPDat %>% group_by(OM.Name,ppnCUsLowerBM) %>% 
  summarise(LRP.50=median(sAg), LRP.95=quantile(sAg,0.95),LRP.05=quantile(sAg,0.05))

  # -- Add rows for all OMs combined
LRPs_combined<-projLRPDat %>% group_by(ppnCUsLowerBM) %>% 
  summarise(LRP.50=median(sAg), LRP.95=quantile(sAg,0.95),LRP.05=quantile(sAg,0.05))
LRPs_combined<-LRPs_combined %>% add_column(OM.Name=rep("Combined",nrow(LRPs_combined)), .before="ppnCUsLowerBM")

# Final LRP table
LRPs<-bind_rows(LRPs_byOM, LRPs_combined)


figDir <- here(cohoDir, "Figures", "ProjectedLRPs")
if (file.exists(figDir) == FALSE){
  dir.create(figDir)
}

# Loop over OMs and create / save plots

p<-0.80

for (i in 1:length(OMsToInclude)) {

  # Plot prop-Agg abundance relationship
  projLRPDat.i<-projLRPDat %>% filter (OM.Name == OMsToInclude[i])

  projDat4Plot <- data.frame(AggSpawners=projLRPDat.i[,"sAg"], ppnCUs=projLRPDat.i[, "ppnCUsLowerBM"])
  LRP4Plot <- as.data.frame(LRPs %>% filter(ppnCUsLowerBM==p & OM.Name==OMsToInclude[i]) %>% select(fit=LRP.50, lwr=LRP.05, upr=LRP.95))

  plotProjected(Data = projDat4Plot, LRP = LRP4Plot,
              plotName = paste("ProjMod", OMsToInclude[i], year, sep ="_"), 
              outDir = figDir, p = p)
  
  # Plot spawner abundance projections, by CU
  
  makeSpawnerPlot<- function(i, projSpwnDat, CUNames) {
    
    plotDat<-projSpwnDat %>% filter(CU==i) %>% group_by(year, expRate) %>% 
      summarise(medSpawn=median(spawners), lwr=quantile(spawners,0.10),upr=quantile(spawners,0.90))
    
    p <- ggplot(data=plotDat, mapping=aes(x=year,y=medSpawn, colour=factor(expRate))) +
      geom_ribbon(data=plotDat, aes(ymin = lwr, ymax = upr, x=year, fill=factor(expRate)), alpha=0.2) +
      geom_line(mapping=aes(x=year, y=medSpawn)) +
      geom_line(data=plotDat %>% filter(year < 18), col="black", size=1) +
      ggtitle(CUNames[i]) +
      xlab("Year") + ylab("Spawners") +
      theme_classic()  
    
  }
  
  projCUSpDat.i<- projCUSpDat %>% filter(OM.Name == OMsToInclude[i])
  
  ps<-lapply(1:length(unique(SRDat$CU_Name)), makeSpawnerPlot, projSpwnDat = projCUSpDat.i,CUNames=unique(SRDat$CU_Name))
  
  pdf(paste(cohoDir,"/Figures/ProjectedLRPs/", OMsToInclude[i], "_CUSpawnerProj.pdf", sep=""), 
      width=9, height=6)
  do.call(grid.arrange,  ps)
  dev.off()

}





#### Plot to compare LRPs


#OMsToInclude<-c("IM.sigER.057", "IM.sigER.0855", "IM.sigER.057.corr0.75","IM.sigER.057.corr0.5")


plotDat2<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% c("IM.sigER.057", "IM.sigER.057.corr0.75", "IM.sigER.057.corr0.5" ))

g <-ggplot(data=plotDat2, mapping=aes(x=OM.Name, y=LRP.50)) +
             geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
             xlab("Operating Model") + ylab("LRP; p = 0.8") +
             scale_x_discrete(labels=c("Base", "EscCorr75%", "EscCorr50%"))
  
  
plotDat3<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% c("IM.sigER.057", "IM.sigER.07125", "IM.sigER.0855"))

g2 <-ggplot(data=plotDat3, mapping=aes(x=OM.Name, y=LRP.50)) +
  geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
  xlab("Operating Model") + ylab("LRP; p = 0.8") +
  scale_x_discrete(labels=c("Base", "ERSig125%", "ERSig150%"))


