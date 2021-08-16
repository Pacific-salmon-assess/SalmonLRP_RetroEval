# ============================================================================
# Calculation of Projected Limit Reference Points for Interior Fraser Coho
# Carrie Holt, Last update: May 12, 2021
# ================================================================================

# Projected LRPs represent the aggregate abundance that is associated with a specified probability
# that a required proportion of CUs will be above their lower benchmarks in projections
#--  E.g., the LRP may represent the aggregate abundance that is projected to have a 50% probability
# that all CUs (100% of CUs) will be above their lower benchmarks

# Projections are done using the LRP branch of the samSim model (https://github.com/Pacific-salmon-assess/samSim)
# The code in this file is divided into the following sections:
#     (1) Read-in Coho data
#     (2) Specify initial parameters & data sets for projections
#     (3) Run base projections (Using 4 different OM models at present)
#     (4) Run sensitivity analysis projections
#     (5) Estimate and save LRPs, and associated plots
#     (6) Plot CU-level spawner abundance projections (Optional)
#     (7) Code to create mcmcOut for Ricker pars from an assumed distn of a and
#         SREP
#     (8) Code to calculate tau for variability in age proportions
#     (9) Code to plot distribution of correlations among CUs/inlets
#     (10) Make histograms of cvER
#     (11) Plots of LRP stabilitization with number of trials
#     (12) Plot LRPs with various plevels
#     (13) Run reconstruction for WCVI CK CUs
#     (14) Plot age proportions in recruitment
#     (15) Make comparison plots among scenarios (NOT CURRENTLY WORKING)

# ===============================================================================


library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(viridis)
#library(tmbstan)
library(here)
library(zoo)
library(corrplot)


setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
wcviCKDir<-paste(rootDir,"/WCVIChinookStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("ProjLRP_Functions.r")
  source("plotFunctions.r")
  source("helperFunctions.r")
  source("get.mv.logistic.tau.r")
}
sourceAll()

# # Load TMB models for fitting Bayesian Ricker stock recruit models;
# #   outputs from these model fits are used to parameterize samSim
#
# compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))
#
# compile("TMB_Files/SR_HierRicker_Surv_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_noLRP"))
#
# compile("TMB_Files/SR_IndivRicker_SurvCap_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_noLRP"))
#
# compile("TMB_Files/SR_HierRicker_SurvCap_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_noLRP"))
#

# ======================================================================
#(1)  Read-in WCVI Chinook data:
# =====================================================================
setwd(wcviCKDir)


wcviCKSRDat <- read.csv("DataIn/Inlet_Sum.csv")
wcviCKSRDat$yr_num <- group_by(wcviCKSRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
wcviCKSRDat$CU_ID <- group_by(wcviCKSRDat, Inlet_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
# ======================================================================
# (2) Specify initial parameters & datasets for projections
# =====================================================================

# # Subset data up to current year
# year <- 2018 # - last year of data for parameterization
# BroodYrLag <- 2 # - number of years between brood year and first age of return (2 years for coho)
#
# # Only use SR data for brood years that have recruited by specified year
# # (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 2 years) from current year)
# SRDat <- CoSRDat %>%  filter(BroodYear <= year-BroodYrLag)
#
# SRDat$yr_num <- group_by(SRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
# SRDat$CU_ID <- group_by(SRDat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
#
#
# # TMB input parameters:
# TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1,
#                       logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1,
#                       gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)
#
#
# TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
#                       Tau_dist = 0.1,
#                       gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)
#
#
# # Prior means come from running "compareRickerModelTypes.r"
# cap_priorMean_HM<-c(10.957092, 5.565526, 11.467815, 21.104274, 14.803877)
#
# TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1,
#                                logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1,
#                                gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
#                                cap_mean=cap_priorMean_HM, cap_sig=sqrt(2))
#
# # Prior means come from running "compareRickerModelTypes.r"
# cap_priorMean_IM<-c(11.153583,  5.714955, 11.535779, 21.379558, 14.889006)
#
# TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1,
#                                gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
#                                cap_mean=cap_priorMean_IM, cap_sig=sqrt(2))



# Create output directories for Projected LRP outputs
figDir <- here(wcviCKDir, "Figures")
if (file.exists(figDir) == FALSE){
  dir.create(figDir)
}
figDir2 <- here(figDir, "ProjectedLRPs")
if (file.exists(figDir2) == FALSE){
  dir.create(figDir2)
}

projOutDir <- here(wcviCKDir, "DataOut")
if (file.exists(projOutDir) == FALSE){
  dir.create(projOutDir)
}
projOutDir2 <- here(projOutDir, "ProjectedLRPs")
if (file.exists(projOutDir2) == FALSE){
  dir.create(projOutDir2)
}


# ===================================================================
# (3) Run Base Projections
# ==================================================================

setwd(codeDir)
devtools::install_github("Pacific-salmon-assess/samSim", ref="LRP")




# Create a correlation matrix from spawner time-series, as a proxy for
# correlation in recruitment residuals assuming no density dependence and
# constant harvest. Only used if recruitment time-series are missing
dum <- wcviCKSRDat %>% dplyr::select(CU_ID, BroodYear, Spawners)
dum <- dum %>% pivot_wider(id_cols=c(CU_ID, BroodYear), names_from=CU_ID,
                           values_from=Spawners) %>% dplyr::select (!BroodYear)
dum <- dum %>% drop_na()
corMat <- cor(dum)


# Plot Bubble plot of correlations
# rownames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))$stkName
# colnames(corMat) <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))$stkName
#
# # png(filename=paste(wcviCKDir, "/Figures/SpawnerCorrelation.png", sep=""), width=4, height=4.5, units="in", res=500)
# corrplot(corMat, method="circle", p.mat=corMat, insig="p-value", type="lower")
# # dev.off()

# Alternatively, create correlation from Ricker residuals from run
# reconstruction. However, reconstruction currently done at CU level. Although I
# have inlet-level raw escapements, without infilling gaps, the time-series of
# comparable data ends up being n=5 years, which is too short for a cor matrix
# (gives very high & variable correlations), and results in not pos. def. errors
# Ricker residuals estimated in
# "WCVI_term_model_revisions_updated_withIneltReconstruction.xlsx"
#wcviRicResids <- data.frame(read.csv("DataIn/WCVI_RickerResids.csv"))
#corMat <- cor(na.omit(wcviRicResids))

# InletRickerResid <- data.frame(read.csv(paste(wcviCKDir, "DataIn/CURickerResid.csv",sep="/")))
# acf(InletRickerResid$SWVI[7:30], na.rm=TRUE)
# acf(InletRickerResid$No.KY, na.rm=TRUE)
# acf(InletRickerResid$NWVI, na.rm=TRUE)
# png(paste(wcviCKDir,"/Figures/acfWCVIRickerResid.png", sep=""), width=5, height=4, units="in", res=500)
# acf(InletRickerResid$WCVI)
# dev.off()
#-------------------------------------------------------------------------------
#TESTING
# SRDat <- SRDat %>% mutate(Recruits=NA) %>% dplyr::select(-c('Age_3_Recruits',
#                                                      'Age_4_Recruits',
#                                                      'STAS_Age_3', 'STAS_Age_4',
#                                                      'ER_Age_3', 'ER_Age_4',
#                                                      'Hatchery'))

#-------------------------------------------------------------------------------

# Create samSim input files for current scenario
setwd(codeDir)






scenarioName <- "baseER"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)

scenarioName <- "baseERn10000"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)

scenarioName <- "cvER0" # should also be run with 50,000 for LRPs

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000,cvER = 0, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)

scenarioName <- "cvER0.17" #should also be run with 50,000 for LRPs

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.17, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)


scenarioName <- "ER0"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0)

scenarioName <- "ER0.05"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.05)
scenarioName <- "ER0.10"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.10)

scenarioName <- "ER0.15"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.15)

scenarioName <- "ER0.2"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.2)

scenarioName <- "ER0.25"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.25)



scenarioName <- "ER0.35"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.35)

scenarioName <- "ER0.4"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.4)

scenarioName <- "ER0.45"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.45)





#
scenarioName <- "recCorSca0"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=0, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
# scenarioName <- "recCorSca0.1"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.1, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
#
# scenarioName <- "recCorSca0.2"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.2, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
# scenarioName <- "recCorSca0.3"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.3, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
#
# scenarioName <- "recCorSca0.4"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.4, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
scenarioName <- "recCorSca0.5"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=0.5, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
# scenarioName <- "recCorSca0.6"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.6, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
# scenarioName <- "recCorSca0.7"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.7, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
# scenarioName <- "recCorSca0.8"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.8, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)
#
# scenarioName <- "recCorSca0.9"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=0.9, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)


scenarioName <- "cvER0n50000" # should also be run with 50,000 for LRPs

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000,cvER = 0, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)

scenarioName <- "cvER0.17n50000" #should also be run with 50,000 for LRPs

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.17, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)

scenarioName <- "agePpnConst"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=TRUE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3)



scenarioName <- "alphaScalar0.75n50000"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3,
                                alphaScalar=0.75, SREPScalar=1)


scenarioName <- "alphaScalar1.5n50000"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3,
                                alphaScalar=1.5, SREPScalar=1)





scenarioName <- "Anarrow"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3,
                                aNarrow=TRUE)

scenarioName <- "AlifeStageModel"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3,
                                alphaScalar="lifeStageModel", SREPScalar=1)


scenarioName <- "annualcvERCU"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=TRUE, biasCorrectProj=TRUE, ER=0.3)

# scenarioName <- "noMCMC"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=F,
#                                 nMCMC=NULL, nProj=3000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.29)
# ==================================================================
# (4) Run Sensitivity Analyses
# ====================================================================

#
# # Create samSim input files for current scenario
# scenarioName <- "IM.cvER1.5"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName, useGenMean = F,
#                                 genYrs = genYrs,  TMB_Inputs=NULL,
#                                 outDir=wcviCKDir, runMCMC=F, nMCMC=NA,
#                                 nProj=10, cvER = 0.456*1.5, recCorScalar=1)
#
#
# # Create samSim input files for current scenario
# scenarioName <- "IM.cvER2.0"
# #etc. I should include variability in input Ricker a, b, sig, covar, as well as cvER


# ===================================================================
#  (5) Estimate and Save LRPs (and associated plots)
# ==================================================================

# Specify threshold to use when calculating LRP
# # Note: may want to loop over probThresholds as well; still needs to be added
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-0.50 # probability theshhold; the LRP is set as the aggregate abundance that has this
# probability that the propCUThreshold is met

# Specify scenarios to calculate LRPs and make plots for.
# These scenarios will be looped over below with a LRP (and LRP plot) saved for each scenario
OMsToInclude<-c(
  # "baseERn10000",
  # "agePpnConst")
  "cvER0",
  "baseERn10000",
  "cvER0.17")
# "cvER0.21ER0",
  # "cvER0.21ER0.43")
# #"cvER0.21n50000_20yrs")
  # "cvER0.21.AlifeStageModel")
  # "cvER0.21.annualcvERCU")
  # "cvER0.42")
# "cvER0.21",
# "cvER0.21.recCorSca0",
# "cvER0.21.recCorSca0.1",
# "cvER0.21.recCorSca0.2",
# "cvER0.21.recCorSca0.3",
# "cvER0.21.recCorSca0.4",
# "cvER0.21.recCorSca0.5",
# "cvER0.21.recCorSca0.6",
# "cvER0.21.recCorSca0.7",
# "cvER0.21.recCorSca0.8",
# "cvER0.21.recCorSca0.9")
# "cvER0.21",
# "cvER0.21.agePpnConst")

# "cvER0",
# "cvER0.21",
# "cvER0.42",
# "cvER0.21.agePpnConst")
  # "cvER0.21",
  # "cvER0.21.alphaScalar0.75",
  # "cvER0.21.alphaScalar1.5",
  # "cvER0.21.Anarrow",
  # "cvER0.21.noMCMC")
  #


# Loop over OM Scenarios
for (i in 1:length(OMsToInclude)) {

  # Read in samSim outputs for OM
  filename<-paste("projLRPDat_",OMsToInclude[i],".csv",sep="")
  projLRPDat<-read.csv(here(wcviCKDir, "SamSimOutputs", "simData",filename))
  CUpars <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))
  projLRPDat<-projLRPDat %>% filter(year > CUpars$ageMaxRec[1]*10)#)max(SRDat$yr_num)+4)

  # Create bins for projected spawner abundances
  minBreak<-0
  maxBreak<-round(max(projLRPDat$sAg),digits=-2)
  binSize<-200 # Note: bin size is currently set here
  breaks<-seq(minBreak, maxBreak,by=binSize)

  # Set bin labels as the mid-point
  projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))

  # Summarize nSims in each bin
  tmp<-projLRPDat  %>%  group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))

  # Filter out bins with < 10 nSims
  tmp2<-projLRPDat %>% group_by(bins)%>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
    add_column(nSims=tmp$nSims) %>% filter(nSims>=100)

  # For each bin, calculate probability that required proportion of CUs above benchmark
  projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
  # For each bin, calculate the difference between the threshold probability and the calculated probability
  projLRPDat$diff<-abs(probThresh-projLRPDat$prob)

  # Save projection summaries used to create plots
  projLRPDat$OM.Name<-OMsToInclude[i]
  if (i == 1) projLRPDat.plot<-projLRPDat
  if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)

  # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
  LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))

  # Create a table of LRP estimates to be saved for each OM model
  if (i ==1) {
    LRP_Ests<-data.frame(OMsToInclude[i], probThresh, propCUThresh, LRP, binSize)
    names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
  } else {
    tmp.df<-data.frame(OMsToInclude[i], probThresh, propCUThresh, LRP, binSize)
    names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
    LRP_Ests<-rbind(LRP_Ests,tmp.df)
  }

  # Plot projected LRP abundance relationship ===============================================================
  png(paste(wcviCKDir,"/Figures/ProjectedLRPs/", OMsToInclude[i],
            "_ProjLRPCurve_prob",probThresh,".png", sep=""), width=5, height=4,
      units="in", res=500)
  # pdf(paste(wcviCKDir,"/Figures/ProjectedLRPs/", OMsToInclude[i], "_ProjLRPCurve_prob",probThresh,".pdf", sep=""),
  #     width=6, height=6)

  plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19,
       xlim=c(0, max( as.numeric(as.character(projLRPDat$bins)), na.rm=T) ),
       cex=0.5, cex.lab=1.5,
       xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
  abline(h=probThresh, lty=2)
  abline(v=LRP, lwd=4, col="orange")

  dev.off()

  # Option to plot histogram of nSims in each Agg Abundance Bin
  #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)

}


# Save LRPs for all OM scenarios
write.csv(LRP_Ests, paste(projOutDir2, "ProjectedLRPs_cvER.csv", sep="/"), row.names=F)
# Save LRP projection summaries used for calculating and plotting LRP (Optional)
write.csv(projLRPDat.plot, paste(projOutDir2, "ProjectedLRP_data_cvER.csv", sep="/"), row.names=F)





# ===================================================================
# (6) Plot CU-level Spawner Abundance Projections (Optional)
# ==================================================================


for (i in 1:length(OMsToInclude)) {


  filename<-paste( "projSpwnDat_",OMsToInclude[i],".csv",sep="")
  spDat.i<-read.csv(here(wcviCKDir,"SamSimOutputs", "simData",filename))
  spDat.i$OM.Name<-OMsToInclude[i]
  # if (i == 1) projCUSpDat<-spDat.i
  # if (i > 1) projCUSpDat<-rbind(projCUSpDat,spDat.i)

  projCUSpDat.i<-spDat.i

  # Function to plot spawner abundance projections, by CU
  makeCUSpawnerProjPlot<- function(i, projSpwnDat, CUNames) {
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

  ps<-lapply(1:length(unique(cuPar$stkName)), makeCUSpawnerProjPlot, projSpwnDat = projCUSpDat.i,CUNames=unique(cuPar$stkName))

  pdf(paste(wcviCKDir,"/Figures/ProjectedLRPs/", OMsToInclude[i], "_CUSpawnerProj.pdf", sep=""),
      width=9, height=6)
  do.call(grid.arrange,  ps)
  dev.off()

}




# ===================================================================
# (7) Code to create mcmcOut for SR parameters for WCVI CK
# ===================================================================

# The mcmc is written to a file that is read-in by ProjRLRP_Functions.r

# SREP files are from github repository, Watershed-Area-Model and are from
# "Watershed-Area-Model/DataOut/WCVI_SMSY_noEnh.csv"

# # Hard coding from Watershed-Area-Model directory
# if (remove.EnhStocks) SREP <- data.frame(read.csv(
#   "c:/github/Watershed-Area-Model/DataOut/WCVI_SMSY_noEnh.csv"))
# if (!remove.EnhStocks) SREP <- data.frame(read.csv(
#   "c:/github/Watershed-Area-Model/DataOut/WCVI_SMSY_wEnh.csv"))

# For now, I have copied the SREP files to the SalmonLRP_RetroEval repository
# If the watershed-area-model is updated, these files will need to be updated

createMCMCout <- TRUE
setwd(wcviCKDir)
alphaScalar <- 1
SREPScalar <- 1
# Only need to run once to create mcmcOut.csv file with a given assumed
# distribution of alpha and SREP
Inlet_Names <- read.csv(paste("samSimInputs/CUPars.csv"))$stkName
CU_inlet <- data.frame(Inlet_Names=Inlet_Names, CU_Names=NA)
CU_inlet[Inlet_Names=="Barkley",2] <- "WCVI South"#"Southwest_Vancouver_Island"
CU_inlet[Inlet_Names=="Clayoquot",2] <- "WCVI South"#"Southwest_Vancouver_Island"
CU_inlet[Inlet_Names=="Kyuquot",2] <- "WCVI Nootka & Kyuquot"#"Nootka_Kyuquot"
CU_inlet[Inlet_Names=="Nootka/Esperanza",2] <- "WCVI Nootka & Kyuquot"#"Nootka_Kyuquot"
CU_inlet[Inlet_Names=="Quatsino",2] <- "WCVI North"#"Northwest_Vancouver_Island"


if(createMCMCout){
  set.seed(1)
  remove.EnhStocks <- TRUE
  nTrials <- 50000
  # Set up matrix of random numbers to use for generating alphas, so that
  # the same random numbers are used for Ricka estimates with bias correction
  # and without bias correction when using alpha to estimata beta (lnA/SREP)
  a_rand <- matrix(runif(nTrials*1.5*length(Inlet_Names)), nrow=nTrials*1.5, ncol=length(Inlet_Names))

  # Pull SREP estimates from Watershed-Area model (see repository "Watershed-
  # Area-Model"). That model included a bias correction for back-transformation
  # from log-space
  if (remove.EnhStocks) SREP <- data.frame(read.csv(
    "DataIn/WCVI_SMSY_noEnh_wBC.csv"))
  if (!remove.EnhStocks) SREP <- data.frame(read.csv(
    "DataIn/WCVI_SMSY_wEnh_wBC.csv"))



  #Get lnaplha

  # lnalpha from Diana Dobson's run reconstruction coded in R/TMB with bias correction
  lnalpha_inlet <- read.csv("samSimInputs/CUPars.csv") %>% dplyr::select(alpha,stkName) %>% rename(inlets=stkName, lnalpha=alpha)
  lnalpha_nBC_inlet <- read.csv("samSimInputs/CUPars_nBC.csv") %>% dplyr::select(alpha,stkName) %>% rename(inlets=stkName, lnalpha_nBC=alpha)
  lnalpha_inlet$lnalpha <- lnalpha_inlet$lnalpha * alphaScalar#1 #(value for life-stage model =1)
  lnalpha_nBC_inlet$lnalpha_nBC <- lnalpha_nBC_inlet$lnalpha_nBC * alphaScalar#1-0.5^2/2#

  Inlet_Names <- lnalpha_inlet$inlets

  #Inlet_Names <- c("Kyuquot", "Clayoquot", "Quatsino", "Barkley", "Nootka/Esperanza")
  SREP <- SREP %>% filter(Stock %in% Inlet_Names) %>% filter(Param=="SREP") %>%
    dplyr::select(!c(X, Param)) %>% rename(SREP=Estimate, inlets=Stock)
  SREP <- SREP %>% mutate(SREP=SREP * SREPScalar) %>%
    mutate(LL=LL * SREPScalar) %>%
    mutate(UL=UL * SREPScalar)


  out <- SREP %>% left_join(lnalpha_inlet, by="inlets") %>% left_join(lnalpha_nBC_inlet, by="inlets")

  #Draw alpha value, then draw logSREP parameters,then calc beta for that draw
  # (lnalpha/SREP)
  for (i in 1:length(Inlet_Names)){
    meanSREP <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(SREP)
    logmeanSREP <- log(meanSREP)
    ULSREP <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(UL)
    logULSREP <- log(ULSREP)
    LLSREP <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(LL)
    logLLSREP <- log(LLSREP)
    sigSREP <- (logmeanSREP-logLLSREP)/1.96
    #sigSREP <- (logULSREP-logmeanSREP)/1.96 #Check should be same
    rSREP <- exp(rnorm(nTrials*1.5, logmeanSREP,sigSREP))

    rsig <-  rsig <- read.csv(paste("samSimInputs/CUPars.csv")) %>%
      filter(stkName==Inlet_Names[i]) %>% dplyr::select(sigma,stk)

    meanlnalpha_nBC <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(lnalpha_nBC)#1-(rsig$sigma^2)/2# (life-stage model)
    meanlnalpha <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(lnalpha)#1# (life-stage model)
    # ULlnalpha <- 2
    # LLlnalpha <- 0
    siglnalpha <- 0.5 # Assuming 95% CIs at 0 and 2, sig ~0.5.#0.25 (narrow)


    # Generate random lnalpha values using same random numbers with and withtout
    # bias correction (but diff for each CU or inlet)
    rlnalpha_nBC <- data.frame(a=qnorm(a_rand[,i], meanlnalpha_nBC, siglnalpha))
    rlnalpha <- data.frame(a=qnorm(a_rand[,i], meanlnalpha, siglnalpha))
    amin <- 0#(meanlnalpha - siglnalpha)# (narrow)- not implemented
    amax <- max(2,alphaScalar*2)#(meanlnalpha + siglnalpha)# (narrow)- not implemented



    # Create a dataframe of alpha (with BC), beta (from alpha w/out BC to
    #stabilize beta with and without BC)
    df <- data.frame( stk=rsig$stk, alpha=rlnalpha$a,
                      beta=rlnalpha_nBC$a/rSREP, sigma= rsig$sigma, SREP=rSREP,
                      stkName=Inlet_Names[i], alpha_nBC = rlnalpha_nBC$a )
    #Remove all rows with Ricker a greater or less than bounds
    df <- df %>% filter(alpha > amin & alpha < amax & alpha_nBC > amin &
                          alpha_nBC < amax) %>% slice(1:nTrials)

    if (i==1) mcmcOut <- df
    if (i>1) mcmcOut <- mcmcOut %>% add_row(df)

  }

  if(alphaScalar==1&SREPScalar==1) write.csv(mcmcOut, paste(wcviCKDir, "SamSimInputs","Ricker_mcmc.csv", sep="/"),#"Ricker_mcmc_narrow.csv",#_lifeStageModel
                          row.names=F)

  if(alphaScalar!=1 | SREPScalar!=1) write.csv(mcmcOut, paste(wcviCKDir, "/SamSimInputs/Ricker_mcmc_alphaScalar",alphaScalar ,"_SREPScalar",SREPScalar,".csv", sep=""),
            row.names=F)

  #plot of alpha and SREP density
  alphaDensity <- mcmcOut %>% ggplot(aes(alpha, colour=factor(stkName))) + #geom_histogram()
    geom_density() +theme(legend.title = element_blank())

  # Histogram might be better to see smoothed curves with 50,000 mcmc
  mcmcOut %>% ggplot(aes(alpha)) + geom_histogram() + facet_wrap(~factor(stkName))

  SREPDensity <- mcmcOut  %>%
    ggplot(aes(SREP, colour=factor(stkName), fill=factor(stkName))) +
    geom_density(alpha=0.1) +theme(legend.title = element_blank()) +xlim(0,30000)

  if(alphaScalar==1 & SREPScalar==1) {
    ggsave(paste(wcviCKDir,"/Figures/AlphaDensity.png",sep=""),#_lifeStageModel#_narrow
           plot = alphaDensity,
           width = 6, height = 4, units = "in")
    ggsave(paste(wcviCKDir,"/Figures/SREPDensity.png",sep=""),
           plot = SREPDensity,#"/Figures/AlphaDensity_narrow.png"
           width = 6, height = 4, units = "in")
  }

  if(alphaScalar!=1)  {
    ggsave(paste(wcviCKDir,"/Figures/AlphaDensity_alphaScalar",alphaScalar,
                 ".png", sep=""),
           plot = alphaDensity,
           width = 6, height = 4, units = "in")
  }
  if(SREPScalar!=1) {
    ggsave(paste(wcviCKDir,"/Figures/SREPDensity_SREPScalar",SREPScalar,
                 ".png", sep=""),
           plot = SREPDensity,
           width = 6, height = 4, units = "in")

  }

}


#Look at Sgens by inlet
# test <- mcmcOut#
 ##test <- read.table("samSimInputs/Ricker_mcmc_alphaScalar0.75_SREPScalar1.csv")
# SGEN <- unlist(purrr::map2 (.x=test$alpha, .y=test$beta, .f=sGenSolver))
# SGEN <- data.frame(SGEN=SGEN, stkName=test$stkName)
# sGenDensity <- SGEN %>% ggplot(aes(SGEN)) + geom_density() +
#   xlim(0,3000) + facet_wrap(~stkName)
# ggsave(paste(wcviCKDir,"/Figures/sGenDensity.png", sep=""),
#        plot=sGenDensity, width=4, height=3, units="in")



# ===================================================================
# (8) Code to estimate uncertainty in age ppns in recruitment by BY
# ===================================================================

calcTau <- FALSE
if(calcTau){
  setwd(wcviCKDir)

  Inlet_Names <- read.csv(paste("samSimInputs/CUPars.csv"))$stkName
  CU_inlet <- data.frame(Inlet_Names=Inlet_Names, CU_Names=NA)
  CU_inlet[Inlet_Names=="Barkley",2] <- "Southwest_Vancouver_Island"
  CU_inlet[Inlet_Names=="Clayoquot",2] <- "Southwest_Vancouver_Island"
  CU_inlet[Inlet_Names=="Kyuquot",2] <- "Nootka_Kyuquot"
  CU_inlet[Inlet_Names=="Nootka/Esperanza",2] <- "Nootka_Kyuquot"
  CU_inlet[Inlet_Names=="Quatsino",2] <- "Northwest_Vancouver_Island"

  CUages <- data.frame(read.csv("DataIn/CUages.csv"))
  CU.tau <- NA
  for (i in 1: length(unique(CUages$CU_Names))){
    CUages.byCU <- CUages %>% filter(CU_Names== unique(CUages$CU_Names)[i]) %>% dplyr::select(-c(Year, CU, CU_Names))
    CU.tau[i] <- "get.mv.logistic.tau"(CUages.byCU)$best.tau
  }
  df <- data.frame(CU_Names=unique(CUages$CU_Names), CU.tau=CU.tau)
  inletTau <- left_join(CU_inlet, df)

  # Use these tau values for `tauCycAge` in samSim- logistic variation in age
  #  structure
  # Inlet_Names                   CU_Names CU.tau
  #           Kyuquot             Nootka_Kyuquot    0.6
  #         Clayoquot Southwest_Vancouver_Island    0.7
  #          Quatsino Northwest_Vancouver_Island    0.7
  #           Barkley Southwest_Vancouver_Island    0.7
  #  Nootka/Esperanza             Nootka_Kyuquot    0.6

}



# ===================================================================
# (9) Code to plot distribution of correlations among CUs/inlets
# ===================================================================


OMsToTest<-c(
  "recCorSca0",
  "recCorSca0.5",
  "baseERn10000")
  # "baseERn10000",
  # "agePpnConst")
    # "cvER0",
    # "baseERn10000",
    # "cvER0.17")

LRPFileName <- "ProjectedLRPs_cvERage.csv"

for (j in 1:length(OMsToTest)) {

  filename<-paste( "projSpwnDat_",OMsToTest[j],".csv",sep="")
  spDat<-read.csv(here(wcviCKDir,"SamSimOutputs", "simData",filename))

  spDat<-as_tibble(spDat)
  spDat<-spDat %>% dplyr::select(-X)

  RecCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))
  SpwnCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))

  for (i in 1:max(spDat$iteration)) {

    recruits.i<-spDat %>% filter(iteration==i) %>%  dplyr::select(-spawners)# & expRate==0.325
    recruits.i<-recruits.i %>% dplyr::select(-expRate, -iteration)
    cor_mat<-recruits.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=recruits) %>% dplyr::select(-year) %>% cor()
    RecCorMat[,,i]<-cor_mat

    spawners.i<-spDat %>% filter(iteration==i) %>%  dplyr::select(-recruits)# & expRate==0.325
    spawners.i<-spawners.i %>% dplyr::select(-expRate, -iteration) %>% filter(year>40)
    cor_mat<-spawners.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=spawners) %>% dplyr::select(-year) %>% cor()
    SpwnCorMat[,,i]<-cor_mat

    if (i ==1) SpwnCorrValues<-SpwnCorMat[,,i][lower.tri(SpwnCorMat[,,i])==TRUE]
    if (i > 1) SpwnCorrValues<-c(SpwnCorrValues,SpwnCorMat[,,i][lower.tri(SpwnCorMat[,,i])==TRUE])
  }

  OM_Name<-rep(OMsToTest[j],max(spDat$iteration))

  if (j ==1) SpwnCorr.df<-data.frame(OM_Name,SpwnCorrValues)
  if (j > 1) {
    tmp<-data.frame(OM_Name,SpwnCorrValues)
    SpwnCorr.df<-rbind(SpwnCorr.df,tmp)
  }

}

# Identify a subset of correlations to jitter on the boxplot
# SpwnCorr.df <- read.csv(paste(wcviCKDir, "/DataOut/ProjectedLRPs/SpwnCorr.df.csv", sep=""))
len <- dim(SpwnCorr.df)[1]
SpwnCorr.df <- SpwnCorr.df %>% dplyr::select(c(OM_Name,SpwnCorrValues))
SpwnCorr.df <- SpwnCorr.df %>% add_column(
  highlight = sample(c(1,NA), len, replace=T,prob=c(0.005*len,0.995*len)))
# SpwnCorr.df <- SpwnCorr.df %>% mutate(
#   highlight = c(sample(c(1,NA), len, replace=T,prob=c(0.01*len,0.99*len)),rep(1,10)))


# Add observed escapement correlations to correlation data frame
# cor_mat <- corMat
# SpwnCorrValues.Obs <- cor_mat[lower.tri(cor_mat)==TRUE]

# Over the most recent 20 years:
Years <- 1995:2000
tri <- matrix(NA, nrow=10,ncol=length(Years))
med <- NA
Inlet_Sum.df.long <- read.csv(paste(wcviCKDir, "/DataIn/Inlet_Sum.csv", sep=""))
for (i in 1:length(Years)){
  Year.Start <- Years[i]
  Inlet_Sum.df <-   Inlet_Sum.df.long %>%
    filter(BroodYear %in% c(Year.Start : (Year.Start + 20) ) ) %>%
    dplyr::select(-c(CU_Name, Inlet_ID, Recruits)) %>%
    pivot_wider(id_cols= c (Inlet_Name, BroodYear),
                names_from=Inlet_Name, values_from=Spawners) %>%
    dplyr::select(-BroodYear) %>% na.omit()

  co <- cor(Inlet_Sum.df)
  tri[,i] <- co[lower.tri(co)==TRUE]
  med[i] <- median(co[lower.tri(co)==TRUE])
}

dum <- as.data.frame(tri)
names(dum) <- as.character(Years)
tri_longer <- dum %>% pivot_longer(col=names(dum), names_to="StartYear",
                                   values_to="Correlation")
tri_longer$StartYear <- factor(tri_longer$StartYear, levels = names(dum),
                               ordered = TRUE)

SpwnCorrValues.Obs <- tri_longer %>% filter(StartYear==2000) %>% pull(Correlation)


#Add observed values to simulated values in df
tmp<-data.frame(OM_Name = "Observed",SpwnCorrValues = SpwnCorrValues.Obs, highlight=1)
SpwnCorr.df<-rbind(SpwnCorr.df,tmp)

if(len > 12000) {ind <- is.na(SpwnCorr.df$highlight)} else ind <- rep(FALSE,len)

# Save LRPs for all OM scenarios
write.csv(SpwnCorr.df, paste(projOutDir2, "SpwnCorr.df.csv", sep="/"), row.names=F)
#SpwnCorr.df <- read.csv(paste(projOutDir2, "SpwnCorr.df.csv", sep="/"))


#LRPs <- read.csv(paste(wcviCKDir,"/DataOut/ProjectedLRPs/", LRPFileName, sep="")) %>% pull(LRP)


factor(SpwnCorr.df$OM_Name,levels = c(
  "recCorSca0",
  "recCorSca0.5",
  # "baseER"),
  "baseERn10000"),
  # "baseERn10000",
  # "agePpnConst"),
  # "cvER0",
  # "baseERn10000",
  # "cvER0.17"),

 ordered = TRUE)


g <- ggplot(SpwnCorr.df,aes(y=SpwnCorrValues,x=as.factor(OM_Name))) +
  geom_boxplot(width=0.5, outlier.shape=NA) + ylim(-0.25,1)+
  geom_jitter(data=SpwnCorr.df[!ind, ],
              position=position_jitter(0.2), col="dark grey", alpha=0.95,
              size=0.2) +
  scale_x_discrete(limits=c("Observed",
                            # "cvER0",
                            # "cvER0.21baseERn10000",
                            # # "cvER0.21.annualcvERCU")
                            # "cvER0.42"),
                            # "baseERn10000",
                            # "agePpnConst"),
                            "recCorSca0",
                            "recCorSca0.5",
                            "baseERn10000"),

                            # "cvER0",
                            # "baseERn10000",
                            # "cvER0.17"),

                       labels=c("Observed",
                            # "0",
                            # "0.085",
                            # #"const\ndeviations\nover years",
                            # # #"annual\ndeviations\nover years")) +
                            # "0.17" )) +
                            "0",
                            "0.5",
                            "1")) +
                            # "Variable\nAge Ppn\nAmong Inlets",
                            # "Constant\nAge Ppn\nAmong Inlets" )) +



  # xlab("cv in Exploitation Rates among inlets") +
  # xlab("Scalar for Ricker Residual Correlation Matrix") +
  xlab("") +
  ylab("Pairwise correlations in spawners among inlets") +
  theme(axis.text=element_text(size=7),  axis.title=element_text(size=8)) #+
  # annotate("text", x = c(1.5,2:4), y = -0.2, label =c("LRP:",LRPs), size=3)
  # annotate("text", x = c(1.4,2:12), y = -0.2, label =c("LRP:",LRPs), size=2)
  # annotate("text", x = c(1.5,2:3), y = -0.2, label =c("LRP:",LRPs), size=5)
  # annotate("text", x = c(1.5,2:5), y = -0.2, label =c("LRP:",LRPs), size=4)
  # annotate("text", x = c(1.5,2:6), y = -0.2, label =c("LRP:",LRPs), size=4)



ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareEscCor-recCorSca.png",sep=""), plot = g, #-recCorSca#-cvER#-Ages
       width = 4, height = 3, units = "in")

# summary(SpwnCorr.df %>% filter(OM_Name=="cvER0.21.cvERSMU0.42.agePpnConst.recCorSca0.1.n100.mcmc") %>% pull(SpwnCorrValues))
# summary(corMat[lower.tri(corMat)])

# ===================================================================
# (10) Make histograms of cvER
# ==================================================================
nTrials <- 100000
canERlabel <- 0.3
# Density of ERs with cvER=0.21 (either interannual or among CUs)
canER <- canERlabel
cvER <- 0.21
sigCanER <- cvER*canER
shape1 <- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
shape2 <-shape1 * (1/canER-1)
out1 <- rbeta(nTrials,shape1,shape2)

# Density of ERs with cvER=0.42
cvER <- 0.42
sigCanER <- cvER*canER
shape1 <- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
shape2 <-shape1 * (1/canER-1)
out2 <- rbeta(nTrials,shape1,shape2)

#Density of ERs with cvER=0.42 interannually and 0.21 among CUs
canER <- out2
cvER <- 0.21
sigCanER <- cvER*canER
shape1 <- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
shape2 <- shape1 * (1/canER-1)


sampBeta<-function(nTrial) {
  x<-rbeta(1,shape1[nTrial],shape2[nTrial])
}

out3 <- sapply(1:nTrials,sampBeta)
out3b <- sapply(1:nTrials,sampBeta)
out3c <- sapply(1:nTrials,sampBeta)
out3d <- sapply(1:nTrials,sampBeta)
out3e <- sapply(1:nTrials,sampBeta)

#Density of ERs with cvER=0.42 interannually and 0.42 among CUs
canER <- out2
cvER <- 0.42
sigCanER <- cvER*canER
shape1 <- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
shape2 <- shape1 * (1/canER-1)


out4<-sapply(1:nTrials,sampBeta)

out <- data.frame( cvER = c(rep("0.21",nTrials), rep("0.42",nTrials),
                            rep("0.21 x 0.42",nTrials), rep("0.42 x 0.42",nTrials)),
                   ExploitationRate = c(out1, out2, out3, out4) )

g1 <- out %>% dplyr::filter(cvER=="0.21") %>%
  ggplot(aes(ExploitationRate, colour = cvER, fill = cvER)) +
  geom_density (alpha = 0.1) +
  geom_vline (xintercept = canERlabel) +
  xlim (0,1) +
  theme(axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18))

g1b <- out %>% dplyr::filter(cvER=="0.21") %>%
  ggplot(aes(ExploitationRate, colour = cvER, fill = cvER)) +
  geom_density (alpha = 0.1, size=3) +
  geom_vline (xintercept = canERlabel) +
  xlim (0,1) + theme(legend.position="none", panel.grid = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank())



ggsave(paste(wcviCKDir,"/Figures/histcvER0.21.png",sep=""), plot = g1,
         width = 8, height = 6, units = "in")
ggsave(paste(wcviCKDir,"/Figures/histcvER0.21ppt.png",sep=""), plot = g1b,
       width = 3, height = 2, units = "in")

g2 <- out %>% dplyr::filter(cvER=="0.21"|cvER=="0.42") %>%
  ggplot(aes(ExploitationRate, colour = cvER, fill = cvER)) +
  geom_density (alpha = 0.1) +
  xlim (0,1) +
  geom_vline (xintercept = canERlabel) +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))

g2b <- out %>% dplyr::filter(cvER=="0.42") %>%
  ggplot(aes(ExploitationRate)) +
  geom_density (alpha = 0.1, size=3, colour="aquamarine3", fill="aquamarine3") +
  geom_vline (xintercept = canERlabel) +
  xlim (0,1) + ylim(0,6.5) +
  theme(legend.position="none", panel.grid = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank())

ggsave(paste(wcviCKDir,"/Figures/histcvER0.42.png",sep=""), plot = g2,
       width = 8, height = 6, units = "in")
ggsave(paste(wcviCKDir,"/Figures/histcvER0.42pptx.png",sep=""), plot = g2b,
       width = 4, height = 3, units = "in")

g3 <- out %>% dplyr::filter(cvER=="0.21"|cvER=="0.42"|cvER=="0.21 x 0.42") %>%
  ggplot(aes(ExploitationRate, colour = cvER, fill = cvER)) +
  geom_density (alpha = 0.1) +
  xlim (0,1) +
  geom_vline (xintercept = canERlabel) +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))

ggsave(paste(wcviCKDir,"/Figures/histcvER0.21x0.42.png",sep=""), plot = g3,
       width = 8, height = 6, units = "in")

g4 <- out %>%
  ggplot(aes(ExploitationRate, colour = cvER, fill = cvER)) +
  geom_density (alpha = 0.1) +
  xlim (0,1) +
  geom_vline (xintercept = canERlabel) +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18))

ggsave(paste(wcviCKDir,"/Figures/histcvER0.42x0.42.png",sep=""), plot = g4,
       width = 8, height = 6, units = "in")

lineAggData <- out %>% dplyr::filter(cvER=="0.42") %>% slice(1:20)  %>% add_column(Year=1:20) %>% add_column(Label="SMU")


inletLineData <- out %>% dplyr::filter(cvER=="0.21 x 0.42") %>% slice(1:20) %>% add_column(Year=1:20) %>% add_column(Label="Inlet")
inletLineDatab <- data.frame(Label="Inletb", Year=1:20, ExploitationRate=out3b[1:20], cvER="0.21 x 0.42")
inletLineDatac <- data.frame(Label="Inletc", Year=1:20, ExploitationRate=out3c[1:20], cvER="0.21 x 0.42")
inletLineDatad <- data.frame(Label="Inletd", Year=1:20, ExploitationRate=out3d[1:20], cvER="0.21 x 0.42")
inletLineDatae <- data.frame(Label="Inlete", Year=1:20, ExploitationRate=out3e[1:20], cvER="0.21 x 0.42")
inletLineData <- inletLineData %>% add_row(inletLineDatab) %>% add_row(inletLineDatac) %>% add_row(inletLineDatad) %>% add_row(inletLineDatae)
lineAggData <- lineAggData %>% add_row(inletLineData) %>% dplyr::select(-cvER)


g5 <- lineAggData %>% filter(Label == "SMU") %>%
  ggplot(aes(Year, ExploitationRate, colour=Label)) +
  scale_colour_manual(values=viridis(3)[1]) +
  ylim(0,0.75)+
  geom_line(size=2) +
  theme(axis.title.x=element_text(18),
        axis.title.y=element_text(18))
g5

ggsave(paste(wcviCKDir,"/Figures/ERtimeseries1.png",sep=""), plot = g5,
       width = 8, height = 6, units = "in")

g6 <- lineAggData %>% filter(Label=="SMU"|Label=="Inlet") %>% ggplot(aes(Year, ExploitationRate, colour=Label)) +
  scale_colour_manual(values=c(viridis(3)[2],viridis(3)[1])) +
  geom_line(size=2) +
  ylim(0,0.75)+
  theme(axis.title.x=element_text(18),
        axis.title.y=element_text(18))

g6

ggsave(paste(wcviCKDir,"/Figures/ERtimeseries2.png",sep=""), plot = g6,
       width = 8, height = 6, units = "in")

g7 <- lineAggData  %>% ggplot(aes(Year, ExploitationRate, colour=Label)) +
  scale_colour_manual(values=c(rep(viridis(3)[2],5),viridis(3)[1])) +
  geom_line(size=2) +
  ylim(0,0.75)+
  theme(axis.title.x=element_text(18),
        axis.title.y=element_text(18))

g7

ggsave(paste(wcviCKDir,"/Figures/ERtimeseries3.png",sep=""), plot = g7,
       width = 8, height = 6, units = "in")


# ===================================================================
# (11) Plots of LRP stabilization with number of trials
# ==================================================================

# Specify threshold to use when calculating LRP
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark

probThresh<-c(0.5, 0.66, 0.9, 0.99)# probability threshold; the LRP is set as the
# aggregate abundance that has this probability that the propCUThreshold is met

OMsToInclude<-c("cvER0.21n5000")



# Read in samSim outputs for OM
filename<-paste("projLRPDat_",OMsToInclude,".csv",sep="")
projLRPDat<-read.csv(here(wcviCKDir, "SamSimOutputs", "simData",filename))
CUpars <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))

nTrials <- 50000#2000

# Loop over probability levels
for (p in 1:length(probThresh)){
  # Loop over nTrials
  for (i in 1:500) {#(i in 1:200) {

    projLRPDat<-read.csv(here(wcviCKDir, "SamSimOutputs", "simData",filename))
    projLRPDat<-projLRPDat %>% filter(year > CUpars$ageMaxRec[1]*100) %>% filter(iteration < (i*100))

    # Create bins for projected spawner abundances
    minBreak<-0
    maxBreak<-round(max(projLRPDat$sAg),digits=-2)
    binSize<-200 # Note: bin size is currently set here
    breaks<-seq(minBreak, maxBreak,by=binSize)

    # Set bin labels as the mid-point
    projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))

    # Summarize nSims in each bin
    tmp<-projLRPDat %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))

    # Filter out bins with < 100 nSims
    minSims <- nTrials/100
    tmp2<-projLRPDat %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
      add_column(nSims=tmp$nSims) %>% filter(nSims>=minSims)

    # For each bin, calculate probability that required proportion of CUs above benchmark
    projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
    # For each bin, calculate the difference between the threshold probability and the calculated probability
    projLRPDat$diff<-abs( probThresh[p] - projLRPDat$prob )

    # Save projection summaries used to create plots
    projLRPDat$OM.Name<-OMsToInclude
    if (i == 1) projLRPDat.plot<-projLRPDat
    if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)

    # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
    LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))

    # Create a table of LRP estimates to be saved for each OM model
    if (i ==1) {
      LRP_Ests_nTrials <- NA
      LRP_Ests_nTrials<-data.frame(OMsToInclude, i*10, probThresh[p], propCUThresh, LRP, binSize)
      names(LRP_Ests_nTrials)<-c("OM", "nTrials","ProbThresh", "PropCURequired", "LRP", "binSize")
    } else {
      tmp.df<-data.frame(OMsToInclude, i*10, probThresh[p], propCUThresh, LRP, binSize)
      names(tmp.df)<-c("OM", "nTrials", "ProbThresh", "PropCURequired", "LRP", "binSize")
      LRP_Ests_nTrials<-rbind(LRP_Ests_nTrials,tmp.df)
    }
    print(c(probThresh[p],i,LRP))
  }

  # Save LRPs
  write.table(LRP_Ests_nTrials, paste(projOutDir2, "/ProjectedLRPscvER0.21_nTrials_p", probThresh[p], ".csv", sep=""),
                row.names=F)

  # Plot stabilization of LRPs with nTrials
  #LRP_nTrials <- as.data.frame(LRP_Ests_nTrials) %>% filter(ProbThresh==0.5)
  LRP_nTrials <- as.data.frame(read.table(paste(projOutDir2, "/ProjectedLRPscvER0.21_nTrials_p",probThresh[p],".csv", sep=""), header=T))
  LRP_nTrials_plot <- ggplot(LRP_nTrials, aes(nTrials,LRP)) + geom_line() #+
  #ylim(14000,19000) #+
  #geom_vline(xintercept=100, linetype="dashed")


  ggsave(paste(wcviCKDir,"/Figures/LRPcvER0.21_nTrials_p",probThresh[p],".png",sep=""), plot = LRP_nTrials_plot,
         width = 6, height = 4, units = "in")


}




# ===================================================================
# (12) Plot LRPs with various plevels
# ==================================================================

# Specify threshold to use when calculating LRP
# # Note: may want to loop over probThresholds as well; still needs to be added
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-c(0.50,0.66)#,0.9, 0.99) # probability theshhold; the LRP is set as the aggregate abundance that has this
# probability that the propCUThreshold is met

# Specify scenarios to calculate LRPs and make plots for.
# These scenarios will be looped over below with a LRP (and LRP plot) saved for each scenario
OMsToInclude<-c(
  # "baseER")
  #"ER0",
  "ER0.05",
  "ER0.10",
  "ER0.15",
  "ER0.2",
  "ER0.25",
  "baseER",
  "ER0.35",
  "ER0.4",
  "ER0.45")
  # "alphaScalar0.75",
  # "baseERn10000",
  # "alphaScalar1.5")
  # "cvER0",
  # "baseER",
  # "cvER0.17")




if(length(OMsToInclude)==1) OMsToIncludeName <- OMsToInclude[1]
if(length(OMsToInclude)==9) OMsToIncludeName <- "ERs"
if(length(OMsToInclude)==3) OMsToIncludeName <- "cvER"#"cvER"#"Alphas"#"cvER"#"

LRP <- NA

for (OM in 1:length(OMsToInclude)){

  # Loop over OM Scenarios
  for (i in 1:length(probThresh)) {

    # Read in samSim outputs for OM
    filename<-paste("projLRPDat_",OMsToInclude[OM],".csv",sep="")
    projLRPDat<-read.csv(here(wcviCKDir, "SamSimOutputs", "simData",filename))
    CUpars <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))
    projLRPDat<-projLRPDat %>% filter(year > CUpars$ageMaxRec[1]*10)#)max(SRDat$yr_num)+4)

    # Create bins for projected spawner abundances
    minBreak<-0
    maxBreak<-round(max(projLRPDat$sAg),digits=-2)
    binSize<-200 # Note: bin size is currently set here
    breaks<-seq(minBreak, maxBreak,by=binSize)

    # Set bin labels as the mid-point
    projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))

    # Summarize nSims in each bin
    tmp<-projLRPDat %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))

    # Filter out bins with < 100 nSims
    tmp2<-projLRPDat %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
      add_column(nSims=tmp$nSims) %>% filter(nSims>=10)

    # For each bin, calculate probability that required proportion of CUs above benchmark
    projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
    # For each bin, calculate the difference between the threshold probability and the calculated probability
    tmp3 <- projLRPDat %>% filter(nSims>100)# Remove bins where there are very few nSims among LRP options
    min <- min(abs(probThresh[i]-tmp3$prob))

    projLRPDat$diff<-abs(probThresh[i]-projLRPDat$prob)

    # Save projection summaries used to create plots
    projLRPDat$OM.Name<-OMsToInclude[OM]
    if (i == 1) projLRPDat.plot<-projLRPDat
    if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)

    # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
    #LRP[i]<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
    LRP[i]<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min]))

    # Create a table of LRP estimates to be saved for each OM model
    if (i ==1) {
      LRP_Ests<-data.frame(OMsToInclude[OM], probThresh[i], propCUThresh, LRP[i], binSize)
      names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
    } else {
      tmp.df<-data.frame(OMsToInclude[OM], probThresh[i], propCUThresh, LRP[i], binSize)
      names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
      LRP_Ests<-rbind(LRP_Ests,tmp.df)
    }

    if(i==1){# Plot projected LRP abundance relationship ===============================================================
      if(OM==1) {
        if(length(OMsToInclude)==1|length(OMsToInclude)==9) {
          plot.width <- 5
          plot.height <- 4
        }
        if(length(OMsToInclude)==3) {
          plot.width <- 5
          plot.height <- 1.5
        }

        png(paste(wcviCKDir,"/Figures/ProjectedLRPs/", OMsToIncludeName,
                  "-ProjLRPCurve-ALLp.png", sep=""), width=plot.width,
            height=plot.height,
            units="in", res=500)
        if(length(OMsToInclude)==9) layout(matrix(c(1:9), 3, 3, byrow = TRUE))
        if(length(OMsToInclude)==3) layout(matrix(c(1:3), 1, 3, byrow = TRUE))

      }# End of if(OM==1) {


    if(length(OMsToInclude)==1){
      plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19,
           xlim=c(0, max( as.numeric(as.character(projLRPDat$bins)),
                          na.rm=T)*1.0 ),
           ylim=c(0,1),
           cex=0.5, cex.lab=1,#1.5,
           xlab="Aggregate Abundance", ylab="Pr (All inlets > Lower Benchmark)")
    }# End of if(length(OMsToInclude)==1){

    if(length(OMsToInclude)==9){
      par(mar=c(2.8,3,0.6,1))
      xMax <- 50000
      if(OM<7){
        xaxt <- "n"#par(xaxt="n")
      }
      if(OM>=7){
        xaxt <- "s"#par(xaxt="s")
      }
    }# End of if(length(OMsToInclude)==9){
    if(length(OMsToInclude)==3){
        par(mar=c(2.8,3,1,1))
        xaxt <- "s"
        xMax <- 70000
      }

      if(length(OMsToInclude)>1){
        plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19,
             xlim=c(0, xMax ),
             ylim=c(0,1),
             cex=0.3, cex.lab=1,#1.5,
             xlab="", ylab="", xaxt=xaxt)

      }


      if(length(OMsToInclude)==9){
        if(OM<7){
          at1 <- seq(0, 50000, 10000)
          axis(side =1,  at=at1, labels = FALSE)
        }


        panel.title <- c("5%", "10%", "15%", "20%", "25%", "30%", "35%", "40%",
                         "45%")
        mtext(text=panel.title[OM], side=3, line=0, at=5000, cex=0.4)

        LRP_50 <- LRP_Ests$LRP[1]#(read.csv(paste(wcviCKDir,
                  #              "/DataOut/ProjectedLRPs/ProjectedLRPs",
                  #              OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
                  # pull(LRP))[1]

        LRP_66 <- (read.csv(paste(wcviCKDir,
                                  "/DataOut/ProjectedLRPs/ProjectedLRPs",
                                  OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
                     pull(LRP))[2]
        text(x=35000, y=0.15, labels=paste("LRP(p=0.5)= ", LRP_50), cex=0.4)
        if(OM<7) text(x=35000, y=0.05, labels=paste("LRP(p=0.66)= ", LRP_66), cex=0.4)

        if(OM==4) {mtext("Probability of all inlets > lower benchmark", side=2,
                        line=1.8,at=0.5, cex=1) }
        if(OM==8) {mtext("Aggregate Abundance", side=1, line=1.8, at=40000,
                         cex=0.7) }

      }# End of if(length(OMsToInclude)==9){


    if(length(OMsToInclude)==3){

      # panel.title <- c("0.75 x productivity", "Base productivty",
      #                  "1.5 x productivity")
      panel.title <- c("cvER = 0", "cvER = 0.085",
                       "cvER = 0.17")
      mtext(text=panel.title[OM], side=3, line=0, at=20000, cex=0.5)

      LRP_50 <- (read.csv(paste(wcviCKDir,
                              "/DataOut/ProjectedLRPs/ProjectedLRPs",
                              OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
                  pull(LRP))[1]
      LRP_66 <- (read.csv(paste(wcviCKDir,
                                "/DataOut/ProjectedLRPs/ProjectedLRPs",
                                OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
                  pull(LRP))[2]

      text(x=40000, y=0.15, labels=paste("LRP(p=0.5)= ", LRP_50), cex=0.4)#if (OM<3): alpha
      text(x=40000, y=0.05, labels=paste("LRP(p=0.66)= ", LRP_66), cex=0.4)# if (OM>1): alpha

      if(OM==1) {mtext("Prob(all inlets)>lower benchmark", side=2,
                       line=1.8,at=0.4, cex=0.55) }
      if(OM==2) {mtext("Aggregate Abundance", side=1, line=1.8, at=40000,
                       cex=0.7) }

    }# End of if(length(OMsToInclude)==3){

    }# End of if(i==1){

    if (length(OMsToInclude == 9)) lrp.lwd <- 1
    if (length(OMsToInclude != 9)) lrp.lwd <- 2
    abline(h=probThresh[i], lty=2, lwd=lrp.lwd)
    if(OMsToInclude[OM]!="alphaScalar1.5") { if (i==1)
      abline(v=LRP[i], col="orange", lwd=lrp.lwd) }
    if(OMsToInclude[OM]!="alphaScalar0.75"&OM<7) { if (i==2)
      abline(v=LRP[i], col=viridis(4, alpha=0.3)[3], lwd=lrp.lwd) }
    abline(h=0.9, lty=2)
    abline(h=0.99, lty=2)

    # if(i==3) abline(v=LRP[i], lwd=4, col=viridis(4, alpha=0.3)[2] )
    # if(i==4) abline(v=LRP[i], lwd=4, col=viridis(4, alpha=0.2)[1] )

    if(i==length(probThresh)) {
      if(OM==length(OMsToInclude)) {
        dev.off()
      }
    }

    # Option to plot histogram of nSims in each Agg Abundance Bin
    #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)

  }

  # Save LRPs for all OM scenarios
  write.csv(LRP_Ests, paste(projOutDir2, "/ProjectedLRPs",  OMsToInclude[OM],
                            "_ALLp.csv", sep=""), row.names=F)
  # Save LRP projection summaries used for calculating and plotting LRP (Optional)
  write.csv(projLRPDat.plot, paste(projOutDir2, "/ProjectedLRP_data", OMsToInclude[OM],
                                   "_Allp.csv", sep=""), row.names=F)


}# End of for OM in 1:length(OMsToInclude)


# ===================================================================
# (13) Run reconstruction for WCIV  CK
# ==================================================================
# This code runs run reconstrution from Diana Dobson, updated by Diana Mchugh
# 1 June 2021"WCVI_term_model_revisions_updated-2021.xlsx" in TMBstan

run.RunReconstruction <- FALSE
setwd(codeDir)

if(run.RunReconstruction){
  SRDat <- data.frame(read.csv(paste(wcviCKDir, "/DataIn/WCVI_SRbyCU.csv", sep="")))
  #SRDat <- data.frame(read.csv("DataIn/WCVI_SRbyCU.csv"))
  SRDat <- SRDat %>% drop_na()
  # stkName <- unique(SRDat$Stock)
  # # # for (i in 1:length(stkName)){}
  # # i <- 1
  # # srDat %>% filter(Stock == stkName[i])
  #
  # for (i in (1:3)){
  #   r <- SRDat %>% filter(CU_ID==0)%>% pull(Recruits)
  #   s <- SRDat %>% filter(CU_ID==0)%>% pull(Spawners)
  #   lm(log(r/s)~s)
  # }
  Mod <- "SR_RickerBasic"

  Scale <- 1000

  data <- list()
  data$biasCor <- as.numeric(TRUE)
  data$S <- SRDat$Spawners/Scale
  data$logR <- log(SRDat$Recruits/Scale)
  data$stk <- as.numeric(SRDat$CU_ID)
  N_Stocks <- length(unique(SRDat$CU_Name))
  #data$N_Stks <- N_Stocks
  #data$yr <- SRDat$BroodYear

  param <- list()
  param$logA <- rep(1, N_Stocks)
  param$logB <- as.vector(log(1/( (SRDat %>% group_by(CU_ID) %>%
                                     summarise(x = quantile(
                                       Spawners, 0.8, na.rm=T)))$x/Scale) ) )
  param$logSigma <- rep(-2, N_Stocks)


  dyn.unload(dynlib(paste("TMB_Files/",Mod, sep="")))
  compile(paste(codeDir, "/TMB_Files/", Mod, ".cpp", sep=""))

  dyn.load(dynlib("TMB_Files/SR_RickerBasic"))

  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
  pl <- obj$env$parList(opt$par)
  pl
  #summary(sdreport(obj), p.value=TRUE)

  ricA.rr <- data.frame(CU_Name=unique(SRDat$CU_Name), logRicA=pl$logA)
  write.csv(ricA.rr, paste(wcviCKDir,"/DataIn/ricArr.csv", sep=""))

}

# ===================================================================
# (14) Plot time-series of age proportions in recruitment
# ==================================================================

cuAges <- read.csv(paste(wcviCKDir, "/DataIn/CUages.csv", sep=""))
cuAges <- cuAges %>% pivot_longer(cols=c("age2", "age3", "age4", "age5"),
                                  names_to="Age", values_to="Proportion")

cuAges$CU_Names. <- factor(cuAges$CU_Names,
                            levels=c("Southwest_Vancouver_Island",
                                     "Nootka_Kyuquot",
                                     "Northwest_Vancouver_Island",
                                     "Westcoast_Vancouver_Island"))

cuAges %>% filter(CU_Names!="Westcoast_Vancouver_Island") %>%
  ggplot(aes(Year, Proportion, group=CU_Names., colour=CU_Names.)) +
  geom_line() + facet_wrap(~Age, ncol=1)



# ===================================================================
# (15) Look at distribution of alphas used in samSim
# ==================================================================
scenarioToLook <- "cvER0.21n5000"
d <- readr::read_rds(paste(wcviCKDir, "/SamSimOutputs/simData/", scenarioToLook,
                           "/", scenarioToLook,"_baseER_cuDat.RData", sep=""))
names(d)
d["medAlpha"][[1]] #  matrix(..., nrow = nTrials, ncol = nCU), median over years
d["varAlpha"][[1]] # matrix(..., nrow = nTrials, ncol = nCU), cv over years

# medAlpha[n, ] <- apply(na.omit(alphaMat[yrsSeq, ]), 2, median)
# varAlpha[n, ] <- apply(na.omit(alphaMat[yrsSeq, ]), 2, cv)


# ===================================================================
# (16) Make Comparison Plots Among Scenarios (NOT CURRENTLY WORKING)
# ==================================================================

# Note: The below code needs to be updated for new projected LRP method (Apr 26, 2021)


# # Plot to compare LRP among different SRR structures =====================
#
# OMsToPlot<-c("IM.Base", "HM.Base", "IMCap.base", "HMCap.base")
#
# p<-1.0
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#              geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#              xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#               scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                                labels=c("IM", "IM.Cap", "HM", "HM.Cap"))
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_byOM_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
# p<-0.8
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap"))
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_byOM_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
# # Plot to show combined OM scenarios =====================================
#
# OMsToPlot<-c("IM.Base","IMCap.base", "CombinedIM")
#
# p<-1.0
#
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
#
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base", "CombinedIM"),
#                    labels=c("IM", "IM.Cap", "IM.Composite"))
#
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_compositeIM_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
# # Plot to show sensitivity analysis to recruitment correlation
#
# #OMsToPlot<-c("IM.Base","IM.80RecCor","IM.60RecCor","IM.40RecCor","IM.10RecCor")
#
# OMsToPlot<-c("IM.Base","IM.60RecCor","IM.40RecCor","IM.10RecCor", "IM-.20RecCor","IM-.40RecCor")
#
# p<-1.0
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base","IM.60RecCor","IM.40RecCor","IM.10RecCor", "IM-.20RecCor","IM-.40RecCor"),
#                    labels=c("Base(MPD)", "60%Corr", "40%Corr","20%Corr","-20%Corr","-40%Corr"))
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_SAnalRecCorr_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
# p<-0.8
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base","IM.80RecCor","IM.60RecCor","IM.40RecCor","IM.10RecCor"),
#                    labels=c("Base(MPD)", "80%Corr", "60%Corr","40%Corr","10%Corr"))
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_SAnalRecCorr_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
#
# # ---- Temporary code to test effect of correlation scalar on projected correlations
#
# MPD_recCormat<-read.csv(paste(wcviCKDir,"/SamSimInputs/IM.base/cohoCorrMat.csv", sep=""), header=F)
#
# OMsToTest<-c("IM.Base","IM.80RecCor","IM.60RecCor","IM.40RecCor","IM.40RecCor","IM-.40RecCor")
#
# #OMsToTest<-c("IM.cvER1.25", "IM.cvER1.5","IM.cvER2.0","IM.cvER2.5","IM.cvER3.0")
#
# filename<-paste( "projSpwnDat_",OMsToTest[6],".csv",sep="")
# spDat<-read.csv(here(wcviCKDir,"SamSimOutputs", "simData",filename))
#
# # Calculate correlation matrix in MPD recruitment residuals ========================
# #resids<-as_tibble(data.frame(stock=data$stk,year=data$yr, resid=obj$report()$R_Resid))
# spDat<-as_tibble(spDat)
# spDat<-spDat%>%dplyr::select(-X)
#
# RecCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))
# SpwnCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))
#
# for (i in 1:max(spDat$iteration)) {
#
#   recruits.i<-spDat %>% filter(iteration==i & expRate==0.125) %>%  dplyr::select(-spawners)
#   recruits.i<-recruits.i %>% dplyr::select(-expRate, -iteration)
#   cor_mat<-recruits.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=recruits) %>% dplyr::select(-year) %>% cor()
#   RecCorMat[,,i]<-cor_mat
#
#   spawners.i<-spDat %>% filter(iteration==i & expRate==0.125) %>%  dplyr::select(-recruits)
#   spawners.i<-spawners.i %>% dplyr::select(-expRate, -iteration)
#   cor_mat<-spawners.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=spawners) %>% dplyr::select(-year) %>% cor()
#   SpwnCorMat[,,i]<-cor_mat
#
# }
#
#
# RecCorMat_Ave<-apply(RecCorMat, c(1,2), mean)
# SpwnCorMat_Ave<-apply(SpwnCorMat, c(1,2), mean)
#
#
#
#
# # Plot to show sensitivity analysis to ER variability
#
# OMsToPlot<-c("IM.Base","IM.cvER1.5","IM.cvER2.0","IM.cvER2.5","IM.cvER3.0")
#
# p<-1.0
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base", "IM.cvER1.5","IM.cvER2.0", "IM.cvER2.5","IM.cvER3.0"),
#                    labels=c("Base", "1.5CV","2.0CV", "2.5CV","3.0CV"))
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_SAnalCvER_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
# p<-0.8
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base", "IM.cvER1.5","IM.cvER2.0", "IM.cvER2.5","IM.cvER3.0"),
#                    labels=c("Base", "1.5CV","2.0CV", "2.5CV","3.0CV"))
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareLRP_SAnalCvER_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
#
# # Plot to compare projected and empirical LRPs
#
# OMsToPlot<-c("IM.Base", "HM.Base", "IMCap.base", "HMCap.base")
#
# p<-1.0
# projPlotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# projPlotDat <- projPlotDat %>% add_column(Method=rep("Projected", 4))
#
# # Empirical LRPs with p = 0.99 and likelihood penalty
# HM.Cap<-c(32962.96235,	28565.95578,	37359.96893)
# IM<-c(24358.94812,	20298.83689,	28419.05935)
# IM.Cap<-c(33787.86878,	29323.28802,	38252.44954)
# HM<-c(19055.00439,	15281.82266,	22828.18612)
#
# empiricalPlotDat<-data.frame(OM.Name=c("HM.Base", "HMCap.base","IM.Base", "IMCap.base"),
#                              ppnCUsLowerBM = rep(1,4),
#                              LRP.50=rep(NA,4),LRP.95=rep(NA,4),LRP.05=rep(NA,4),
#                              Method = "Empirical")
#
# empiricalPlotDat[1,3:5] <- HM
# empiricalPlotDat[2,3:5] <- HM.Cap
# empiricalPlotDat[3,3:5] <- IM
# empiricalPlotDat[4,3:5] <- IM.Cap
#
# empiricalPlotDat<-as_tibble(empiricalPlotDat)
#
# plotDat<-bind_rows(projPlotDat, empiricalPlotDat)
#
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50, colour=as.factor(Method))) +
#   geom_point(position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0,position=position_dodge(width=0.5)) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap")) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   labs(colour = "Method")
#
#
#
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareProjVsEmp_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
# # ---- with p = 0.80
# p<-0.80
# projPlotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# projPlotDat <- projPlotDat %>% add_column(Method=rep("Proj", 4))
#
#
# # Empirical LRPs with p = 0.80 and likelihood penalty
# HM.Cap<-c(22246.52014, 19686.91885, 24806.12142)
# IM<-c(16702.85366, 14317.14271, 19088.56461)
# IM.Cap<-c(22765.51884, 20171.41927,25359.61841)
# HM<-c(13624.39163,	11523.63963,	15725.14363)
#
# empiricalPlotDat<-data.frame(OM.Name=c("HM.Base", "HMCap.base","IM.Base", "IMCap.base"),
#                              ppnCUsLowerBM = rep(1,4),
#                              LRP.50=rep(NA,4),LRP.95=rep(NA,4),LRP.05=rep(NA,4),
#                              Method = "Empirical")
#
# empiricalPlotDat[1,3:5] <- HM
# empiricalPlotDat[2,3:5] <- HM.Cap
# empiricalPlotDat[3,3:5] <- IM
# empiricalPlotDat[4,3:5] <- IM.Cap
#
# empiricalPlotDat<-as_tibble(empiricalPlotDat)
#
# plotDat<-bind_rows(projPlotDat, empiricalPlotDat)
#
# g<-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50, colour=as.factor(Method))) +
#   geom_point(position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0,position=position_dodge(width=0.5)) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap")) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   labs(colour = "Method")
#
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareProjVsEmp_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
# # Empirical regression using p = 0.999
# p<-1.0
# projPlotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# projPlotDat <- projPlotDat %>% add_column(Method=rep("Projected", 4))
#
# # Empirical LRPs with p = 0.99 and likelihood penalty
# HM.Cap<-c(32962.96235,	28565.95578,	37359.96893)
# IM<-c(24358.94812,	20298.83689,	28419.05935)
# IM.Cap<-c(33787.86878,	29323.28802,	38252.44954)
# HM<-c(19055.00439,	15281.82266,	22828.18612)
#
# empiricalPlotDat<-data.frame(OM.Name=c("HM.Base", "HMCap.base","IM.Base", "IMCap.base"),
#                              ppnCUsLowerBM = rep(1,4),
#                              LRP.50=rep(NA,4),LRP.95=rep(NA,4),LRP.05=rep(NA,4),
#                              Method = "Empirical")
#
# empiricalPlotDat[1,3:5] <- HM
# empiricalPlotDat[2,3:5] <- HM.Cap
# empiricalPlotDat[3,3:5] <- IM
# empiricalPlotDat[4,3:5] <- IM.Cap
#
# empiricalPlotDat<-as_tibble(empiricalPlotDat)
#
# plotDat<-bind_rows(projPlotDat, empiricalPlotDat)
#
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50, colour=as.factor(Method))) +
#   geom_point(position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0,position=position_dodge(width=0.5)) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap")) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   labs(colour = "Method")
#
#
#
# ggsave(paste(wcviCKDir,"/Figures/ProjectedLRPs/compareProjVsEmp_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")
#
#
#
# # Plot correlation matrix to show scenarios
# library(ggcorrplot)
# cormat<-read.csv(paste(wcviCKDir,"/SamSimInputs/IM.base/cohoCorrMat.csv", sep=""), header=F)
# rownames(cormat) <- unique(SRDat$CU_Name)
# colnames(cormat) <- unique(SRDat$CU_Name)
# ggcorrplot(cormat,
#            hc.order = TRUE,
#            type = "lower",
#            outline.color = "white",
#            lab=T)
#
#
#
#
# # Plot ER distribution
#
# canER<-0.125
# cvER<-0.456
# sigCanER<-cvER*canER
#
# shape1<- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
# shape2<-shape1 * (1/canER-1)
#
# sampBeta<-rbeta(1000,shape1,shape2)
