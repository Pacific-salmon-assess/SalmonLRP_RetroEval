# Code by K.Holt & B.Davis
# Started on: February 7, 2020

# This file contains the code required to explore Interior Fraser Coho escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions

# ==============================================================================================
# STEPS FOR RUNNING ANALYSES IN THIS FILE 
# ===================================================
#
# (1) Read-in and format IF Coho data
# (2) Call functions to plot data
# (3) Set-up TMB inputs for logistic regression fits
# (4) RETROSPECTIVE ANALYSIS: Run annual retrospective analyses
#   (4.1) Print saved logistic model fit diagnostics for any of the above runs
#   (4.2) Create plot to compare retrospective time series of LRP estimates
# (5) EFFECT OF MISSING CUS: Run all available combinations of number of CUs
#   (5.1) Make Plots for Aggregate Status as a function of nCUS 
# (6) Create figures with two or more probabilities on a logistic fit for a specified year


# ===================================================================================================

#s<-FALSE
useFrenchCaptions<-TRUE

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(tmbstan)
library(here)


setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
  source("LRdiagnostics.r")
}

sourceAll()

# Fit SR models for each CU and save results in DataOut folder
source(paste(cohoDir,"compareRickerModelTypes_SRonly.r",sep="/"))

# Load TMB models

compile("TMB_Files/SR_HierRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_LowAggPrior"))

compile("TMB_Files/SR_HierRicker_SurvCap_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_SurvCap_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_LowAggPrior"))

compile("TMB_Files/SR_HierRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv"))

compile("TMB_Files/SR_IndivRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv"))

compile("TMB_Files/SR_HierRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap"))

compile("TMB_Files/SR_IndivRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap"))

# Only needed for leave-one-out logistic regression diagnostics
compile("TMB_Files/LRP_BasicLogistic_Only_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/LRP_BasicLogistic_Only_LowAggPrior"))

# Only needed for leave-one-out logistic regression diagnostics
compile("TMB_Files/LRP_Logistic_Only_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/LRP_Logistic_Only_LowAggPrior"))

compile("TMB_Files/LRP_Logistic_Only.cpp")
dyn.load(dynlib("TMB_Files/LRP_Logistic_Only"))

# ======================================================================
# (1) Read-in and format IF Coho data:  
# =====================================================================
setwd(cohoDir)

# CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
#   # Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
#   colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
#   colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
#   colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
#   colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")

CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")

# Change column names to yr, CU_Name, Escp, Subpop_Name, CU_ID
  CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(MU_Name=MU_Name, yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
  tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
  CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)


# Summary of CU-level escapements based on Natural.Spawners  
CoEscpDat <- as_tibble(CoEscpDat_bySubpop) %>% group_by(MU_Name, CU_Name, CU_ID, yr) %>% select(MU_Name, CU_Name, CU_ID, yr, Escp) %>%summarize(Escp=sum(Escp))
CoEscpDat <- as.data.frame(CoEscpDat)
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"

  
setwd(codeDir)

# Run annual restrospective analyses over various levels of p, restrict escapement dataset to 1998+ ============================================================
# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)
CoEscpDat_bySubpop <- CoEscpDat_bySubpop %>% filter(yr >= 1998)

# Roll up escpaments, and get Gen Mean of that
AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))



# ======================================================================================================================
# (2) Call functions to plot data:
# ===================================================================================================================

# Ensure 'Figures' output directory exists
figureDir<-paste(cohoDir,"/Figures", sep="")
if (!file.exists(figureDir)){
  dir.create(figureDir)
}

 plot_CU_Escp_Over_Time(CoEscpDat, cohoDir, plotName="coho-CU-EscpSeries", samePlot = F, withSgen=TRUE, addGenMean=T,
                        SgenFileName="ModelFits/AllEsts_Indiv_Ricker_Surv")

 plot_CU_Escp_Over_Time(Dat=CoEscpDat, Dir=cohoDir, plotName="coho-CU-EscpSeries-Combined", samePlot = TRUE, withSgen=FALSE, addGenMean=FALSE,
                       SgenFileName="ModelFits/AllEsts_Indiv_Ricker_Surv")
 
 plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="coho-Subpop-EscpSeries", samePlot = F, addGenMean=T)


 plot_CU_Escp_withStatus(CoEscpDat, cohoDir, plotName="coho-CU-EscpSeries-wStatus", SgenFileName="ModelFits/AllEsts_Indiv_Ricker_Surv")
 
 
 plot_Subpop_Escp_withStatus(Dat=CoEscpDat_bySubpop, Dir=cohoDir, plotName="coho-Subpop-EscpSeries-wStatus",useFrenchCaptions=useFrenchCaptions)
 
 
# ==================================================================================================================
# (3) Set-up TMB inputs for logistic regression fits
# =====================================================================================================================

TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                      extra_eval_iter=FALSE,biasCorrect=TRUE)

# Prior means come from running "compareRickerModelTypes_onlySR.r", with bias correction
 ## Note: using expansion by 1.4 for prior cap (not 1.5, like Arbeider et al., Korman et al)
# Current estimates (based on expansion of 1.4):
 cap_priorMean_IM <- c(11.661606, 4.220915, 13.584558, 20.156285, 17.128963)
 # With expansion of 1.5:
# cap_priorMean_IM<-c(12.494575, 4.522409, 14.554897, 21.596031, 18.352463)
# With 1.35
# cap_priorMean_IM<-c(11.245117, 4.070168, 13.099408, 19.436428, 16.517216)



TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_IM, cap_sig=sqrt(2),
                               extra_eval_iter=FALSE,biasCorrect=TRUE)

# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                      logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1, 
                      extra_eval_iter=FALSE,biasCorrect=TRUE)


# Prior means come from running "compareRickerModelTypes_onlySR.r"
## Note: using expansion by 1.40 for prior cap (not 1.5, like Arbeider et al., Korman et al)
cap_priorMean_HM<-  c(10.834784, 4.854988, 12.752003, 19.777548, 15.856521)

## Note: using expansion by 1.35 for prior cap (not 1.5, like Arbeider et al., Korman et al)
#cap_priorMean_HM<-c(10.447827, 4.681596, 12.296574, 19.071207, 15.290217)

TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1,
                               logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1,
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_HM, cap_sig=sqrt(2),
                               extra_eval_iter=FALSE, biasCorrect=T)


TMB_Inputs_Subpop <- list(Scale = 1000, extra_eval_iter=FALSE)



# # Test projected LRP
# devtools::install_github("Pacific-salmon-assess/samSim", ref="LRP")
# ps <- c(0.6, 0.8, 1.0)
# 
# for(pp in 1:length(ps)){
# 
#   runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
#                BMmodel = "SR_IndivRicker_Surv", LRPmodel="Proj", integratedModel=F,
#                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Proj.IndivRickerSurv_",ps[pp]*100, sep=""),
#                bootstrapMode = F, plotLRP=T)
# }


# ===============================================================================
# (4) RETROSPECTIVE ANALYSIS: Run annual retrospective analyses
# ====================================================================================

# Note: if .cpp files are already compiled (.dll and .o files are present in SalmonLRP_RetroEval/Code/TMB_Files ), 
# code below will run and produce figures, but will give this error message: 
# "Error in .Call("FreeADFunObject", ptr, PACKAGE = DLL) : "FreeADFunObject" not available for .Call() for package "<TMB model name>""
# Deleting the .dll and .o files in /TMB_Files folder and recompiling .cpp files gets rid of this error. 

# Note: BroodYrLag = the number of years to subtract from the current year to get to the most recent brood year

# Loop over p values and run annual retrospective analyses for each level of p
  ps <- c(0.5, 0.66, 0.9, 0.99)
  for(pp in 1:length(ps)){
  
     ## Run with Bernoulli LRP model with individual model Ricker
     runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2010, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
                    BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=F, codeDir=codeDir)
           
     # ## Run with Bernoulli LRP model with hierarchical Ricker
     # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
     #                BMmodel = "SR_HierRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
     #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bern.HierRickerSurv_",ps[pp]*100, sep=""),
     #                bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
     
     ## Run with Bernoulli LRP model with individual model Ricker, with prior on capacity
     runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2010, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
                    BMmodel = "SR_IndivRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurvCap_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=F,codeDir=codeDir)

     # Run with Bernoulli LRP model with hierarchical Ricker, with prior on capacity
     # CW : I am getting an error when running this, since we will not be using these for now, I'll comment it out
     #runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
     #            BMmodel = "SR_HierRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
     #            useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bern.HierRickerSurvCap_",ps[pp]*100, sep=""),
     #            bootstrapMode = F, plotLRP=T,runLogisticDiag=T)

     ## Run using distribibutional benchmark for CUs based on 50% of subpopulations within a CU > 1000 fish   
     runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2010, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
                    BMmodel = "ThreshAbund_Subpop1000_ST", LRPmodel="BernLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.SPopAbundThreshST_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=F, codeDir=codeDir)

     # ## CW: Run using multidimensional rapid assessments results  (in progress)  
     # runAnnualRetro(EscpDat=CoEscpDat_bySubpop, #do not know what to insert here
     #                SRDat=NULL, startYr=2015, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
     #                BMmodel = "RapidAssessment_Ricker",#working on this option
     #                 LRPmodel="BernLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
     #                useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.Multidim_Ricker",ps[pp]*100, sep=""),
     #                bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
     # 
     # runAnnualRetro(EscpDat=CoEscpDat_bySubpop, #do not know what to insert here
     #                SRDat=NULL, startYr=2015, endYr=2020, BroodYrLag=4, genYrs=3, p = ps[pp],
     #                BMmodel = "RapidAssessment_Ricker_Cap",#working on this option
     #                 LRPmodel="BernLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
     #                useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.Multidim_Ricker_Cap",ps[pp]*100, sep=""),
     #                bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
     
     
  }




# =============================================================================================
# (4.1) Print and save logistic model fit diagnostics for any of the above runs
# =============================================================================================

# To look at individual model results:
#load("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/AnnualRetrospective/Bern.IndivRickerSurv_50/logisticFitDiagStats_2020.rda")
#load("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/AnnualRetrospective/Bern.IndivRickerSurvCap_50/logisticFitDiagStats_2020.rda")
#load("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/AnnualRetrospective/Bern.SPopAbundThreshST_50/logisticFitDiagStats_2020.rda")

# names(logisticDiagStats)
# logisticDiagStats$pBoxTidwell
# logisticDiagStats$DevResid
# max(abs(logisticDiagStats$DevResid))
# logisticDiagStats$p1
#  etc ...


# To save over all models and probabilities:

compile("TMB_Files/LRP_Logistic_Only_Bernoulli.cpp")
dyn.load(dynlib("TMB_Files/LRP_Logistic_Only_Bernoulli"))


year <- 2020
useBern_Logistic <- TRUE
ps <- c(0.5, 0.66, 0.9, 0.99)

# Create vectors of diagnostic stats to fill for each model
modelName<-NA
pBoxTidwell<-NA
maxDevResid<-NA
AR1.dev<-NA
minSampSize<-NA
sampSize<-NA
p_WaldB0<-NA
p_WaldB1<-NA
quasiR2<-NA
pGoodnessOfFit<-NA
hitRatio<-NA
hitRatio_LOO<-NA

for (pp in 1:length(ps)) {

  p<-ps[pp]
  
  # Model 1 =====================
  mm<-1
  modelName[mm]<-"IM"
  # Specify model fit to evaluate 
  BMmodel<-"SR_IndivRicker_Surv"
  RunName <- paste("Bern.IndivRickerSurv_",p*100, sep="")
  TMB_Inputs<-TMB_Inputs_IM
  outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")
  
  # Extract logistic regression data from .rda file
  load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
  pBoxTidwell[mm]<-round(logisticDiagStats$pBoxTidwell,3)
  maxDevResid[mm]<-round(max(abs(logisticDiagStats$DevResid)),2)
  AR1.dev[mm]<-round(logisticDiagStats$AR1.dev,2)
  minSampSize[mm]<-round(logisticDiagStats$minSampleSize,0)
  sampSize[mm]<-round(logisticDiagStats$sampleSize,0)
  p_WaldB0[mm]<-round(logisticDiagStats$signTable[logisticDiagStats$signTable$Param=="B_0","P.value"],3)
  p_WaldB1[mm]<-round(logisticDiagStats$signTable[logisticDiagStats$signTable$Param=="B_1","P.value"],3)
  quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
  pGoodnessOfFit[mm]<-round(logisticDiagStats$pDRT,3)
  hitRatio[mm]<-round(logisticDiagStats$hitRatio,2)
  # Run leave-one-out logistic model diagnostic:
  hitRatio_LOO[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                                RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)
  
  # Model 2 =====================
  mm<-2
  modelName[mm]<-"IM_Cap"
  # Specify model fit to evaluate
  BMmodel<-"SR_IndivRicker_SurvCap"
  RunName <- paste("Bern.IndivRickerSurvCap_",p*100, sep="")
  TMB_Inputs<-TMB_Inputs_IM_priorCap
  outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")
  
  # Extract logistic regression data from .rda file
  load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
  pBoxTidwell[mm]<-round(logisticDiagStats$pBoxTidwell,3)
  maxDevResid[mm]<-round(max(abs(logisticDiagStats$DevResid)),2)
  AR1.dev[mm]<-round(logisticDiagStats$AR1.dev,2)
  minSampSize[mm]<-round(logisticDiagStats$minSampleSize,0)
  sampSize[mm]<-round(logisticDiagStats$sampleSize,0)
  p_WaldB0[mm]<-round(logisticDiagStats$signTable[logisticDiagStats$signTable$Param=="B_0","P.value"],3)
  p_WaldB1[mm]<-round(logisticDiagStats$signTable[logisticDiagStats$signTable$Param=="B_1","P.value"],3)
  quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
  pGoodnessOfFit[mm]<-round(logisticDiagStats$pDRT,3)
  hitRatio[mm]<-round(logisticDiagStats$hitRatio,2)
  # Run leave-one-out logistic model diagnostic:
  hitRatio_LOO[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                                RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)
  
  
  # # Model 3 =====================
  mm<-3
  modelName[mm]<-"SubPop"
  # Specify model fit to evaluate
  RunName <- paste("Bern.SPopAbundThreshST_",p*100, sep="")
  TMB_Inputs<-TMB_Inputs_Subpop
  outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")
  
  # Extract logistic regression data from .rda file
  load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
  pBoxTidwell[mm]<-round(logisticDiagStats$pBoxTidwell,3)
  maxDevResid[mm]<-round(max(abs(logisticDiagStats$DevResid)),2)
  AR1.dev[mm]<-round(logisticDiagStats$AR1.dev,2)
  minSampSize[mm]<-round(logisticDiagStats$minSampleSize,0)
  sampSize[mm]<-round(logisticDiagStats$sampleSize,0)
  p_WaldB0[mm]<-round(logisticDiagStats$signTable[logisticDiagStats$signTable$Param=="B_0","P.value"],3)
  p_WaldB1[mm]<-round(logisticDiagStats$signTable[logisticDiagStats$signTable$Param=="B_1","P.value"],3)
  quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
  pGoodnessOfFit[mm]<-round(logisticDiagStats$pDRT,3)
  hitRatio[mm]<-round(logisticDiagStats$hitRatio,2)
  # Run leave-one-out logistic model diagnostic:
  hitRatio_LOO[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                                RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)
  
  
  # Summarize model fit diagnostics
  
  diagStats_byModel<-data.frame(model = modelName, pBoxTidwell, maxDevResid, AR1.dev, minSampSize, sampSize,
                                p_WaldB0, p_WaldB1, quasiR2, pGoodnessOfFit,                 
                                hitRatio, hitRatio_LOO)
  
  write.csv(diagStats_byModel,file=paste(cohoDir,"/DataOut/LogisticDiagStatsByModel_p=",p,"_",year,".csv",sep=""))
  
}


# ============================================================================================
# (4.2) Create plot to compare retrospective time series of LRP estimates from different methods
# ===========================================================================================

EscpDat.yy <- CoEscpDat %>% filter(yr <= 2020) 

modelFitList<-c("Bern.IndivRickerSurv",
                "Bern.IndivRickerSurvCap",
                "Bern.SPopAbundThreshST")

# Label names
L_Names<-c("Sgen-Ricker", "Sgen-priorCap", "IFCRT")  

if (useFrenchCaptions==TRUE) L_Names<-c("Ggén-Ricker", "Ggén-aprioriCap", "ERCFI")  

#pList<-c(50, 66, 90, 99)

pList<-c(50)

# Compare among probability thresholds:
#plotAnnualRetro_CompareProbs(Dat=EscpDat.yy, Names = modelFitList, pList=pList, L_Names = L_Names, outDir=cohoDir, useGenMean = T, genYrs = 3)

# Compare among methods:
plotAnnualRetro_CompareMethods(Dat=EscpDat.yy, Names = modelFitList, pList=pList, L_Names = L_Names, outDir=cohoDir, useGenMean = T, genYrs = 3, useFrenchCaptions)




# =================================================================================================================
# (5) EFFECT OF MISSING CUs: Run all available combinations of number of CUs 
# ================================================================================================================

nCUs<-c(5,4,3)
ps<-c(0.50, 0.66, 0.90,0.99)

for (pp in 1:length(ps)){

  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2020, BroodYrLag=4,
               genYrs=3, p = ps[pp], BMmodel="SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
               useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
               runLogisticDiag=F)
  
  # runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2020, BroodYrLag=4,
  #              genYrs=3, p = ps[pp], BMmodel="SR_HierRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
  #              useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bern.HierRickerSurv_",ps[pp]*100, sep=""),
  #              runLogisticDiag=F)
  
  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2020, BroodYrLag=4,
               genYrs=3, p = ps[pp], BMmodel="SR_IndivRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
               useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurvCap_",ps[pp]*100, sep=""),
               runLogisticDiag=F)
  
  # runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2020, BroodYrLag=4,
  #              genYrs=3, p = ps[pp], BMmodel="SR_HierRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
  #              useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bern.HierRickerSurvCap_",ps[pp]*100, sep=""),
  #              runLogisticDiag=F)
    
  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2020, BroodYrLag=4,
               genYrs=3, p = ps[pp], BMmodel="ThreshAbund_Subpop1000_ST", LRPmodel="BernLogistic", LRPfile="LRP_Logistic_Only", integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.SPopAbundThreshST_",ps[pp]*100, sep=""),
               runLogisticDiag=F)
  
}



# =======================================================================
# (5.1) Make Plots for Aggregate Status as a function of nCUS 
# =======================================================================

probThresh<-0.5

yearList<-2017:2020
plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_IndivRicker_Surv",p=probThresh, Dir=cohoDir,
                     inputPrefix="Bern.IndivRickerSurv_50",plotAveLine=TRUE, ymax=6)

yearList<-2017:2020
plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_IndivRicker_SurvCap",p=probThresh, Dir=cohoDir,
                     inputPrefix="Bern.IndivRickerSurvCap_50",plotAveLine=TRUE, ymax=5)

# yearList<-2015:2020
# plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv",p=probThresh, Dir=cohoDir,
#                      inputPrefix="Bern.HierRickerSurv_50",plotAveLine=TRUE, ymax=6)

yearList<-2017:2020
plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=probThresh, Dir=cohoDir,
                     inputPrefix="Bern.SPopAbundThreshST_50", plotAveLine=TRUE,ymax=8)



# =============================================================================================
# (6) Create figures with two or more probabilities on a logistic fit for specified years
# =============================================================================================

# Plot 1: Example plot for Methods Section. Show all 4 p-values without error bars on one plot as an example ======

BMmodel <- "SR_HierRicker_Surv"
TMB_Inputs<-TMB_Inputs_HM

year<-2019
BroodYrLag<-4
genYrs<-3

plotMultiP <-c(0.50, 0.66, 0.90, 0.99)
# Specify which years to use
EscpDat<-CoEscpDat
SRDat<-CoSRDat
Dat <- SRDat %>%  filter(BroodYear <= (year - BroodYrLag))
EscpDat.yy <- EscpDat %>% filter(yr <= (year - BroodYrLag))

Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")])) # This line is creating problems when I make CU_ID in data being read in

LRP_List<-list()

colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

for (pp in 1:length(plotMultiP)) {
  
  LRP_List[[pp]] <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = TRUE, 
                                   useGenMean = FALSE, genYrs = genYrs, p = plotMultiP[pp],  TMB_Inputs)
  
  
  if (pp == 1) {
    # Create plot using first p value
    example_LRP_plot <- ggplot(data=LRP_List[[1]]$Preds, mapping=aes(x=xx,y=fit)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr, x=xx), fill = "grey90") +
      geom_line(mapping=aes(x=xx, y=fit), col="black", size=1) +
      geom_line(mapping=aes(x=xx, y=upr), col="grey85") +
      geom_line(mapping=aes(x=xx, y=lwr), col="grey85") +
      geom_point(data=LRP_List[[1]]$Logistic_Data,
       aes(x = xx, y = yy)) +
      #geom_point(data=LRP_List[[1]]$Logistic_Data[which.max(LRP_List[[1]]$Logistic_Data$yr),], 
      #  aes(x = xx, y = yy), shape=4, stroke = 2) +
      xlab("Aggregate Spawner Abundance") + ylab("Pr(All CUs > Lower Benchmark)") +
      coord_cartesian(ylim=c(0,1)) +
      theme_classic()
  }
  
  # Add LRP estimates for each p-value =============================
  example_LRP_plot <- example_LRP_plot +
    geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
    geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5)
  #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) 
  
}

example_LRP_plot <- example_LRP_plot + geom_line(mapping=aes(x=xx, y=fit), col="black", size=1)

# Save plot
ggsave(paste(cohoDir, "/Figures/","methods-Example-LogisticLRP.png",sep=""), plot = example_LRP_plot,
       width = 4, height = 3, units = "in")  



# Plot 2: IM model Logistic fit in 2020 with all p-values ==================================================

BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs<-TMB_Inputs_IM

year<-2020
BroodYrLag<-4
genYrs<-3

plotMultiP <-c(0.50, 0.66, 0.90, 0.99)
# Specify which years to use
EscpDat<-CoEscpDat
SRDat<-CoSRDat
Dat <- SRDat %>%  filter(BroodYear <= (year - BroodYrLag))
EscpDat.yy <- EscpDat %>% filter(yr <= (year - BroodYrLag))

Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")])) # This line is creating problems when I make CU_ID in data being read in

LRP_List<-list()

colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

if (useFrenchCaptions==FALSE) {
  Xlab<-"Aggregate Spawner Abundance"
  Ylab<-"Pr(All CUs > Lower Benchmark)"
}

if (useFrenchCaptions==TRUE) {
  Xlab<-"Abondance agrégée des géniteurs"
  Ylab<-"Pr(toutes les UC > PRI)"
}

for (pp in 1:length(plotMultiP)) {
  
  LRP_List[[pp]] <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = TRUE, 
                                   useGenMean = FALSE, genYrs = genYrs, p = plotMultiP[pp],  TMB_Inputs)
  
  if (pp == 1) {
    # Create plot using first p value
    LRP_plot <- ggplot(data=LRP_List[[1]]$Preds, mapping=aes(x=xx,y=fit)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr, x=xx), fill = "grey90") +
      geom_line(mapping=aes(x=xx, y=fit), col="black", size=1) +
      geom_line(mapping=aes(x=xx, y=upr), col="grey85") +
      geom_line(mapping=aes(x=xx, y=lwr), col="grey85") +
      geom_point(data=LRP_List[[1]]$Logistic_Data, aes(x = xx, y = yy)) +
      geom_point(data=LRP_List[[1]]$Logistic_Data[which.max(LRP_List[[1]]$Logistic_Data$yr),],
       aes(x = xx, y = yy), shape=4, stroke = 2) +
      annotate("rect", xmin = LRP_List[[pp]]$LRP$lwr, xmax = LRP_List[[pp]]$LRP$upr, ymin=0, ymax=plotMultiP[pp], fill=colList[pp], alpha = .2) +
      geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
      geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5) +
      #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) +
      xlab(Xlab) + ylab(Ylab) +
      coord_cartesian(ylim=c(0,1)) +
      theme_classic()
  } else {
    # Add LRP estimates for each p-value =============================
    LRP_plot <- LRP_plot +
      geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
      geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5)
    #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) 
  }
  
}

if (useFrenchCaptions==FALSE) plotTitle1<-"Logistic:Sgen-Ricker"
if (useFrenchCaptions==TRUE) plotTitle1<-"Logistique : Ggén-Ricker"

LRP_plot1 <- LRP_plot + ggtitle(plotTitle1) +
      geom_line(mapping=aes(x=xx, y=fit), col="black", size=1)

if (useFrenchCaptions==FALSE) plotName<-"coho-IM2020-LogisticLRP.png"
if (useFrenchCaptions==TRUE) plotName<-"coho-IM2020-LogisticLRP-FN.png"

# Save plot
ggsave(paste(cohoDir, "/Figures/",plotName,sep=""), plot = LRP_plot1,
       width = 4, height = 3, units = "in")  



# Plot 3: IM.cap model Logistic fit in 2020 with all p-values ==================================================


BMmodel <- "SR_IndivRicker_SurvCap"
TMB_Inputs<-TMB_Inputs_IM_priorCap

year<-2020
BroodYrLag<-4
genYrs<-3

plotMultiP <-c(0.50, 0.66, 0.90, 0.99)
# Specify which years to use
EscpDat<-CoEscpDat
SRDat<-CoSRDat
Dat <- SRDat %>%  filter(BroodYear <= (year - BroodYrLag))
EscpDat.yy <- EscpDat %>% filter(yr <= (year - BroodYrLag))

Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")])) # This line is creating problems when I make CU_ID in data being read in

LRP_List<-list()

colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

for (pp in 1:length(plotMultiP)) {
  
  LRP_List[[pp]] <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = TRUE, 
                                   useGenMean = FALSE, genYrs = genYrs, p = plotMultiP[pp],  TMB_Inputs)
  
  if (pp == 1) {
    # Create plot using first p value
    LRP_plot <- ggplot(data=LRP_List[[1]]$Preds, mapping=aes(x=xx,y=fit)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr, x=xx), fill = "grey90") +
      geom_line(mapping=aes(x=xx, y=fit), col="black", size=1) +
      geom_line(mapping=aes(x=xx, y=upr), col="grey85") +
      geom_line(mapping=aes(x=xx, y=lwr), col="grey85") +
      geom_point(data=LRP_List[[1]]$Logistic_Data, aes(x = xx, y = yy)) +
      geom_point(data=LRP_List[[1]]$Logistic_Data[which.max(LRP_List[[1]]$Logistic_Data$yr),],
       aes(x = xx, y = yy), shape=4, stroke = 2) +
      annotate("rect", xmin = LRP_List[[pp]]$LRP$lwr, xmax = LRP_List[[pp]]$LRP$upr, ymin=0, ymax=plotMultiP[pp], fill=colList[pp], alpha = .2) +
      geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
      geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5) +
      #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) +
      xlab("Aggregate Spawner Abundance") + ylab("Pr(All CUs > Lower Benchmark)") +
      coord_cartesian(ylim=c(0,1)) +
      theme_classic()
  } else {
    # Add LRP estimates for each p-value =============================
    LRP_plot <- LRP_plot +
      geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
      geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5)
    #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) 
  }
  
}

if (useFrenchCaptions==FALSE) plotTitle2<-"Logistic:Sgen-priorCap"
if (useFrenchCaptions==TRUE) plotTitle2<-"Logistique : Ggén-aprioriCap"

LRP_plot2 <- LRP_plot + ggtitle(plotTitle2) +
geom_line(mapping=aes(x=xx, y=fit), col="black", size=1)

# Save plot

if (useFrenchCaptions==FALSE) plotName<-"coho-IMCap2020-LogisticLRP.png"
if (useFrenchCaptions==TRUE) plotName<-"coho-IMCap2020-LogisticLRP-FN.png"

ggsave(paste(cohoDir, "/Figures/",plotName,sep=""), plot = LRP_plot2,
       width = 4, height = 3, units = "in")  





# Plot 4: Subpopulation threshold model Logistic fit in 2020 with all p-values ==================================================

BMmodel <- "ThreshAbund_Subpop1000_ST"
TMB_Inputs<-TMB_Inputs_Subpop

year<-2020
BroodYrLag<-4
genYrs<-3

plotMultiP <-c(0.50, 0.66, 0.90, 0.99)

# Specify which years to use
EscpDat<-CoEscpDat_bySubpop
EscpDat.yy <- EscpDat %>% filter(yr <= (year - BroodYrLag))

Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")])) # This line is creating problems when I make CU_ID in data being read in

# Calculate geometric means for subpopulation escapements (where, subpopulation is sum of tributaries) 
EscpDat.cu <- EscpDat %>% group_by(CU_Name, CU_ID, Subpop_Name, yr)  %>% summarise(SubpopEscp=sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(SubpopEscp, genYrs, gm_mean, fill = NA, align="right"))
# Add a column indicating whether geometric mean escapement is > 1000 fish threshold for each subpopulation 
#Above1000<-ifelse(EscpDat.cu$Gen_Mean >= 1000, 1, 0)
Above1000<-ifelse(EscpDat.cu$SubpopEscp >= 1000, 1, 0)
EscpDat.cu<-EscpDat.cu %>% add_column(Above1000)

# Calculate the number of subpopulations that are above 1000 fish threshold in each CU
LBM_status_byCU <- EscpDat.cu %>% group_by(CU_Name, CU_ID, yr) %>% 
  summarise(Escp=sum(SubpopEscp), N = length(unique(Subpop_Name)),N_grThresh=sum(Above1000))

# Add a column indicating whether >= 50% of subpopulations were above 1000 fish threshold in each CU (1 = yes, 0 = no) (short-term recovery goal)
HalfGrThresh<-ifelse(LBM_status_byCU$N_grThresh >=  ceiling(LBM_status_byCU$N/2),1,0)
LBM_status_byCU <- LBM_status_byCU %>% add_column(AboveBenchmark=HalfGrThresh)


LRP_List<-list()

colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")


for (pp in 1:length(plotMultiP)) {
  
  LRP_List[[pp]] <- Run_LRP(Dat=LBM_status_byCU, Mod = "LRP_Logistic_Only", useBern_Logistic = TRUE, 
                            useGenMean = FALSE, genYrs = genYrs, p = plotMultiP[pp],  TMB_Inputs)
  
  if (pp == 1) {
    # Create plot using first p value
    LRP_plot <- ggplot(data=LRP_List[[1]]$Preds, mapping=aes(x=xx,y=fit)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr, x=xx), fill = "grey90") +
      geom_line(mapping=aes(x=xx, y=fit), col="black", size=1) +
      geom_line(mapping=aes(x=xx, y=upr), col="grey85") +
      geom_line(mapping=aes(x=xx, y=lwr), col="grey85") +
      geom_point(data=LRP_List[[1]]$Logistic_Data, aes(x = xx, y = yy)) +
      geom_point(data=LRP_List[[1]]$Logistic_Data[which.max(LRP_List[[1]]$Logistic_Data$yr),],
       aes(x = xx, y = yy), shape=4, stroke = 2) +
      annotate("rect", xmin = LRP_List[[pp]]$LRP$lwr, xmax = LRP_List[[pp]]$LRP$upr, ymin=0, ymax=plotMultiP[pp], fill=colList[pp], alpha = .2) +
      geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
      geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5) +
      #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) +
      xlab(Xlab) + ylab(Ylab) +
      coord_cartesian(ylim=c(0,1)) +
      theme_classic()
  } else {
    # Add LRP estimates for each p-value =============================
    LRP_plot <- LRP_plot +
      geom_line(dat=data.frame(x=rep(LRP_List[[pp]]$LRP$fit,2),y=c(0,plotMultiP[pp])),aes(x=x,y=y), color=colList[pp],size=1) +
      geom_hline(yintercept= plotMultiP[pp], linetype="dotted", color="black", size = 0.5)
    #geom_line(dat=data.frame(x=c(0,LRP_List[[pp]]$LRP$fit),y=rep(plotMultiP[pp],2)),aes(x=x,y=y), color=colList[pp], linetype = "dotted",size=0.8) 
  }
  
}

if (useFrenchCaptions==FALSE) plotTitle3<-"Logistic:IFCRT LRP"
if (useFrenchCaptions==TRUE) plotTitle3<-"Logistique : ERCFI PRL"

LRP_plot3 <- LRP_plot + ggtitle(plotTitle3) +
geom_line(mapping=aes(x=xx, y=fit), col="black", size=1)


if (useFrenchCaptions==FALSE) plotName<-"coho-ThreshAb2020-LogisticLRP.png"
if (useFrenchCaptions==TRUE) plotName<-"coho-ThreshAb2020-LogisticLRP-FN.png"

# Save plot
ggsave(paste(cohoDir, "/Figures/",plotName,sep=""), plot = LRP_plot3,
       width = 4, height = 3, units = "in")  



LRP_plotall<-plot_grid(LRP_plot1,LRP_plot2,LRP_plot3)

if (useFrenchCaptions==FALSE) plotName<-"coho-all-LogisticLRP.png"
if (useFrenchCaptions==TRUE) plotName<-"coho-all-LogisticLRP-FN.png"

ggsave(paste(cohoDir, "/Figures/",plotName,sep=""), plot = LRP_plotall,
       width = 8, height = 7)  

# =================================================================================================================
# (6b) NOT YET WORKING: Create table of 2020 estimates
# ================================================================================================================

plotMultiP <- c(0.50,0.66,0.90,0.99) 
methodList<-c("Bern.IndivRickerSurv", "Bern.IndivRickerSurvCap", "Bern.SPopAbundThreshST")

for (pp in 1:length(plotMultiP)) {
  
  
  
}




# ============================================================================================
# (7) STILL TO BE TESTED:Create plot to compare status as a function of nCUs for different LRP options
# ===========================================================================================
  

modelFitList<-c("Bern.IndivRickerSurv",
                "Bern.HierRickerSurv",
                "Bern.IndivRickerSurvCap",
                "Bern.HierRickerSurvCap",
                "Bern.SPopAbundThreshST")


plotAggStatus_byNCUs_Compare(estYear=2020, nCUList=c(5,4,3), Names=modelFitList, p=0.50, Dir=cohoDir,plotAveLine=TRUE)  


# Notes on comparing LW results run (2020-10-26) with results on dropbox https://www.dropbox.com/sh/otiku88jc2eu8cx/AACagLQd85blX7jc8-yBYs4Na?dl=0 
#   1. DataOut/ModelFits/AllEsts_Hier_Ricker_Surv.csv: 
#         sigmaA estimate is 0.321763 vs. -1.133939234 on dropbox
#         other values look the same
#   2. DataOut/AnnualRetrospective/annualRetro_LRPs.csv:
#         values are identical




