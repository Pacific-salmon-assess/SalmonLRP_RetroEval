
#useFrenchCaptions<-FALSE
useFrenchCaptions<-TRUE

biasCorrectEst<-TRUE

# Make simplified version of TMB run code for demonstration
library(TMB)
library(dplyr)
library(zoo)
library(ggplot2)
library(gridExtra)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

# Create directory to save model fitting estimates
outDir <- paste(cohoDir, "/DataOut", sep="")
if (file.exists(outDir) == FALSE){
  dir.create(outDir)
}

# Create directory to save model fitting estimates
outDir <- paste(outDir, "/ModelFits", sep="")
if (file.exists(outDir) == FALSE){
  dir.create(outDir)
}

# Create a directory to save model fitting figures
figDir <- paste(cohoDir, "/Figures/ModelFits", sep="")
if (file.exists(figDir) == FALSE){
  dir.create(figDir)
}


setwd(codeDir)

compile("TMB_Files/SR_HierRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_noLRP"))

compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))

compile("TMB_Files/SR_HierRicker_SurvCap_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_noLRP"))

compile("TMB_Files/SR_IndivRicker_SurvCap_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_noLRP"))


sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
}

sourceAll()

# Read-in data
setwd(cohoDir)

CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
  # Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
  colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
  colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
  colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
  colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")

# CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")
# 
#   CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
#   tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
#   CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)

setwd(codeDir)
  

# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
#CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)

# Prep data frame
CoSRDat$yr_num <- CoSRDat$BroodYear - min(CoSRDat$BroodYear)
CoSRDat$CU_ID <- group_indices(CoSRDat, CU_ID) - 1
# 
# CoEscpDat$yr_num <- group_indices(CoEscpDat, yr) - 1
# CoEscpDat<- CoEscpDat %>% right_join(unique(CoSRDat[,c("CU_ID", "CU_Name")]))


SRDat<-CoSRDat
#EscDat <- CoEscpDat
#Bern_Logistic <- T
useGenMean <- F
genYrs <- 3
p<-0.5

Scale<-1000


# Create vector of spawner abundances to plot predicted SR relationships from 
Pred_Spwn <- rep(seq(0,30000/Scale,length=301), 5)
Pred_Spwn_CU <- c(rep(0,301), rep(1,301), rep(2,301), rep(3,301), rep(4,301))


# *************************************************************************************
# Fit Indiv_Ricker_Surv model
# ********************************************************************************

TMB_Inputs<- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                      extra_eval_iter=FALSE,biasCorrect=biasCorrectEst)

BMmodel <- "SR_IndivRicker_Surv"

Mod <- paste(BMmodel,"noLRP",sep="_")

# Set-up for call to TMB
Scale <- TMB_Inputs$Scale

data <- list()
data$Bayes <- 1 # Indicate that this is an MLE, not Bayesian
if (is.null(biasCorrectEst)) { 
  data$BiasCorrect <- 0 } else {
    # Indicate whether log-normal bias correction should be applied in estimation
    if (biasCorrectEst == TRUE) data$BiasCorrect <-1
    if (biasCorrectEst == FALSE) data$BiasCorrect <-0
  }
data$S <- SRDat$Spawners/Scale
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Sgen_sig <- TMB_Inputs$Sgen_sig # set variance to be used for likelihood for estimating Sgen

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
data$Tau_dist <- TMB_Inputs$Tau_dist
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)
#muSurv <- SRDat %>% group_by(CU_ID) %>%
#  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
#data$muLSurv <- log(muSurv$muSurv)

# Base mu survival on mean of age 3 survival (not weighted by historic age at return)
muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
data$muLSurv <- muLSurv$muLSurv

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean
data$gamma_sig <- TMB_Inputs$gamma_sig

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  data$logMuA_mean <- TMB_Inputs$logMuA_mean
  data$logMuA_sig <- TMB_Inputs$logMuA_sig
  data$Tau_A_dist <- TMB_Inputs$Tau_A_dist
  param$logMuA <- TMB_Inputs$logA_Start
  param$logSigmaA <- 1
}

if(BMmodel %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv")) {
  param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
}

if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
  param$cap <- TMB_Inputs$cap_mean
  data$cap_mean<-TMB_Inputs$cap_mean
  data$cap_sig<-TMB_Inputs$cap_sig
}

# range of spawner abundance to predict recruitment from
#data$Pred_Spwn <- rep(seq(0,max(data$S)*1.1,length=100), N_Stocks) # vectors of spawner abundance to use for predicted recruits, one vector of length 100 for each stock


data$Pred_Spwn<-Pred_Spwn
data$stk_predS<-Pred_Spwn_CU


# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks))) # Fix Sgen

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
  #  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
  #obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE)
}

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
# 

# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
if (BMmodel == "SR_HierRicker_Surv" | BMmodel == "SR_HierRicker_SurvCap") {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE)
}


# Create upper bounds vector that is same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.01)
lower[names(lower) =="cap"] <- 0
lower<-unname(lower)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower)

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)



# put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests$Mod <- Mod
All_Ests$CU_ID[!(All_Ests$Param %in% c("logMuA", "logSigmaA", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests$CU_ID[All_Ests$Param=="Rec_Preds"]<-data$stk_predS
All_Ests <- left_join(All_Ests, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests$Estimate[All_Ests$Param == "logSigma"] <- exp(All_Ests$Estimate[All_Ests$Param == "logSigma"] )
All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
All_Ests$Estimate[All_Ests$Param == "logSigmaA"] <- exp(All_Ests$Estimate[All_Ests$Param == "logSigmaA"] )
All_Ests$Param[All_Ests$Param == "logSigmaA"] <- "sigmaA"
All_Ests$Estimate[All_Ests$Param == "logB"] <- exp(All_Ests$Estimate[All_Ests$Param == "logB"] )
All_Ests$Param[All_Ests$Param == "logB"] <- "B"
All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
All_Ests[All_Ests$Param %in% c("Sgen", "Smsy","SRep"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "Smsy", "SRep")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds_Rec <- All_Ests %>% filter(Param == "Rec_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Rec_Preds"))) 

write.csv(All_Ests,paste(outDir,"/AllEsts_Indiv_Ricker_Surv.csv", sep=""))







# *************************************************************************************
# Fit Indiv_Ricker_Surv_PriorCap model
# ********************************************************************************

# Extract capacity estimates from Indiv_Ricker_Surv model fit to use as priors on carrying capacity
CU_Names <- unique(SRDat[, "CU_Name"])
cap<-rep(NA,N_Stocks)

for (i in 1:N_Stocks) {
 logA<-All_Ests[All_Ests$Param =="logA" & All_Ests$CU_Name == CU_Names[i],"Estimate"]
 A<-All_Ests[All_Ests$Param =="A" & All_Ests$CU_Name == CU_Names[i],"Estimate"]
 B<-All_Ests[All_Ests$Param =="B" & All_Ests$CU_Name == CU_Names[i],"Estimate"] * Scale
 cap[i] <- log(logA)/B
}

cap_priorMean<- cap*1.40

print("Prior means for capacity in Indiv model, by CU:")
print(cap_priorMean)

pngNameSrep <- "coho-SrepPriorDist.png"
if (useFrenchCaptions == TRUE) pngNameSrep <- "coho-SrepPriorDist-FN.png"

png(paste(figDir, "/",pngNameSrep, sep=""), width=480, height=300)
# Plot prior distributions, by CU
xx<-seq(1,30,length=1000)

# manually set x-axes for now by specifying min and max for each CU
x.plotMin<-c(5000, 800, 5000, 13000, 8000)
x.plotMax<-c(18000, 12500, 18000, 28000,22000)
par(mfrow=c(2,3), mar=c(4,3,2,2), oma=c(1,1,1,1))
# Loop over CUs to plot
for (i in 1:5) {
  if (CU_Names[i] == "Middle_Fraser") plotName<-"Middle Fraser"
  if (CU_Names[i] == "Fraser_Canyon") plotName<-"Fraser_Canyon"
  if (CU_Names[i] == "North_Thompson") plotName<-"North Thompson"
  if (CU_Names[i] == "South_Thompson") plotName<-"South Thompson"
  if (CU_Names[i] == "Lower_Thompson") plotName<-"Lower Thompson"
  
  if (useFrenchCaptions == TRUE) {
    if (CU_Names[i] == "Middle_Fraser") plotName<-"Moyen Fraser"
    if (CU_Names[i] == "Fraser_Canyon") plotName<-"Canyon du Fraser"
    if (CU_Names[i] == "North_Thompson") plotName<-"Thompson Nord"
    if (CU_Names[i] == "South_Thompson") plotName<-"Thompson Sud"
    if (CU_Names[i] == "Lower_Thompson") plotName<-"Thompson inférieure"
  }
    
  
  Xlab<-"SRep"
  if(useFrenchCaptions == TRUE) Xlab<-"GRem"
  
  yy<-dnorm(xx,mean=cap_priorMean[i],sqrt(2))
  plot(xx*1000,yy*1000, typ="l", lwd=2, 
       ylab="", xlab=Xlab, xlim=c(x.plotMin[i], x.plotMax[i]), axes=F, cex.lab=1.5)
  abline(v=cap[i]*1000,lty=2, col="red")
  title(main=plotName, cex.main=1.8)
  axis(side = 1,labels=T, cex.axis=1.5)
  axis(side = 2,labels=F, cex.axis=1.5)
  box()
}
dev.off()



# Compile TMB model

# Inputs we are using:
TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                      logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                      cap_mean=cap_priorMean, cap_sig=sqrt(2),
                      extra_eval_iter=FALSE,biasCorrect=biasCorrectEst)


BMmodel <- "SR_IndivRicker_SurvCap"


Mod <- paste(BMmodel,"noLRP",sep="_")

# Set-up for call to TMB
Scale <- TMB_Inputs$Scale

data <- list()
data$Bayes <- 1 # Indicate that this is an MLE, not Bayesian
if (is.null(biasCorrectEst)) { 
  data$BiasCorrect <- 0 } else {
    # Indicate whether log-normal bias correction should be applied in estimation
    if (biasCorrectEst == TRUE) data$BiasCorrect <-1
    if (biasCorrectEst == FALSE) data$BiasCorrect <-0
  }
data$S <- SRDat$Spawners/Scale
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Sgen_sig <- TMB_Inputs$Sgen_sig # set variance to be used for likelihood for estimating Sgen

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
data$Tau_dist <- TMB_Inputs$Tau_dist
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)
#muSurv <- SRDat %>% group_by(CU_ID) %>%
#  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
#data$muLSurv <- log(muSurv$muSurv)

# Base mu survival on mean of age 3 survival (not weighted by historic age at return)
muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
data$muLSurv <- muLSurv$muLSurv

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean
data$gamma_sig <- TMB_Inputs$gamma_sig

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  data$logMuA_mean <- TMB_Inputs$logMuA_mean
  data$logMuA_sig <- TMB_Inputs$logMuA_sig
  data$Tau_A_dist <- TMB_Inputs$Tau_A_dist
  param$logMuA <- TMB_Inputs$logA_Start
  param$logSigmaA <- 1
}

if(BMmodel %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv")) {
  param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
}

if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
  param$cap <- TMB_Inputs$cap_mean
  data$cap_mean<-TMB_Inputs$cap_mean
  data$cap_sig<-TMB_Inputs$cap_sig
}

# range of spawner abundance to predict recruitment from
data$Pred_Spwn<-Pred_Spwn
data$stk_predS<-Pred_Spwn_CU


# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks))) # Fix Sgen

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
  #  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
  #obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE)
}

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
 

# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
if (BMmodel == "SR_HierRicker_Surv" | BMmodel == "SR_HierRicker_SurvCap") {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE)
}


# Create upper bounds vector that is same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.01)
lower[names(lower) =="cap"] <- 0
lower<-unname(lower)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower)

All_Ests_cap <- data.frame(summary(sdreport(obj)))
All_Ests_cap$Param <- row.names(All_Ests_cap)

# put together readable data frame of values
All_Ests_cap$Param <- sapply(All_Ests_cap$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests_cap$Mod <- Mod
All_Ests_cap$CU_ID[!(All_Ests_cap$Param %in% c("logMuA", "logSigmaA","gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests_cap$CU_ID[All_Ests_cap$Param=="Rec_Preds"]<-data$stk_predS
All_Ests_cap <- left_join(All_Ests_cap, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests_cap$Estimate[All_Ests_cap$Param == "logSigma"] <- exp(All_Ests_cap$Estimate[All_Ests_cap$Param == "logSigma"] )
All_Ests_cap$Param[All_Ests_cap$Param == "logSigma"] <- "sigma"
All_Ests_cap$Estimate[All_Ests_cap$Param == "logSigmaA"] <- exp(All_Ests$Estimate[All_Ests$Param == "logSigmaA"] )
All_Ests_cap$Param[All_Ests_cap$Param == "logSigmaA"] <- "sigmaA"
All_Ests_cap$Param[All_Ests_cap$Param == "SMSY"] <- "Smsy"
All_Ests_cap[All_Ests_cap$Param == "B",] <- All_Ests_cap %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests_cap[All_Ests_cap$Param %in% c("Sgen", "Smsy", "SRep", "cap"), ] <-  All_Ests_cap %>% filter(Param %in% c("Sgen", "Smsy", "SRep","cap")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds_Rec_cap <- All_Ests_cap %>% filter(Param == "Rec_Preds")
All_Ests_cap <- All_Ests_cap %>% filter(!(Param %in% c( "logSgen", "Rec_Preds"))) 


write.csv(All_Ests_cap,paste(outDir,"/AllEsts_Indiv_Ricker_Surv_priorCap.csv", sep=""))




# *************************************************************************************
# Fit Hier_Ricker_Surv model
# ********************************************************************************

# What we're using:
TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                   logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1, 
                   extra_eval_iter=FALSE,biasCorrect=biasCorrectEst)


BMmodel <- "SR_HierRicker_Surv"

# Removing most recent data year due to non-convergence using data up to return year 2020
#SRDat <- CoSRDat %>% filter(BroodYear <= 2014)

#SRDat <- CoSRDat %>% filter(BroodYear <= 2015)

SRDat <- CoSRDat

Mod <- paste(BMmodel,"noLRP",sep="_")

# Set-up for call to TMB
Scale <- TMB_Inputs$Scale

data <- list()
data$Bayes <- 1 # Indicate that this is an MLE, not Bayesian
if (is.null(biasCorrectEst)) { 
  data$BiasCorrect <- 0 } else {
    # Indicate whether log-normal bias correction should be applied in estimation
    if (biasCorrectEst == TRUE) data$BiasCorrect <-1
    if (biasCorrectEst == FALSE) data$BiasCorrect <-0
  }
data$S <- SRDat$Spawners/Scale
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Sgen_sig <- TMB_Inputs$Sgen_sig # set variance to be used for likelihood for estimating Sgen

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
data$Tau_dist <- TMB_Inputs$Tau_dist
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)
muSurv <- SRDat %>% group_by(CU_ID) %>%
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

# Base mu survival on mean of age 3 survival (not weighted by historic age at return)
#muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
#data$muLSurv <- muLSurv$muLSurv

data$gamma_mean <- TMB_Inputs$gamma_mean
data$gamma_sig <- TMB_Inputs$gamma_sig

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  data$logMuA_mean <- TMB_Inputs$logMuA_mean
  data$logMuA_sig <- TMB_Inputs$logMuA_sig
  data$Tau_A_dist <- TMB_Inputs$Tau_A_dist
  param$logMuA <- TMB_Inputs$logA_Start
  param$logSigmaA <- TMB_Inputs$logMuA_sig
}

if(BMmodel %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv")) {
  param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
}

if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
  param$cap <- TMB_Inputs$cap_mean
  data$cap_mean<-TMB_Inputs$cap_mean
  data$cap_sig<-TMB_Inputs$cap_sig
}

data$Pred_Spwn<-Pred_Spwn
data$stk_predS<-Pred_Spwn_CU


# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks))) # Fix Sgen

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
  #  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
  #obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE)
}

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
if (BMmodel == "SR_HierRicker_Surv" | BMmodel == "SR_HierRicker_SurvCap") {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE)
}


# Create upper bounds vector that is same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.01)
lower[names(lower) =="cap"] <- 0
lower<-unname(lower)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower)

All_Ests_HM <- data.frame(summary(sdreport(obj)))
All_Ests_HM$Param <- row.names(All_Ests_HM)



# put together readable data frame of values
All_Ests_HM$Param <- sapply(All_Ests_HM$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests_HM$Mod <- Mod
All_Ests_HM$CU_ID[!(All_Ests_HM$Param %in% c("logMuA", "logSigmaA", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests_HM$CU_ID[All_Ests_HM$Param=="Rec_Preds"]<-data$stk_predS
All_Ests_HM <- left_join(All_Ests_HM, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests_HM$Estimate[All_Ests_HM$Param == "logSigma"] <- exp(All_Ests_HM$Estimate[All_Ests_HM$Param == "logSigma"] )
All_Ests_HM$Param[All_Ests_HM$Param == "logSigma"] <- "sigma"
All_Ests_HM$Estimate[All_Ests_HM$Param == "logSigmaA"] <- exp(All_Ests_HM$Estimate[All_Ests_HM$Param == "logSigmaA"] )
All_Ests_HM$Param[All_Ests_HM$Param == "logSigmaA"] <- "sigmaA"
All_Ests_HM$Estimate[All_Ests_HM$Param == "logB"] <- exp(All_Ests_HM$Estimate[All_Ests_HM$Param == "logB"] )
All_Ests_HM$Param[All_Ests_HM$Param == "logB"] <- "B"
All_Ests_HM[All_Ests_HM$Param == "B",] <- All_Ests_HM %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests_HM$Param[All_Ests_HM$Param == "SMSY"] <- "Smsy"
All_Ests_HM[All_Ests_HM$Param %in% c("Sgen", "Smsy","SRep"), ] <-  All_Ests_HM %>% filter(Param %in% c("Sgen", "Smsy", "SRep")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds_Rec_HM <- All_Ests_HM %>% filter(Param == "Rec_Preds")
All_Ests_HM <- All_Ests_HM %>% filter(!(Param %in% c( "logSgen", "Rec_Preds"))) 

write.csv(All_Ests_HM,paste(outDir,"/AllEsts_Hier_Ricker_Surv.csv", sep=""))





# *************************************************************************************
# Fit Hier_Ricker_Surv model with prior cap
# ********************************************************************************


# Extract capacity estimates from Indiv_Ricker_Surv model fit to use as priors on carrying capacity
CU_Names <- unique(SRDat[, "CU_Name"])
cap<-rep(NA,N_Stocks)

for (i in 1:N_Stocks) {
  logA<-All_Ests_HM[All_Ests_HM$Param =="logA" & All_Ests_HM$CU_Name == CU_Names[i],"Estimate"]
  A<-All_Ests_HM[All_Ests_HM$Param =="A" & All_Ests_HM$CU_Name == CU_Names[i],"Estimate"]
  B<-All_Ests_HM[All_Ests_HM$Param =="B" & All_Ests_HM$CU_Name == CU_Names[i],"Estimate"] * Scale
  cap[i] <- log(logA)/B
}

cap_priorMean<- cap*1.40

print("Prior means for capacity in Hier model, by CU:")
print(cap_priorMean)



png(paste(figDir, "/SrepPriorDist_HM_priorCap.png", sep=""), width=480, height=300)
# Plot prior distributions, by CU
xx<-seq(1,30,length=1000)
# manually set x-axes for now by specifying min and max for each CU
x.plotMin<-c(5000, 800, 5000, 13000, 8000)
x.plotMax<-c(18000, 12500, 18000, 28000,22000)
par(mfrow=c(2,3), mar=c(4,3,1,1))
# Loop over CUs to plot
for (i in 1:5) {
  yy<-dnorm(xx,mean=cap_priorMean[i],sqrt(2))
  plot(xx*1000,yy*1000, main=CU_Names[i], typ="l", lwd=2, 
       ylab="", xlab="SRep", xlim=c(x.plotMin[i], x.plotMax[i]), axes=F)
  abline(v=cap[i]*1000,lty=2, col="red")
  axis(side = 1,labels=T)
  axis(side = 2,labels=F)
  box()
}
dev.off()

# Compile TMB model

# Inputs we are using:
TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                     logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                     gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                     cap_mean=cap_priorMean, cap_sig=sqrt(2),
                     extra_eval_iter=FALSE,biasCorrect=biasCorrectEst)
  

BMmodel <- "SR_HierRicker_SurvCap"


Mod <- paste(BMmodel,"noLRP",sep="_")

# Set-up for call to TMB
Scale <- TMB_Inputs$Scale

data <- list()
data$Bayes <- 1 # Indicate that this is an MLE, not Bayesian
if (is.null(biasCorrectEst)) { 
  data$BiasCorrect <- 0 } else {
    # Indicate whether log-normal bias correction should be applied in estimation
    if (biasCorrectEst == TRUE) data$BiasCorrect <-1
    if (biasCorrectEst == FALSE) data$BiasCorrect <-0
  }
data$S <- SRDat$Spawners/Scale
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Sgen_sig <- TMB_Inputs$Sgen_sig # set variance to be used for likelihood for estimating Sgen

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
data$Tau_dist <- TMB_Inputs$Tau_dist
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)
muSurv <- SRDat %>% group_by(CU_ID) %>%
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

# Base mu survival on mean of age 3 survival (not weighted by historic age at return)
#muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
#data$muLSurv <- muLSurv$muLSurv

data$gamma_mean <- TMB_Inputs$gamma_mean
data$gamma_sig <- TMB_Inputs$gamma_sig

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  data$logMuA_mean <- TMB_Inputs$logMuA_mean
  data$logMuA_sig <- TMB_Inputs$logMuA_sig
  data$Tau_A_dist <- TMB_Inputs$Tau_A_dist
  param$logMuA <- TMB_Inputs$logA_Start
  param$logSigmaA <- TMB_Inputs$logMuA_sig
}

if(BMmodel %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv")) {
  param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
}

if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
  param$cap <- TMB_Inputs$cap_mean
  data$cap_mean<-TMB_Inputs$cap_mean
  data$cap_sig<-TMB_Inputs$cap_sig
}

data$Pred_Spwn<-Pred_Spwn
data$stk_predS<-Pred_Spwn_CU


# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks))) # Fix Sgen

if (BMmodel %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")) {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
  #  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
  #obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE)
}

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
if (BMmodel == "SR_HierRicker_Surv" | BMmodel == "SR_HierRicker_SurvCap") {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE)
}


# Create upper bounds vector that is same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.01)
lower[names(lower) =="cap"] <- 0
lower<-unname(lower)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower)

All_Ests_HMcap <- data.frame(summary(sdreport(obj)))
All_Ests_HMcap$Param <- row.names(All_Ests_HMcap)



# put together readable data frame of values
All_Ests_HMcap$Param <- sapply(All_Ests_HMcap$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests_HMcap$Mod <- Mod
All_Ests_HMcap$CU_ID[!(All_Ests_HMcap$Param %in% c("logMuA", "logSigmaA", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests_HMcap$CU_ID[All_Ests_HMcap$Param=="Rec_Preds"]<-data$stk_predS
All_Ests_HMcap <- left_join(All_Ests_HMcap, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests_HMcap$Estimate[All_Ests_HMcap$Param == "logSigma"] <- exp(All_Ests_HMcap$Estimate[All_Ests_HMcap$Param == "logSigma"] )
All_Ests_HMcap$Param[All_Ests_HMcap$Param == "logSigma"] <- "sigma"
All_Ests_HMcap$Estimate[All_Ests_HMcap$Param == "logSigmaA"] <- exp(All_Ests_HMcap$Estimate[All_Ests_HMcap$Param == "logSigmaA"] )
All_Ests_HMcap$Param[All_Ests_HMcap$Param == "logSigmaA"] <- "sigmaA"
All_Ests_HMcap$Estimate[All_Ests_HMcap$Param == "logB"] <- exp(All_Ests_HMcap$Estimate[All_Ests_HMcap$Param == "logB"] )
All_Ests_HMcap$Param[All_Ests_HMcap$Param == "logB"] <- "B"
All_Ests_HMcap[All_Ests_HMcap$Param == "B",] <- All_Ests_HMcap %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests_HMcap$Param[All_Ests_HMcap$Param == "SMSY"] <- "Smsy"
All_Ests_HMcap[All_Ests_HMcap$Param %in% c("Sgen", "Smsy","SRep"), ] <-  All_Ests_HMcap %>% filter(Param %in% c("Sgen", "Smsy", "SRep")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds_Rec_HMcap <- All_Ests_HMcap %>% filter(Param == "Rec_Preds")
All_Ests_HMcap <- All_Ests_HMcap %>% filter(!(Param %in% c( "logSgen", "Rec_Preds"))) 

write.csv(All_Ests_HMcap,paste(outDir,"/AllEsts_Hiercap_Ricker_Surv.csv", sep=""))






















# ========================================================================================================
# Create plots of fit SR curves
# ================================================================================================

CU_list<-unique(SRDat[, "CU_Name"])
CUID_list<-unique(SRDat[, "CU_ID"])
nCUs<-length(CU_list)


# Create dataframes with predictions to plot
for (i in 1:nCUs) {
 
#   # create data frame for hier_surv model
#   plotDat.CU<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
#                     Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
#                     Pred_Rec=Preds_Rec[Preds_Rec$CU_ID == CUID_list[i], "Estimate"] * Scale,
#                     Pred_Rec_SE=Preds_Rec[Preds_Rec$CU_ID == CUID_list[i], "Std..Error"] * Scale)
#   
#     if (i == 1) plotDat<-plotDat.CU
#     if (i > 1) plotDat<-rbind(plotDat, plotDat.CU)
#     
#     # create data frame for hier_surv_prior_cap model
#     plotDat.CU_cap<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
#                            Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
#                            Pred_Rec=Preds_Rec_cap[Preds_Rec_cap$CU_ID == CUID_list[i], "Estimate"] * Scale,
#                            Pred_Rec_SE=Preds_Rec_cap[Preds_Rec_cap$CU_ID == CUID_list[i], "Std..Error"] * Scale)
#     
#     if (i == 1) plotDat_cap<-plotDat.CU_cap
#     if (i > 1) plotDat_cap<-rbind(plotDat_cap, plotDat.CU_cap)
  
    # create data frame for IM_surv model
    plotDat.CU_IM<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
                           Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
                           Pred_Rec=Preds_Rec[Preds_Rec$CU_ID == CUID_list[i], "Estimate"] * Scale,
                           Pred_Rec_SE=Preds_Rec[Preds_Rec$CU_ID == CUID_list[i], "Std..Error"] * Scale)
    
    if (i == 1) plotDat_IM<-plotDat.CU_IM
    if (i > 1) plotDat_IM<-rbind(plotDat_IM, plotDat.CU_IM)
    
    # create data frame for IM_surv_prior_cap model
    plotDat.CU_IM_cap<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
                              Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
                              Pred_Rec=Preds_Rec_cap[Preds_Rec_cap$CU_ID == CUID_list[i], "Estimate"] * Scale,
                              Pred_Rec_SE=Preds_Rec_cap[Preds_Rec_cap$CU_ID == CUID_list[i], "Std..Error"] * Scale)
    
    if (i == 1) plotDat_IM_cap<-plotDat.CU_IM_cap
    if (i > 1) plotDat_IM_cap<-rbind(plotDat_IM_cap, plotDat.CU_IM_cap)
    
    
    
} # end of nCU loop to create data frames
  

plotDat_IM$Pred_Rec_lwr <- plotDat_IM$Pred_Rec  - (plotDat_IM$Pred_Rec_SE * 1.96)
plotDat_IM$Pred_Rec_upr <- plotDat_IM$Pred_Rec + (plotDat_IM$Pred_Rec_SE * 1.96)

plotDat_IM_cap$Pred_Rec_lwr <- plotDat_IM_cap$Pred_Rec  - (plotDat_IM_cap$Pred_Rec_SE * 1.96)
plotDat_IM_cap$Pred_Rec_upr <- plotDat_IM_cap$Pred_Rec + (plotDat_IM_cap$Pred_Rec_SE * 1.96)
 

# Define ggplot function to be applied to each CU:
makeSRplots<-function(i,plotDat,plotDat_cap, SRDat) {
  
  Ylab = "Recruits"
  Xlab = "Spawners"
  
  if (useFrenchCaptions==TRUE) {
    Ylab = "Recrues"
    Xlab = "Géniteurs"
  }
  
  dat<-plotDat %>% filter(CU_ID == i)
  dat_cap<-plotDat_cap %>% filter(CU_ID == i)
  SR_dat <- SRDat %>% filter(CU_ID == i)
  
  xMax<-max(c(SR_dat$Recruits,SR_dat$Spawners))*1.2
  yMax<-max(c(SR_dat$Recruits,SR_dat$Spawners))*1.4
  
  if (i == 1) yMax<-max(c(SR_dat$Recruits,SR_dat$Spawners))*2

  CUname<-unique(SRDat$CU_Name[SRDat$CU_ID==i])
  
  # Specify CU names for french translation
  if (useFrenchCaptions == TRUE) {
    if(i == 0) CUname<-"Moyen Fraser"
    if(i == 1) CUname<-"Canyon du Fraser"
    if(i == 2) CUname<-"Thompson inférieure"
    if(i == 3) CUname<-"Thompson Nord"
    if(i == 4) CUname<-"Thompson Sud"
  }
  
  # Create plot
  p <- ggplot(data=dat, mapping=aes(x=Pred_Spwn,y=Pred_Rec)) +
    # add hier surv model
    geom_ribbon(aes(ymin = Pred_Rec_lwr, ymax = Pred_Rec_upr, x=Pred_Spwn), fill = "black", alpha=0.1) +
    geom_line(mapping=aes(x=Pred_Spwn, y=Pred_Rec), col="black", size=1) +
    geom_line(mapping=aes(x=Pred_Spwn, y=Pred_Rec_upr), col="black", lty=2) +
    geom_line(mapping=aes(x=Pred_Spwn, y=Pred_Rec_lwr), col="black", lty=2) +
    # add hier surv model with prior cap
    geom_ribbon(data=dat_cap, aes(ymin = Pred_Rec_lwr, ymax = Pred_Rec_upr, x=Pred_Spwn), fill = "blue", alpha=0.1) +
    geom_line(data=dat_cap, mapping=aes(x=Pred_Spwn, y=Pred_Rec), col="blue", size=1) +
    geom_line(data=dat_cap, mapping=aes(x=Pred_Spwn, y=Pred_Rec_upr), col="blue", lty=2) +
    geom_line(data=dat_cap, mapping=aes(x=Pred_Spwn, y=Pred_Rec_lwr), col="blue", lty=2) +
    # add data
    geom_point(data=SR_dat, mapping=aes(x=Spawners,y=Recruits) ) +
    # add replacement line
    geom_line(mapping=aes(x=Pred_Spwn, y=Pred_Spwn), col="red") +
    # add title, labels, theme
    ggtitle(CUname) +
    xlab(Xlab) + ylab(Ylab) +
    xlim(0, xMax) + ylim(0, yMax) +
    theme_classic()

}

# Create multi-panel plots of SR fits for IM and IM_cap model =====================

pngName<-ifelse(useFrenchCaptions==FALSE, "coho-compare-SRFits-IM.png", "coho-compare-SRFits-IM-FN.png")

ps<-lapply(CUID_list, makeSRplots, plotDat=plotDat_IM, plotDat_cap=plotDat_IM_cap, SRDat=SRDat)
png(paste(figDir, pngName, sep="/"))
do.call(grid.arrange,  ps)
dev.off()
 



# ========================================================================================================
# Create plots to compare parameter estimates
# ================================================================================================
# 
# 
# # Define ggplot function to be applied to each CU:
# makeParEstPlots<-function(i, est_hier, est_hierCap, est_IM, est_IMCap, parName) {
#   
#   model<-rep(NA,4)
#   parEst<-rep(NA,4)
#   upr<-rep(NA,4)
#   lwr<-rep(NA,4)
#   
#   for (mod in 1:4) {
#     if (mod == 1) dat<-est_hier %>% filter(CU_ID == i, Param == parName)
#     if (mod == 2) dat<-est_IM %>% filter(CU_ID == i, Param == parName)
#     if (mod == 3) dat<-est_hierCap %>% filter(CU_ID == i, Param == parName)
#     if (mod == 4) dat<-est_IMCap %>% filter(CU_ID == i, Param == parName)
#     
#     parEst[mod]<-dat$Estimate
#     upr[mod]<-dat$Estimate + dat$Std..Error * 1.96
#     lwr[mod]<-dat$Estimate - dat$Std..Error * 1.96
#   }
#   
#   model[1]<-"HM"
#   model[2]<-"IM"
#   model[3]<-"HM.HiSRep"
#   model[4]<-"IM.HiSRep"
#   
#   plot.dat<-data.frame(model,parEst,lwr,upr)
#   
#   # assign plot order:
#   plot.dat$model<-factor(plot.dat$model, levels = c("IM", "HM", "IM.HiSRep", "HM.HiSRep"))
#   
#   p<- ggplot(data=plot.dat, mapping=aes(x=model, y=parEst)) +
#     geom_errorbar(aes(x=model,ymax=upr,ymin=lwr), width=0,colour="black") +
#     geom_point(mapping=aes(x=model, y=parEst), col="black", size=2) +
#     xlab("") + ylab(parName) +
#     ylim(min(plot.dat$lwr), max(plot.dat$upr)) +v
#     theme_classic() +
#     ggtitle(CU_list[CUID_list == i])
#   
# }
# 
# est_hier<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_Hier_Ricker_Surv.csv", sep=""))
# est_hierCap<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_Hier_Ricker_Surv_priorCap.csv", sep=""))
# est_IM<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_IM_Ricker_Surv.csv", sep=""))
# est_IMCap<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_IM_Ricker_Surv_priorCap.csv", sep=""))
# 
# # Set-up plot to compare LRP estimates among models
# LRP_ests<-data.frame(rbind(est_hier[est_hier$Param=="Agg_LRP",],
#                            est_hierCap[est_hierCap$Param=="Agg_LRP",],
#                            est_IM[est_IM$Param=="Agg_LRP",],
#                            est_IMCap[est_IMCap$Param=="Agg_LRP",]))
# 
# LRP_ests$ModName[LRP_ests$Mod == "SR_HierRicker_Surv"]<-"HM"
# LRP_ests$ModName[LRP_ests$Mod == "SR_IndivRicker_Surv"]<-"IM"
# LRP_ests$ModName[LRP_ests$Mod == "SR_HierRicker_SurvCap"]<-"HM.HiSRep"
# LRP_ests$ModName[LRP_ests$Mod == "SR_IndivRicker_SurvCap"]<-"IM.HiSRep"
# 
# 
# LRP_ests$Upper<-LRP_ests$Estimate + (1.96*LRP_ests$Std..Error)
# LRP_ests$Lower<-LRP_ests$Estimate - (1.96*LRP_ests$Std..Error)
# 
# # assign plot order:
# LRP_ests$ModName<-factor(LRP_ests$ModName, levels = c("IM", "HM", "IM.HiSRep", "HM.HiSRep"))
# 
# 
# png(paste(cohoDir,"/Figures/", "compareSRpar_LRP_Ests.png", sep=""))
# p.lrp<- ggplot(data=LRP_ests, mapping=aes(x=ModName, y=Estimate)) +
#   geom_errorbar(aes(x=ModName,ymax=Upper,ymin=Lower), width=0,colour="black") +
#   geom_point(mapping=aes(x=ModName, y=Estimate), col="black", size=2) +
#   xlab("") + ylab("Aggregate LRP") +
#   ylim(min(LRP_ests$Lower), max(LRP_ests$Upper)) +
#   theme_classic()
# p.lrp
# dev.off()
# 
# 
# # Set-up Sgen plot
# ps<-lapply(CUID_list, makeParEstPlots, est_hier=est_hier, est_hierCap=est_hierCap, est_IM=est_IM, est_IMCap=est_IMCap, parName = "Sgen")
# # Add LRP to Sgen plot
# #ps[[6]]<- p.lrp
# 
# png(paste(cohoDir,"/Figures/", "compareSRpar_Sgen_Ests.png", sep=""))
# do.call(grid.arrange,  ps)
# dev.off()
# 
# ps<-lapply(CUID_list, makeParEstPlots, est_hier=est_hier, est_hierCap=est_hierCap, est_IM=est_IM, est_IMCap=est_IMCap, parName = "A")
# png(paste(cohoDir,"/Figures/", "compareSRpar_A_Ests.png", sep=""))
# do.call(grid.arrange,  ps)
# dev.off()
# 
# ps<-lapply(CUID_list, makeParEstPlots, est_hier=est_hier, est_hierCap=est_hierCap, est_IM=est_IM, est_IMCap=est_IMCap, parName = "logA")
# png(paste(cohoDir,"/Figures/", "compareSRpar_logA_Ests.png", sep=""))
# do.call(grid.arrange,  ps)
# dev.off()
# 
# 
# 
