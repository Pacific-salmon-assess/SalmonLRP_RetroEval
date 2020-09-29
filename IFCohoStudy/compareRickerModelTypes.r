
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

setwd(codeDir)

compile("TMB_Files/ThreshAbund_Subpop1000.cpp")
dyn.load(dynlib("TMB_Files/ThreshAbund_Subpop1000"))

compile("TMB_Files/SR_HierRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv"))

compile("TMB_Files/SR_IndivRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv"))

compile("TMB_Files/SR_HierRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap"))

compile("TMB_Files/SR_IndivRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap"))


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

CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")

  CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
  tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
  CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)

setwd(codeDir)
  

# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)

# Prep data frame
CoSRDat$yr_num <- CoSRDat$BroodYear - min(CoSRDat$BroodYear)
CoSRDat$CU_ID <- group_indices(CoSRDat, CU_ID) - 1

CoEscpDat$yr_num <- group_indices(CoEscpDat, BroodYear) - 1
CoEscpDat<- CoEscpDat %>% right_join(unique(CoSRDat[,c("CU_ID", "CU_Name")]))


SRDat<-CoSRDat
EscDat <- CoEscpDat
Bern_Logistic <- F
useGenMean <- F
genYrs <- 3
p<-0.80

Scale<-1000


# Create vector of spawner abundances to plot predicted SR relationships from 
Pred_Spwn <- rep(seq(0,30000/Scale,length=301), 5)
Pred_Spwn_CU <- c(rep(0,301), rep(1,301), rep(2,301), rep(3,301), rep(4,301))


# *************************************************************************************
# Fit Hier_Ricker_Surv model
# ********************************************************************************

# TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                    logMuA_sig = 2, Tau_dist = 0.01, Tau_A_dist = 0.01, 
#                    gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 0.5)

# Exact inputs from Arbeider et al:
TMB_Inputs <- list(Scale = Scale, logA_Start = 1, logMuA_mean = 1, 
                   logMuA_sig = sqrt(2), Tau_dist = 0.01, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)



Mod <- "SR_HierRicker_Surv"

# Set-up for call to TMB
#Scale <- TMB_Inputs$Scale

data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)

# also give model year for which will fit logistic model
# only give S_Compare values when all obs for all CUs
# give stocks associated with S_Compare
# which years do we have obs for all CUs?
Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)

Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1

if(useGenMean == T){
  GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
    filter(is.na(Gen_Mean) == F) %>% ungroup() %>% mutate(yr_num = yr - min(yr))
  
  data$LM_S <- GenMean_DF$Gen_Mean / Scale
  data$LM_yr <- GenMean_DF$yr_num
  data$LM_stk <- GenMean_DF$CU_ID
  Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == F) 
  data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
} else {
  data$LM_S <- Logistic_Dat$Escp / Scale
  data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
  data$LM_yr <- Logistic_Dat$yr_num
  data$LM_stk <- Logistic_Dat$CU_ID
}
# range of agg abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
# range of spawner abundance to predict recruitment from
data$Pred_Spwn <- Pred_Spwn 
data$stk_predS <- Pred_Spwn_CU 

data$p <- p

# set variance to be used for likelihood for estimating Sgen
data$Sgen_sig <- TMB_Inputs$Sgen_sig

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
param$B_0 <- 2
param$B_1 <- 0.1
param$logMuA <- TMB_Inputs$logA_Start
param$logSigmaA <- 1
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)

muSurv <- SRDat %>% group_by(CU_ID) %>% 
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean 
data$gamma_sig <- TMB_Inputs$gamma_sig
data$logMuA_mean <- TMB_Inputs$logMuA_mean 
data$logMuA_sig <- TMB_Inputs$logMuA_sig
data$Tau_A_dist <- TMB_Inputs$Tau_A_dist 


# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)


# Create upper and lower bounds vectors that are same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower) 

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
# summary(sdreport(obj))

# Phase 3 fit logistic model ======================
# Hold other estimates constant

map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
# logMuA and logSigmaA disappear after mask, so need to save here
HyperParams <- data.frame(summary(sdreport(obj))) %>% 
  mutate(Param = row.names(summary(sdreport(obj)))) %>% 
  filter(Param %in% c("logMuA", "logSigmaA"))
obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)


upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower)

opt <- tryCatch(
  {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper, lower=lower)},
  error=function(cond) {
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    return(NA)
  }
)   

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests <- full_join(All_Ests, HyperParams)

# put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests$Mod <- Mod
All_Ests$CU_ID[!(All_Ests$Param %in% c("logMuA", "logSigmaA", "B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests$CU_ID[All_Ests$Param=="Rec_Preds"]<-data$stk_predS
All_Ests <- left_join(All_Ests, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests$Estimate[All_Ests$Param == "logSigma"] <- exp(All_Ests$Estimate[All_Ests$Param == "logSigma"] )
All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
All_Ests$Estimate[All_Ests$Param == "logB"] <- exp(All_Ests$Estimate[All_Ests$Param == "logB"] )
All_Ests$Param[All_Ests$Param == "logB"] <- "B"
All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
All_Ests[All_Ests$Param %in% c("Sgen", "Smsy","Agg_LRP", "SRep"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "Smsy", "Agg_LRP", "SRep")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds <- All_Ests %>% filter(Param == "Logit_Preds")
Preds_Rec <- All_Ests %>% filter(Param == "Rec_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds", "Rec_Preds"))) 

write.csv(All_Ests,paste(cohoDir,"/DataOut/ModelFits/AllEsts_Hier_Ricker_Surv.csv", sep=""))

# ========================================================================================

# *************************************************************************************
# Fit Hier_Ricker_Surv_PriorCap model
# ********************************************************************************

# Extract capacity estimates from Hier_Ricker_Surv model fit to use as priors on carrying capacity
CU_Names <- unique(SRDat[, "CU_Name"])
cap<-rep(NA,N_Stocks)

for (i in 1:N_Stocks) {
 logA<-All_Ests[All_Ests$Param =="logA" & All_Ests$CU_Name == CU_Names[i],"Estimate"]
 A<-All_Ests[All_Ests$Param =="A" & All_Ests$CU_Name == CU_Names[i],"Estimate"]
 B<-All_Ests[All_Ests$Param =="B" & All_Ests$CU_Name == CU_Names[i],"Estimate"] * Scale
 cap[i] <- log(logA)/B
}

cap_priorMean<- cap*1.5

print("Prior means for capacity in Hier model, by CU:")
print(cap_priorMean)

# Compile TMB model

# TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                    logMuA_sig = 2, Tau_dist = 0.01, Tau_A_dist = 0.01, 
#                    gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 0.5,
#                    cap_mean=cap_priorMean, cap_sig=0.5)

# Exact inputs from Arbeider et al:
TMB_Inputs <- list(Scale = Scale, logA_Start = 1, logMuA_mean = 1, 
                   logMuA_sig = sqrt(2), Tau_dist = 0.01, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                  # gamma_mean = 0.375, gamma_sig = 0.02, S_dep = 1000, Sgen_sig = 1, # forcing informative prior on gamma
                   cap_mean=cap_priorMean, cap_sig=sqrt(2))


Mod <- "SR_HierRicker_SurvCap"

# Set-up for call to TMB
#Scale <- TMB_Inputs$Scale

data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)

# also give model year for which will fit logistic model
# only give S_Compare values when all obs for all CUs
# give stocks associated with S_Compare
# which years do we have obs for all CUs?
Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)

Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1

if(useGenMean == T){
  GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
    filter(is.na(Gen_Mean) == F) %>% ungroup() %>% mutate(yr_num = yr - min(yr))
  
  data$LM_S <- GenMean_DF$Gen_Mean / Scale
  data$LM_yr <- GenMean_DF$yr_num
  data$LM_stk <- GenMean_DF$CU_ID
  Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == F) 
  data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
} else {
  data$LM_S <- Logistic_Dat$Escp / Scale
  data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
  data$LM_yr <- Logistic_Dat$yr_num
  data$LM_stk <- Logistic_Dat$CU_ID
}
# range of aggregate abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
# range of spawner abundance to predict recruitment from
data$Pred_Spwn <- Pred_Spwn 
data$stk_predS <- Pred_Spwn_CU 

data$p <- p

# set variance to be used for likelihood for estimating Sgen
data$Sgen_sig <- TMB_Inputs$Sgen_sig

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$cap <- TMB_Inputs$cap_mean
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
param$B_0 <- 2
param$B_1 <- 0.1
param$logMuA <- TMB_Inputs$logA_Start
param$logSigmaA <- 1
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)

muSurv <- SRDat %>% group_by(CU_ID) %>% 
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean 
data$gamma_sig <- TMB_Inputs$gamma_sig
data$logMuA_mean <- TMB_Inputs$logMuA_mean 
data$logMuA_sig <- TMB_Inputs$logMuA_sig
data$Tau_A_dist <- TMB_Inputs$Tau_A_dist 
data$cap_mean<-TMB_Inputs$cap_mean
data$cap_sig<-TMB_Inputs$cap_sig


# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)


# Create upper and lower bounds vectors that are same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower) 

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower)
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
# summary(sdreport(obj))

# Phase 3 fit logistic model ======================
# Hold other estimates constant

map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
# logMuA and logSigmaA disappear after mask, so need to save here
HyperParams <- data.frame(summary(sdreport(obj))) %>% 
  mutate(Param = row.names(summary(sdreport(obj)))) %>% 
  filter(Param %in% c("logMuA", "logSigmaA"))
obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)

# Create upper and lower bounds vectors that are same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower) 

opt <- tryCatch(
  {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper, lower=lower)},
  error=function(cond) {
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    return(NA)
  }
)   


# Create Table of outputs
All_Ests_cap <- data.frame(summary(sdreport(obj)))
All_Ests_cap$Param <- row.names(All_Ests_cap)
All_Ests_cap <- full_join(All_Ests_cap, HyperParams)

# put together readable data frame of values
All_Ests_cap$Param <- sapply(All_Ests_cap$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests_cap$Mod <- Mod
All_Ests_cap$CU_ID[!(All_Ests_cap$Param %in% c("logMuA", "logSigmaA", "B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests_cap$CU_ID[All_Ests_cap$Param=="Rec_Preds"]<-data$stk_predS
All_Ests_cap <- left_join(All_Ests_cap, unique(SRDat[, c("CU_ID", "CU_Name")]))



# don't want logged or scaled param values, so need to convert
All_Ests_cap$Estimate[All_Ests_cap$Param == "logSigma"] <- exp(All_Ests_cap$Estimate[All_Ests_cap$Param == "logSigma"] )
All_Ests_cap$Param[All_Ests_cap$Param == "logSigma"] <- "sigma"
All_Ests_cap$Param[All_Ests_cap$Param == "SMSY"] <- "Smsy"
All_Ests_cap[All_Ests_cap$Param == "B",] <- All_Ests_cap %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests_cap[All_Ests_cap$Param %in% c("Sgen", "Smsy", "SRep", "cap", "Agg_LRP"), ] <-  All_Ests_cap %>% filter(Param %in% c("Sgen", "Smsy", "SRep","cap", "Agg_LRP")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds <- All_Ests_cap %>% filter(Param == "Logit_Preds")
Preds_Rec_cap <- All_Ests_cap %>% filter(Param == "Rec_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds", "Rec_Preds"))) 

write.csv(All_Ests_cap,paste(cohoDir,"/DataOut/ModelFits/AllEsts_Hier_Ricker_Surv_priorCap.csv", sep=""))


# *************************************************************************************
# Fit Individual Ricker_Surv model
# ********************************************************************************

TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.01,
                      gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 1)


Mod <- "SR_IndivRicker_Surv"


data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)

# also give model year for which will fit logistic model
# only give S_Compare values when all obs for all CUs
# give stocks associated with S_Compare
# which years do we have obs for all CUs?
Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)

Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1

if(useGenMean == T){
  GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
    filter(is.na(Gen_Mean) == F) %>% ungroup() %>% mutate(yr_num = yr - min(yr))
  
  data$LM_S <- GenMean_DF$Gen_Mean / Scale
  data$LM_yr <- GenMean_DF$yr_num
  data$LM_stk <- GenMean_DF$CU_ID
  Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == F) 
  data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
} else {
  data$LM_S <- Logistic_Dat$Escp / Scale
  data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
  data$LM_yr <- Logistic_Dat$yr_num
  data$LM_stk <- Logistic_Dat$CU_ID
}
# range of agg abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
# range of spawner abundance to predict recruitment from
data$Pred_Spwn <- Pred_Spwn 
data$stk_predS <- Pred_Spwn_CU 

data$p <- p

# set variance to be used for likelihood for estimating Sgen
data$Sgen_sig <- TMB_Inputs$Sgen_sig

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
param$B_0 <- 2
param$B_1 <- 0.1
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)

muSurv <- SRDat %>% group_by(CU_ID) %>% 
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean 
data$gamma_sig <- TMB_Inputs$gamma_sig

# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, map=map2)


# Create upper and lower bounds vectors that are same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower) 

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
# summary(sdreport(obj))

# Phase 3 fit logistic model ======================
# Hold other estimates constant

obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE)

# Create upper & lower bounds vector that is same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower)

opt <- tryCatch(
  {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper, lower=lower)},
  error=function(cond) {
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    return(NA)
  }
)   

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests$Mod <- Mod
All_Ests$CU_ID[!(All_Ests$Param %in% c("B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests$CU_ID[All_Ests$Param=="Rec_Preds"]<-data$stk_predS
All_Ests <- left_join(All_Ests, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests$Estimate[All_Ests$Param == "logSigma"] <- exp(All_Ests$Estimate[All_Ests$Param == "logSigma"] )
All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
All_Ests$Estimate[All_Ests$Param == "logB"] <- exp(All_Ests$Estimate[All_Ests$Param == "logB"] )
All_Ests$Param[All_Ests$Param == "logB"] <- "B"
All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
All_Ests[All_Ests$Param %in% c("Sgen", "Smsy","Agg_LRP", "SRep"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "Smsy", "Agg_LRP", "SRep")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds <- All_Ests %>% filter(Param == "Logit_Preds")
Preds_Rec_IM <- All_Ests %>% filter(Param == "Rec_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds", "Rec_Preds"))) 

write.csv(All_Ests,paste(cohoDir,"/DataOut/ModelFits/AllEsts_IM_Ricker_Surv.csv", sep=""))




# *************************************************************************************
# Fit Individual Ricker_Surv model with Prior Cap
# ********************************************************************************

# Extract capacity estimates from Indiv_Ricker_Surv model fit to use as priors on carrying capacity
CU_Names <- unique(SRDat[, "CU_Name"])
cap<-rep(NA,N_Stocks)

for (i in 1:N_Stocks) {
  A<-All_Ests[All_Ests$Param =="A" & All_Ests$CU_Name == CU_Names[i],"Estimate"]
  B<-All_Ests[All_Ests$Param =="B" & All_Ests$CU_Name == CU_Names[i],"Estimate"] * Scale
  cap[i] <- log(A)/B
}

cap_priorMean<- cap*1.5

print("Prior means for capacity in individual models, by CU:")
print(cap_priorMean)


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



# Exact inputs from Arbeider et al:
TMB_Inputs <- list(Scale = Scale, logA_Start = 1, Tau_dist = 0.01, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                   cap_mean=cap_priorMean, cap_sig=sqrt(2))


Mod <- "SR_IndivRicker_SurvCap"


data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)

# also give model year for which will fit logistic model
# only give S_Compare values when all obs for all CUs
# give stocks associated with S_Compare
# which years do we have obs for all CUs?
Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)

Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1

if(useGenMean == T){
  GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
    filter(is.na(Gen_Mean) == F) %>% ungroup() %>% mutate(yr_num = yr - min(yr))
  
  data$LM_S <- GenMean_DF$Gen_Mean / Scale
  data$LM_yr <- GenMean_DF$yr_num
  data$LM_stk <- GenMean_DF$CU_ID
  Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == F) 
  data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
} else {
  data$LM_S <- Logistic_Dat$Escp / Scale
  data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
  data$LM_yr <- Logistic_Dat$yr_num
  data$LM_stk <- Logistic_Dat$CU_ID
}
# range of agg abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
# range of spawner abundance to predict recruitment from
data$Pred_Spwn <- Pred_Spwn 
data$stk_predS <- Pred_Spwn_CU 

data$p <- p

# set variance to be used for likelihood for estimating Sgen
data$Sgen_sig <- TMB_Inputs$Sgen_sig

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$cap <- TMB_Inputs$cap_mean
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
param$B_0 <- 2
param$B_1 <- 0.1
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)

muSurv <- SRDat %>% group_by(CU_ID) %>% 
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean 
data$gamma_sig <- TMB_Inputs$gamma_sig
data$cap_mean<-TMB_Inputs$cap_mean
data$cap_sig<-TMB_Inputs$cap_sig

# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, map=map2)


# Create upper and lower bounds vectors that are same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower) 

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper, lower=lower )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
# summary(sdreport(obj))

# Phase 3 fit logistic model ======================
# Hold other estimates constant

obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE)

# Create upper & lower bounds vector that is same length and order as start vector that will be given to nlminb
upper<-unlist(obj$par)
upper[1:length(upper)]<-Inf
upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
upper<-unname(upper)

lower<-unlist(obj$par)
lower[1:length(lower)]<--Inf
lower[names(lower) =="logSgen"] <- log(0.001)
lower<-unname(lower)

opt <- tryCatch(
  {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper, lower=lower)},
  error=function(cond) {
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    return(NA)
  }
)   

# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests$Mod <- Mod
All_Ests$CU_ID[!(All_Ests$Param %in% c("B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod","Rec_Preds"))] <- rep(0:(N_Stocks-1)) 
All_Ests$CU_ID[All_Ests$Param=="Rec_Preds"]<-data$stk_predS
All_Ests <- left_join(All_Ests, unique(SRDat[, c("CU_ID", "CU_Name")]))

# don't want logged or scaled param values, so need to convert
All_Ests$Estimate[All_Ests$Param == "logSigma"] <- exp(All_Ests$Estimate[All_Ests$Param == "logSigma"] )
All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
All_Ests$Estimate[All_Ests$Param == "logB"] <- exp(All_Ests$Estimate[All_Ests$Param == "logB"] )
All_Ests$Param[All_Ests$Param == "logB"] <- "B"
All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
All_Ests[All_Ests$Param %in% c("Sgen", "Smsy","Agg_LRP", "SRep"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "Smsy", "Agg_LRP", "SRep")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds <- All_Ests %>% filter(Param == "Logit_Preds")
Preds_Rec_IM_cap <- All_Ests %>% filter(Param == "Rec_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds", "Rec_Preds"))) 

write.csv(All_Ests,paste(cohoDir,"/DataOut/ModelFits/AllEsts_IM_Ricker_Surv_priorCap.csv", sep=""))


# ========================================================================================================
# Create plots of fit SR curves
# ================================================================================================

CU_list<-unique(SRDat[, "CU_Name"])
CUID_list<-unique(SRDat[, "CU_ID"])
nCUs<-length(CU_list)


# Create dataframes with predictions to plot
for (i in 1:nCUs) {

  # create data frame for hier_surv model
  plotDat.CU<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
                    Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
                    Pred_Rec=Preds_Rec[Preds_Rec$CU_ID == CUID_list[i], "Estimate"] * Scale,
                    Pred_Rec_SE=Preds_Rec[Preds_Rec$CU_ID == CUID_list[i], "Std..Error"] * Scale)
  
    if (i == 1) plotDat<-plotDat.CU
    if (i > 1) plotDat<-rbind(plotDat, plotDat.CU)
    
    # create data frame for hier_surv_prior_cap model
    plotDat.CU_cap<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
                           Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
                           Pred_Rec=Preds_Rec_cap[Preds_Rec_cap$CU_ID == CUID_list[i], "Estimate"] * Scale,
                           Pred_Rec_SE=Preds_Rec_cap[Preds_Rec_cap$CU_ID == CUID_list[i], "Std..Error"] * Scale)
    
    if (i == 1) plotDat_cap<-plotDat.CU_cap
    if (i > 1) plotDat_cap<-rbind(plotDat_cap, plotDat.CU_cap)
  
    # create data frame for IM_surv model
    plotDat.CU_IM<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
                           Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
                           Pred_Rec=Preds_Rec_IM[Preds_Rec_IM$CU_ID == CUID_list[i], "Estimate"] * Scale,
                           Pred_Rec_SE=Preds_Rec_IM[Preds_Rec_IM$CU_ID == CUID_list[i], "Std..Error"] * Scale)
    
    if (i == 1) plotDat_IM<-plotDat.CU_IM
    if (i > 1) plotDat_IM<-rbind(plotDat_IM, plotDat.CU_IM)
    
    # create data frame for IM_surv_prior_cap model
    plotDat.CU_IM_cap<-data.frame(CU_ID = Pred_Spwn_CU[Pred_Spwn_CU == CUID_list[i]],
                              Pred_Spwn = Pred_Spwn[Pred_Spwn_CU == CUID_list[i]] * Scale,
                              Pred_Rec=Preds_Rec_IM_cap[Preds_Rec_IM_cap$CU_ID == CUID_list[i], "Estimate"] * Scale,
                              Pred_Rec_SE=Preds_Rec_IM_cap[Preds_Rec_IM_cap$CU_ID == CUID_list[i], "Std..Error"] * Scale)
    
    if (i == 1) plotDat_IM_cap<-plotDat.CU_IM_cap
    if (i > 1) plotDat_IM_cap<-rbind(plotDat_IM_cap, plotDat.CU_IM_cap)
    
    
    
  }
  
plotDat$Pred_Rec_lwr <- plotDat$Pred_Rec  - (plotDat$Pred_Rec_SE * 1.96)
plotDat$Pred_Rec_upr <- plotDat$Pred_Rec + (plotDat$Pred_Rec_SE * 1.96)

plotDat_cap$Pred_Rec_lwr <- plotDat_cap$Pred_Rec  - (plotDat_cap$Pred_Rec_SE * 1.96)
plotDat_cap$Pred_Rec_upr <- plotDat_cap$Pred_Rec + (plotDat_cap$Pred_Rec_SE * 1.96)

plotDat_IM$Pred_Rec_lwr <- plotDat_IM$Pred_Rec  - (plotDat_IM$Pred_Rec_SE * 1.96)
plotDat_IM$Pred_Rec_upr <- plotDat_IM$Pred_Rec + (plotDat_IM$Pred_Rec_SE * 1.96)

plotDat_IM_cap$Pred_Rec_lwr <- plotDat_IM_cap$Pred_Rec  - (plotDat_IM_cap$Pred_Rec_SE * 1.96)
plotDat_IM_cap$Pred_Rec_upr <- plotDat_IM_cap$Pred_Rec + (plotDat_IM_cap$Pred_Rec_SE * 1.96)
 

# Define ggplot function to be applied to each CU:
makeSRplots<-function(i,plotDat,plotDat_cap, SRDat) {
  
  Ylab = "Recruits"
  Xlab = "Spawners"
  dat<-plotDat %>% filter(CU_ID == i)
  dat_cap<-plotDat_cap %>% filter(CU_ID == i)
  SR_dat <- SRDat %>% filter(CU_ID == i)
  
  xyMax<-max(c(SR_dat$Recruits,SR_dat$Spawners))*1.2

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
    ggtitle(unique(SRDat$CU_Name[SRDat$CU_ID==i])) +
    xlab(Xlab) + ylab(Ylab) +
    xlim(0, xyMax) + ylim(0, xyMax) +
    theme_classic()
  
}

# Create multi-panel plots of SR fits =====================

ps<-lapply(CUID_list, makeSRplots, plotDat=plotDat, plotDat_cap=plotDat_cap, SRDat=SRDat)
pdf(paste(cohoDir,"/Figures/", "compareSR_fits_hierw.cap.pdf", sep=""))
do.call(grid.arrange,  ps)
dev.off()

ps<-lapply(CUID_list, makeSRplots, plotDat=plotDat_IM, plotDat_cap=plotDat_IM_cap, SRDat=SRDat)
pdf(paste(cohoDir,"/Figures/", "compareSR_fits_IMw.cap.pdf", sep=""))
do.call(grid.arrange,  ps)
dev.off()
 


# ========================================================================================================
# Create plots to compare parameter estimates
# ================================================================================================


# Define ggplot function to be applied to each CU:
makeParEstPlots<-function(i, est_hier, est_hierCap, est_IM, est_IMCap, parName) {
  
  model<-rep(NA,4)
  parEst<-rep(NA,4)
  upr<-rep(NA,4)
  lwr<-rep(NA,4)
  
  for (mod in 1:4) {
    if (mod == 1) dat<-est_hier %>% filter(CU_ID == i, Param == parName)
    if (mod == 2) dat<-est_IM %>% filter(CU_ID == i, Param == parName)
    if (mod == 3) dat<-est_hierCap %>% filter(CU_ID == i, Param == parName)
    if (mod == 4) dat<-est_IMCap %>% filter(CU_ID == i, Param == parName)
    
    parEst[mod]<-dat$Estimate
    upr[mod]<-dat$Estimate + dat$Std..Error * 1.96
    lwr[mod]<-dat$Estimate - dat$Std..Error * 1.96
  }
  
  model[1]<-"HM"
  model[2]<-"IM"
  model[3]<-"HM.HiSRep"
  model[4]<-"IM.HiSRep"
  
  plot.dat<-data.frame(model,parEst,lwr,upr)
  
  # assign plot order:
  plot.dat$model<-factor(plot.dat$model, levels = c("IM", "HM", "IM.HiSRep", "HM.HiSRep"))
  
  p<- ggplot(data=plot.dat, mapping=aes(x=model, y=parEst)) +
    geom_errorbar(aes(x=model,ymax=upr,ymin=lwr), width=0,colour="black") +
    geom_point(mapping=aes(x=model, y=parEst), col="black", size=2) +
    xlab("") + ylab(parName) +
    ylim(min(plot.dat$lwr), max(plot.dat$upr)) +
    theme_classic() +
    ggtitle(CU_list[CUID_list == i])
  
}

est_hier<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_Hier_Ricker_Surv.csv", sep=""))
est_hierCap<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_Hier_Ricker_Surv_priorCap.csv", sep=""))
est_IM<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_IM_Ricker_Surv.csv", sep=""))
est_IMCap<-read.csv(paste(cohoDir,"/DataOut/ModelFits/AllEsts_IM_Ricker_Surv_priorCap.csv", sep=""))


# Set-up plot to compare LRP estimates among models
LRP_ests<-data.frame(rbind(est_hier[est_hier$Param=="Agg_LRP",],
                           est_hierCap[est_hierCap$Param=="Agg_LRP",],
                           est_IM[est_IM$Param=="Agg_LRP",],
                           est_IMCap[est_IMCap$Param=="Agg_LRP",]))

LRP_ests$ModName[LRP_ests$Mod == "SR_HierRicker_Surv"]<-"HM"
LRP_ests$ModName[LRP_ests$Mod == "SR_IndivRicker_Surv"]<-"IM"
LRP_ests$ModName[LRP_ests$Mod == "SR_HierRicker_SurvCap"]<-"HM.HiSRep"
LRP_ests$ModName[LRP_ests$Mod == "SR_IndivRicker_SurvCap"]<-"IM.HiSRep"


LRP_ests$Upper<-LRP_ests$Estimate + (1.96*LRP_ests$Std..Error)
LRP_ests$Lower<-LRP_ests$Estimate - (1.96*LRP_ests$Std..Error)

# assign plot order:
LRP_ests$ModName<-factor(LRP_ests$ModName, levels = c("IM", "HM", "IM.HiSRep", "HM.HiSRep"))


pdf(paste(cohoDir,"/Figures/", "compareSRpar_LRP_Ests.pdf", sep=""))
p.lrp<- ggplot(data=LRP_ests, mapping=aes(x=ModName, y=Estimate)) +
  geom_errorbar(aes(x=ModName,ymax=Upper,ymin=Lower), width=0,colour="black") +
  geom_point(mapping=aes(x=ModName, y=Estimate), col="black", size=2) +
  xlab("") + ylab("Aggregate LRP") +
  ylim(min(LRP_ests$Lower), max(LRP_ests$Upper)) +
  theme_classic()
dev.off()


# Set-up Sgen plot
ps<-lapply(CUID_list, makeParEstPlots, est_hier=est_hier, est_hierCap=est_hierCap, est_IM=est_IM, est_IMCap=est_IMCap, parName = "Sgen")
# Add LRP to Sgen plot
#ps[[6]]<- p.lrp

pdf(paste(cohoDir,"/Figures/", "compareSRpar_Sgen_Ests.pdf", sep=""))
do.call(grid.arrange,  ps)
dev.off()

ps<-lapply(CUID_list, makeParEstPlots, est_hier=est_hier, est_hierCap=est_hierCap, est_IM=est_IM, est_IMCap=est_IMCap, parName = "A")
pdf(paste(cohoDir,"/Figures/", "compareSRpar_A_Ests.pdf", sep=""))
do.call(grid.arrange,  ps)
dev.off()

ps<-lapply(CUID_list, makeParEstPlots, est_hier=est_hier, est_hierCap=est_hierCap, est_IM=est_IM, est_IMCap=est_IMCap, parName = "logA")
pdf(paste(cohoDir,"/Figures/", "compareSRpar_logA_Ests.pdf", sep=""))
do.call(grid.arrange,  ps)
dev.off()
