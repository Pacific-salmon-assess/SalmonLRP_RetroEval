# Code used to run sensitivity analyses of prior distribution assumptions for Integrated LRP model
# Interior Fraser Coho case study
# --> where, the stock recruitment model is Hierarchical Ricker with survival co-variate (model 1 from Arbeider et al. 2019)


library(TMB)
library(dplyr)

# read in data
priorDir<-getwd()
setwd('..')
cohoDir<-getwd()
setwd('..')
tmpDir<-getwd()
codeDir<-paste(tmpDir,"/Code",sep="")
setwd(cohoDir)
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")
setwd(codeDir)

# Compile and load LRP estimation model written in TMB
compile("TMB_Files/SR_HierRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv"))


sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
}
sourceAll()


# Change header names to match generic data headers
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"


# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)

# Prep dataframes to work with functions
CoSRDat$yr_num <- CoSRDat$BroodYear - min(CoSRDat$BroodYear)
CoSRDat$CU_ID <- group_indices(CoSRDat, CU_ID) - 1
CoEscpDat$yr_num <- group_indices(CoEscpDat, BroodYear) - 1
CoEscpDat<- CoEscpDat %>% right_join(unique(CoSRDat[,c("CU_ID", "CU_Name")]))


### Plot functions

# plot LRPs
plot_LRPs <- function(LRPs, SDs, legend_text, name) {
  pdf(paste(priorDir, "/Figures/LRP_Ests",name, ".pdf", sep=""), height = 6, width = 10)
  par(mfrow = c(1,2), mar=c(3,3,1,1), oma=c(2,2,2,2), mgp=c(2, 0.7,0))
  xlims <- c(0, length(LRPs) +1) 
  ylims <- c(min(LRPs - 1.96*SDs), max(LRPs + 1.96*SDs))
  plot(1:length(LRPs), LRPs, pch=19, xlab = "Prior Scenarios", 
       ylab = "Agg LRP", ylim = ylims, xlim = xlims)
  segments(x0 = 1:length(LRPs), x1=1:length(LRPs), y0=LRPs-1.96*SDs, y1 = LRPs+1.96*SDs)
  
  plot.new()
  legend("left", pch = as.character(unique(Ests$Scenario)), legend = legend_text, bty="n", cex=0.8, xpd=NA)
  
  dev.off()
}



PlotEsts <- function(Ests, param = "alpha", legend_text, name = "test") {
  nScens <- length(unique(Ests$Scenario))
  pdf(paste(priorDir, "/Figures/", param, "_Ests_",name,".pdf", sep=""))
  par(mfrow = c(3,2), mar=c(3,3,1,1), oma=c(2,2,2,2), mgp=c(2, 0.7,0))
  for(ss in 1:length(Stks)){
    Sdat <- Ests %>% filter(Stk == Stks[[ss]])
    Param_Ests <- Sdat %>% pull(param)
    SDs <- Sdat %>% pull(paste(param, "sd", sep="_"))
    ylims <- c(min(Param_Ests - 1.96*SDs), max(Param_Ests + 1.96*SDs))
    # Temporarily add hardwired ylims to standardize y-axes among stocks:
    if (param == "alpha") ylims <- c(1.8, 3.4)
    
    xlims <- c(0, nScens +1) 
    plot(1:nScens, Param_Ests, pch=19, xlab = "Prior Scenarios", 
         ylab = paste(param, "Estimate"), ylim = ylims, xlim = xlims)
    segments(x0 = 1:nScens, x1=1:nScens, y0=Param_Ests-1.96*SDs, y1 = Param_Ests+1.96*SDs)
    
    # Temporarily add Sgen ests from previous papers:
    if (param == "Sgen") {
      Korman.Sgens.means<-c(1069, 618, 1561, 2670, 1997)
      WSP.Sgens.medians<-c(1585, 741, 1405, 2546, 2377)
      WSP.Sgens.means<-c(1650, 772, 1489, 2603, 2511)
      
      abline(h=Korman.Sgens.means[ss], col="red", lty=2, lwd=2)
      abline(h=WSP.Sgens.means[ss], col="blue", lty=2, lwd=2)
      
    }
    
    if (param == "alpha") {
      abline(h=exp(1), col="red", lty=2, lwd=2)
    }
    
    mtext(side = 3, text = Stks[ss])
  }
  plot.new()
  legend("left", pch = as.character(1:nScens), legend = legend_text, bty="n", cex=0.8, xpd=NA)
  
  dev.off()
}



# Arbeider et al. pars
TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                   logMuA_sig = 2, Tau_dist = 0.01, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 0.5)


# Final parameter set selected by Brooke based on previous sensitivity analysis:
# TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 2.5, 
#                    logMuA_sig = 2, Tau_dist = 0.1, Tau_A_dist = 0.1, 
#                    gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


# Call base function to fit integrated LRP model, where "Run_Ricker_LRP" comes from "LRPFunctions.r"
xx <- Run_Ricker_LRP(SRDat = CoSRDat, EscDat = CoEscpDat, BM_Mod = "SR_HierRicker_Surv", 
               Bern_Logistic = F, useGenMean = F, genYrs=3, p=0.95, TMB_Inputs = TMB_Inputs)


# Run sensitivity analysis on prior choice =============================================
# Specify parameter values to be tested
gamma_sigs <- c(100,10,1)
IG_params <- c( 0.01, 0.1, 1)
logMuA_mean <- c( 1, 2.5, 3)
logMuA_sig  <- c(1,2,3)


outputs <- list()
Priors <- data.frame(Scenario = numeric(), Gamma_Sig = numeric(), Tau_p = numeric(), 
                     Tau_a = numeric(), logMuA= numeric(), logMuA_sig = numeric())
Ests <- data.frame(Scenario = numeric(), Stk = character(), alpha = numeric(), alpha_sd = numeric(), Sgen = numeric(), Sgen_sd = numeric())
Stks <- unique(CoSRDat$CU_Name)
N_Stks <- length(Stks)
i <- 1

# Loop over parameter values to test
  for(tau_a in 1:length(IG_params)){
    for(tau_p in 1:length(IG_params)){
    for(gs in 1:length(gamma_sigs)){
      for(lmm in 1:length(logMuA_mean)){
        for(lms in 1:length(logMuA_sig)){
      outputs[[i]] <- list()
      outputs[[i]]$params <- data.frame(Scenario = i, Gamma_Sig = gamma_sigs[gs], Tau_p = IG_params[tau_p], 
                                        Tau_a = IG_params[tau_a], logMuA_mean = logMuA_mean[lmm], logMuA_sig = logMuA_sig[lms])
      Priors <- rbind(Priors, outputs[[i]]$params)
      
      # Set-up TMB inputs for sensitivity scenario:      
      TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = logMuA_mean[lmm], logMuA_sig = logMuA_sig[lms], 
                          Tau_dist = IG_params[tau_p], Tau_A_dist = IG_params[tau_a], 
                         gamma_mean = 0, gamma_sig = gamma_sigs[gs], S_dep = 1000, Sgen_sig = 1)
    
      # Call function to fit integrated LRP model for sensitivity scenario:
      outputs[[i]]$mod <- Run_Ricker_LRP(SRDat = CoSRDat, EscDat = CoEscpDat, BM_Mod = "SR_HierRicker_Surv", 
                                         Bern_Logistic = F, useGenMean = F, genYrs=3, p=0.95, TMB_Inputs = TMB_Inputs)
      
      
      alphas <-  outputs[[i]]$mod[[1]] %>% filter(Param == "logA") 
      betas <-  outputs[[i]]$mod[[1]] %>% filter(Param == "B") 
      sigmas <- outputs[[i]]$mod[[1]] %>% filter(Param == "sigma")
      Sgens <- outputs[[i]]$mod[[1]] %>% filter(Param == "Sgen") 
      new.rows <- data.frame(Scenario = rep(i, N_Stks), Stk = Stks, 
                             alpha = alphas$Estimate, alpha_sd = alphas$"Std..Error",
                             beta = betas$Estimate, beta_sd = betas$"Std..Error",
                             sigma = sigmas$Estimate, sigma_sd = sigmas$"Std..Error",
                             Sgen = Sgens$Estimate, Sgen_sd = Sgens$"Std..Error")
      Ests <- rbind(Ests, new.rows)
      i <- i+1
        }
      }
    }
  }
}

    
# Specify different sets of scenarios to be plotted:

#    1) Effect of Taus ===================================================================================

Scens <- Priors %>% filter(logMuA_mean == 1, logMuA_sig == 2, Gamma_Sig == 100) %>% pull(Scenario)
legend_text_tau = NULL

for(i in 1:length(Scens)){
#   legend_text_tau[i] <- paste("Tau~gamma(", outputs[[Scens[i]]]$params$Tau_p, ",", outputs[[Scens[i]]]$params$Tau_p, "), ",
 #                          "Tau_a~gamma(", outputs[[Scens[i]]]$params$Tau_a, ",",outputs[[Scens[i]]]$params$Tau_a , "), ", sep="")
  # switched legend order
  legend_text_tau[i] <- paste("Tau_a~gamma(", outputs[[Scens[i]]]$params$Tau_a, ",", outputs[[Scens[i]]]$params$Tau_a, "), ",
                             "Tau~gamma(", outputs[[Scens[i]]]$params$Tau_p, ",",outputs[[Scens[i]]]$params$Tau_p , "), ", sep="")
}
# Plot LRPs ==============
LRPs_tau <- NULL
SDs_tau <- NULL
for(i in 1:length(Scens)){
  LRPs_tau[i] <- outputs[[Scens[i]]]$mod[[1]] %>% filter(Param == "Agg_LRP") %>% pull(Estimate)
  SDs_tau[i] <- outputs[[Scens[i]]]$mod[[1]] %>% filter(Param == "Agg_LRP") %>% pull(Std..Error)
}

 
plot_LRPs(LRPs_tau, SDs_tau, legend_text_tau, name = "Comp_tau")


# Plot CU-level alpha and Sgens =========================

Ests_Tau <- Ests %>% filter(Scenario %in% Scens)


PlotEsts(Ests = Ests_Tau, param = "alpha", legend_text_tau, name = "Comp_Tau")
PlotEsts(Ests = Ests_Tau, param = "Sgen", legend_text_tau, name = "Comp_Tau")
PlotEsts(Ests = Ests_Tau, param = "beta", legend_text_tau, name = "Comp_Tau")
PlotEsts(Ests = Ests_Tau, param = "sigma", legend_text_tau, name = "Comp_Tau")

# Save sensitivity analysis results for tau in a csv file

out.df<-data.frame(legend_text_tau,Ests_Tau[,-1])

names(out.df)[1:2]<-c("Scenario", "CU")

write.csv(out.df,file=paste(priorDir, "/Figures/Tau_SensAnal_ParamEsts.csv", sep=""))

 
#    2) effect of gammas =============================================================================

Scens <- Priors %>% filter(logMuA_mean == 1, logMuA_sig == 2, Tau_p == 0.01, Tau_a == 0.01) %>% pull(Scenario)
legend_text_gam = NULL
for(i in 1:length(Scens)){
  legend_text_gam[i]  <- paste("gamma~norm(0," , outputs[[Scens[i]]]$params$Gamma_Sig, ")", sep="") 
}
# plot LRPs
LRPs_gam <- NULL
SDs_gam <- NULL
for(i in 1:length(Scens)){
  LRPs_gam[i] <- outputs[[Scens[i]]]$mod[[1]] %>% filter(Param == "Agg_LRP") %>% pull(Estimate)
  SDs_gam[i] <- outputs[[Scens[i]]]$mod[[1]] %>% filter(Param == "Agg_LRP") %>% pull(Std..Error)
}

plot_LRPs(LRPs_gam, SDs_gam, legend_text_gam, name = "Comp_gamma_sig")

Ests_gam = Ests %>% filter(Scenario %in% Scens)

PlotEsts(Ests = Ests_gam, param = "alpha", legend_text_gam, name = "Comp_gamma_sig")
PlotEsts(Ests = Ests_gam, param = "Sgen", legend_text_gam, name = "Comp_gamma_sig")



#    3) effect of alphas

Scens <- Priors %>% filter( Tau_p == 0.01, Tau_a == 0.01, Gamma_Sig == 1) %>% pull(Scenario)
legend_text_alpha = NULL
for(i in 1:length(Scens)){
  legend_text_alpha[i] <- paste("Mu_Alpha~norm(", outputs[[Scens[i]]]$params$logMuA_mean,",",
                              outputs[[Scens[i]]]$params$logMuA_sig, ")", sep="") 
}

# plot LRPs
LRPs_alpha <- NULL
SDs_alpha <- NULL
for(i in 1:length(Scens)){
  LRPs_alpha[i] <- outputs[[Scens[i]]]$mod[[1]] %>% filter(Param == "Agg_LRP") %>% pull(Estimate)
  SDs_alpha[i] <- outputs[[Scens[i]]]$mod[[1]] %>% filter(Param == "Agg_LRP") %>% pull(Std..Error)
}

plot_LRPs(LRPs_alpha, SDs_alpha, legend_text_alpha, name = "Comp_alpha3")

Ests_alpha = Ests %>% filter(Scenario %in% Scens)

PlotEsts(Ests = Ests_alpha, param = "alpha", legend_text_alpha, name = "Comp_alpha3")
PlotEsts(Ests = Ests_alpha, param = "Sgen", legend_text_alpha, name = "Comp_alpha3")

setwd(priorDir)
