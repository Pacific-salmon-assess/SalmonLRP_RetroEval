#-------------------------------------------------------------------------------
# Functions to run diagnostics on logistic regression to estimate LRPs for 
# Pacific salmon
# Date last revised: 13 Oct 2021
# Created by: Carrie Holt
#-------------------------------------------------------------------------------

# Function 1: LRdiagnostic(), performs diagnostics using the entire time-series 
# of data

# Function 2: LOO_LRdiagnostic(), performs leave-one-out cross-validation by 
# iteratively re-running LRP estimation of subsets of data, leaving one data 
# point out each time, and estimating classification accuracy on removed data 
# point
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Assumptions of logistic regression
#-------------------------------------------------------------------------------
# 1. The relationship between aggregate abundance and logit(proportions) or log-
# odds is linear.

# 2. The observations are independent of each other (i.e., residuals are not 
# autocorrelated).

# 3. There are no influential outliers. 

# 4. The sample size is sufficiently large. Logistic regression assumes that the 
# sample size of the data set if large enough to draw valid conclusions from the 
# fitted logistic regression model.

# 5. There is no multicollinearity among predictors (NA, only 1 predictor)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function 1: LRdiagnostic(), performs diagnostics using the entire time-series 
# of data (Function 2 is below)

# Steps:
# 1. Box-Tidwell test to assess linearity (Assumption 1). 

# 2. Estimate Pearson residuals and deviance residuals, and plot residuals 
# against fitted values. 

# 3. Plot autocorrelation among residuals. Are residuals autocorrelated? 
# (Assumption 2)

# 4. Are deviance residuals >2, i.e., are they outliers? (Assumption 3)  

# 5. Evaluate if sample size is sufficient. (Assumption 4)

# 6. Evaluate statistical significance of model coefficients using Wald's Test.

# 7. Evaluate goodness-of-fit based on ratio of Deviance to the null model 
# (~likelihood-ratio test).

# 8. Evaluate quasi-Rsquared. What is the ratio of the fit to the null model?

# 9. Evaluate hit ratio. What is the classification accuracy of the LRP based on 
# the logistic regression?

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Arguments:

# SMUlogisticData = Dataframe containing 3 columns: Years, SMU_Esc (annual 
# aggregate escapement to the SMU), and ppn (vector of 0 and 1's indicating if 
# all CUs are above the lower benchmark (1) or not (0); also accepts  proportion 
# of CUs that are above their lower benchmark in each year).

# nCU = number of CUs within the SMU

# All_Ests = Dataframe containing parameters of logistic regression, B_0 and B_1
#   with p values (All_Ests <- data.frame(summary(sdreport(obj), p.value=TRUE)))

# p = the probability that all CUs must be > their lower benchmark when
#   defining the LRP for the Bernoulli logistic regression (or proportion of CUs 
#   for binomial logistic regression)

# Bern_logistic = TRUE/FALSE. Is Bernoulli logistic regression used? (default = 
#   TRUE or 1). If FALSE (or 0), binomial regression used 

# dir = name of directory where plots should be saved followed by /

# plotname = filename for residual plots

#----------------------------------------------------------------------------
## Carrie's inputs for testing
# 
# source("R/WCVILRPs.R")
# zz <- Get.LRP(remove.EnhStocks = TRUE)
# All_Ests <- zz$out$All_Ests
# p <- zz$LRPppn
# Bern_logistic <- as.numeric(FALSE) #For this dummy example on WCVI; usually
#                                    #set to TRUE
# dir <- "DataOut/"
# plotname <- "WCVI_ResidPlots"
# 
# SMUlogisticData <- zz$out$Logistic_Data %>% rename(Years=yr, SMU_Esc =xx ,
#                                                   ppn=yy)
# nCU <- length(names(zz$CU_Status))
# 
# input <- list(SMUlogisticData=SMUlogisticData, nCU=nCU,
#            All_Ests=All_Ests, p=p, Bern_logistic=Bern_logistic, dir="",
#            plotname="test")
# save(input, file="DataIn/Input_LRdiagnostics.rda")
## Not used: nLL<-obj$report()$ans


#----------------------------------------------------------------------------
# Libraries

library(tidyverse)
library(ggplot2)
library(viridis)
library(patchwork)
library(TMB)
library(here)
library(stats)



codeDir <- file.path(dirname(here()), "Code") # get code directory in project parent folder
r_to_source <- dir(codeDir, pattern="\\.R|\\.r", full.names=TRUE) # get r files to source from parent code directory
r_to_source <- r_to_source[-grep("LRdiagnostics.R",r_to_source)]
# function to source all files required
sourceAll <- function(files){
  walk(files, function(x) {source(x, chdir=TRUE)})
}
sourceAll(files= r_to_source) # source all r files from parent directory



# For LOO need to source file from the Watershed-Area-Model repo
# source(here::here("R", "WCVILRPs.R")) # Needed for Step 10 (LOO), specific to 

#----------------------------------------------------------------------------
# Input data
## K.Holt commenting out the following to get code to work with coho case study =========
# caseStudy <- "SCChum"#"WCVIchinook"#"SCChum"# 
# 
# if(caseStudy=="WCVIchinook") {
#   load(here::here("DataIn", "Input_LRdiagnostics.rda"))
#   
#   SMUlogisticData <- input$SMUlogisticData
#   nCU <- input$nCU
#   All_Ests <- input$All_Ests
#   p <- input$p
#   Bern_logistic <- input$Bern_logistic
#   dir <- "DataOut/"#input$dir
#   plotname <- "WCVICK_LRdiag"#input$plotname 
# }
# 
# # ISC Chum
# if(caseStudy=="SCChum") {
#   SCCDir <- file.path(dirname(here()), "SCChumStudy") # get code directory in project parent folder
#   load(paste(SCCDir,"/DataIn/logisticFit_2018.rda" , sep=""))
#   
#   # load(here::here("DataIn", "logisticFit_2018.rda"))#Input_logisticFit_ISC2018v2.rda"))
#   
#   SMUlogisticData <- LRP_Mod$Logistic_Data %>% rename(ppn=yy, SMU_Esc=xx, 
#                                                       Years=yr)
#   nCU <- 4
#   All_Ests <- LRP_Mod$All_Ests
#   p <- 0.5#input$p
#   Bern_logistic <- TRUE#input$Bern_logistic
#   dir <- ""#input$dir
#   plotname <- "SCChumTest"#input$plotname 
#   assumed_added_scalar <- 100 # Need to check this with Luke. I scaled by 10^5 for ISChum, but I think he scaled by 10^3, hence by additional scalar here
#   
# }

# ===================================================================
#-------------------------------------------------------------------------------
# Returns:
# - list of outputs from steps 1-9
# - png of plots of residuals and their autocorrelation stored in directory, dir
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Function:

LRdiagnostics <- function(SMUlogisticData, nCU, All_Ests, p, Bern_logistic, dir, 
                          plotname, caseStudy=NULL, codeDir){
  
  #-----------------------------------------------------------------------------
  # Step 1. Box-Tidwell test to assess linearity between aggregate abundances
  #       and log-odds of all CUs being above their lower benchmarks. 
  #-----------------------------------------------------------------------------
  
  # This code generates inputs for logistic regression and runs the model with 
  # an additional interaction term, aggregate abundances x ln(aggregate 
  # abundance). A significant interaction term indicates that the relationship 
  # is not linear. This code requires that a "Logistic_LRPs_BoxTidwell.cpp" file 
  # with a logistic regression that includes this interaction term has been 
  # compiled and is in the directory "TMB_Inputs". A significant Box-Tidwell 
  # statistic indicates a lack of linearity in the relationship between 
  # aggregate abundances and log-odds (Fox et al. 2016, p. 326-328).

 
  data <- list()
  data$N_Stks <- nCU
  digits <- count.dig(SMUlogisticData$SMU_Esc)
  ScaleSMU <- min(10^(digits -1 ), na.rm=T)
  
  data$LM_Agg_Abund <- SMUlogisticData$SMU_Esc/ScaleSMU
  data$LM_Agg_AbundxLn <- SMUlogisticData$SMU_Esc/ScaleSMU * 
    log(SMUlogisticData$SMU_Esc/ScaleSMU)
  data$N_Above_BM <- SMUlogisticData$ppn * data$N_Stks
  data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund)*1.5, 0.1)
  data$p <- p
  data$Penalty <- as.numeric(FALSE)# Penalty is not relevant with B_2 term.
  
  data$Bern_logistic <- as.numeric(Bern_logistic)
  
  data$B_penalty_mu <- mean(c(200,9962))/ScaleSMU
  data$B_penalty_sig <- 2400/ScaleSMU
  
  param <- list()
  param$B_0 <- -2
  param$B_1 <- 0.1
  param$B_2 <- 0.1

  if (caseStudy == "IF Coho") {
    #codeDir <- file.path(here(), "Code") # K.Holt: if running from IFCohoCohoStudy folder, the here() will not work
    #running from IFCohoCohoStudy
    tmpdir <- getwd()

    codeDir <-paste(strsplit(tmpdir, basename(tmpdir)),"Code", sep="")
    #codeDir <- "C:/github/SalmonLRP_RetroEval/Code"
  }

  
  dyn.load(dynlib(paste(codeDir, "/TMB_Files/Logistic_LRPs_BoxTidwell", sep="")))
  
  obj <- MakeADFun(data, param, DLL="Logistic_LRPs_BoxTidwell", silent=TRUE)
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, 
                control = list(eval.max = 1e5, iter.max = 1e5))
  pl <- obj$env$parList(opt$par) 
  #summary(sdreport(obj), p.value=TRUE)
  All_Ests_BT <- data.frame(summary(sdreport(obj), p.value=TRUE))
  
  All_Ests_BT$Param <- row.names(All_Ests_BT)
  All_Ests_BT$Param <- sapply(All_Ests_BT$Param, function(x) 
    (unlist(strsplit(x, "[.]"))[[1]]))
  
  pBoxTidwell <- All_Ests_BT %>% filter(Param=="B_2") %>% pull(Pr...z.2..)
  pBoxTidwell <- round(pBoxTidwell,2)    
  names(pBoxTidwell) <- c("Box-Tidwell p-value")
  
  pBoxTidwell
  
  #-----------------------------------------------------------------------------
  # Step 2. Estimate Pearson residuals and deviance residuals. 
  #-----------------------------------------------------------------------------
  
  # Get predicted probability of all CUs being above their lower benchmark (or
  # proportion for binomial regression)
  B_0 <- All_Ests %>% filter(Param=="B_0") %>% pull(Estimate)
  B_1 <- All_Ests %>% filter(Param=="B_1") %>% pull(Estimate)
  predPpnAboveBM <- inv_logit(B_0 + B_1*data$LM_Agg_Abund)
  if (caseStudy=="SCChum"){
    predPpnAboveBM <- inv_logit(B_0 + B_1*data$LM_Agg_Abund*assumed_added_scalar)
  }
  
  
  # Pearson residuals: Eq3.15 https://data.princeton.edu/wws509/notes/c3s8
  # setting n=1 (number of trials at each observation of x)
  
  PearResid <- (SMUlogisticData$ppn - predPpnAboveBM) / sqrt(predPpnAboveBM * 
                                                         (1 - predPpnAboveBM)) 

  # Deviance residual:  Fox et al. 2016, p. 412. setting n=1, number of trials 
  # at each observation of x)
  # To avoid NANs, use ifelse statement in the eqn, as suggested here:
  # https://www.datascienceblog.net/post/machine-learning/interpreting_
  # generalized_linear_models/#:~:text=For%20type%20%3D%20%22pearson%22%2
  # 0%2C,%CB%86f(xi)
  # See also Eq3.16 https://data.princeton.edu/wws509/notes/c3s8 for old version
  
  
  binom.resid <- function(y, mu) {
    y * log( ifelse(y== 0, 1, mu/y)) + (1-y) * log( 
      ifelse(y==1 ,1,(1-mu)/(1-y) ) )  
  }

  DevResid <- sign(SMUlogisticData$ppn - predPpnAboveBM ) * 
    sqrt( -2 * binom.resid(y=SMUlogisticData$ppn, mu=predPpnAboveBM) ) 
  
  # binom.resid.old <- function(y, mu) {
  #   y * log( ifelse(y== 0, 1, y/mu)) + (1-y) * log( 
  #     ifelse(y==1 ,1,(1-y)/(1-mu) ) )  
  # }
  # 
  # DevResid.old <- sign(obsPpnAboveBM - predPpnAboveBM ) * 
  #   sqrt( 2 * binom.resid.old(y=obsPpnAboveBM, mu=predPpnAboveBM) ) 
  
  
  ## Testing. Residuals match output from R using glm objects 
  # Must temporarily set penalty to FALSE before running:
  # Fit_Mod <- glm( ppn ~ SMU_Esc , family = quasibinomial, data=SMUlogisticData)
  # B_0 <- Fit_Mod$coef[1]
  # B_1 <- Fit_Mod$coef[2]
  # obsPpnAboveBM <- SMUlogisticData$ppn
  # predPpnAboveBM <- inv_logit(B_0 + B_1*SMUlogisticData$SMU_Esc)
  # residuals(Fit_Mod, type="pearson")
  # residuals(Fit_Mod, type="deviance")
  
  # Put data for diagnostics in a dataframe for plotting
  diagData <- data.frame(predPppnAboveBM = predPpnAboveBM, 
                         PearResid = PearResid, DevResid=DevResid)
  
  p1 <- ggplot(diagData, aes(predPpnAboveBM, PearResid)) +
    geom_point(size=3) + 
    #geom_smooth(method="lm", formula=y~x) + 
    geom_smooth(method="loess", formula=y~x, span=1) + 
    geom_hline(linetype="dashed", yintercept=0) +
    xlab("Predicted proportions") + ylab("Pearson's Residuals") +
    ggtitle("Pearson's Residuals") +
    theme_classic() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(size = 20)
    ) 
  
  p2 <- ggplot(diagData, aes(predPpnAboveBM, DevResid)) +
    geom_point(size=3) + 
    #geom_smooth(method="lm", formula=y~x) + 
    geom_smooth(method="loess", formula=y~x, span=1) + 
    geom_hline(linetype="dashed", yintercept=0) +
    xlab("Predicted proportions") + ylab("Deviance Residuals") +
    ggtitle("Deviance Residuals") +
    theme_classic() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(size = 20)
    ) 
  p1+p2
  ggsave(p2, file=paste(dir, plotname, ".png", sep=""))
  
  #-----------------------------------------------------------------------------
  # Step 3:
  #  Plot autocorrelation among residuals. Are residuals autocorrelated? 
  #   (Assumption 2)
  #-----------------------------------------------------------------------------
  
  # Extract AR-1 values
  AR1.dev<-acf(DevResid)$acf[2]
  AR1.pear<-acf(PearResid)$acf[2]
  
  # See ggplot.cor function in "helperFunctions.r"
  p3 <- ggplot.corr(data=PearResid, title="Pearsons's residuals") 

  p4 <- ggplot.corr(data=DevResid, title="Deviance Residuals") 
 
  p3+p4
  ggsave(p4, file=paste(dir, plotname, "acf.png", sep=""), width=4.5, height=3.5, units="in")
  
  
  #-----------------------------------------------------------------------------
  # Step 4:
  #  Are any residuals outliers? (i.e., are deviance residuals >2, Assumption 3)
  #   Note, we cannot evaluate influence test statistics such as Cook's 
  #   distance because the hat matrix is not provided in TMB outputs
  #-----------------------------------------------------------------------------
  
  # Observations with a deviance residual in excess of 2 may indicate lack of 
  # fit. 
  # https://data.princeton.edu/wws509/notes/c3s8
  # https://www.statisticshowto.com/what-is-a-standardized-residuals/
  
  # See figures which show that deviance residuals for logistic regression >1
  # https://stats.idre.ucla.edu/stata/webbooks/logistic/chapter3/lesson-3-logistic-regression-diagnostics-2/
  
  # Less strict alternative: deviance residual in excess of 3 indicate outlier, 
  # Ahmad et al. 2011 (https://www.researchgate.net/publication/260368584_Diagnostic_for_Residual_Outliers_using_Deviance_Component_in_Binary_Logistic_Regression)
  
  # maximum absolute residuals values > 2?
  round(max(abs(PearResid)),2)
  
  round(max(abs(DevResid)), 2)
  
  
  #-----------------------------------------------------------------------------
  # Step 5:
  #  Evaluate if sample size is sufficient
  #-----------------------------------------------------------------------------
  
  # As a rule of thumb, Peduzzi et al.(1996) suggests a 
  # minimum of 10 cases with the least frequent outcome for each explanatory 
  # variable (1 in this case). For example, if the expected probabilities are 
  # 0.50 and 0.50 (for 0 and 1, respectively), then the minimum sample size of 
  #a t least (10*1) / 0.50 = 20 to avoid biases in model coefficients. 
  # https://www.statology.org/assumptions-of-logistic-regression/
  
  # Frequency of outcomes (coded to accept proportional as well as 0/1 data)
  
  # Original Freq code:
  Freq <- c(sum(floor(SMUlogisticData$ppn))/length(SMUlogisticData$ppn),
            sum(ceiling(SMUlogisticData$ppn))/length(SMUlogisticData$ppn))
  

  # Added by K.Holt for Fraser Coho case study, in which data is binary (success / failure), not proportion format
  if (caseStudy == "IF Coho") {

    Freq <- c(sum(ceiling(SMUlogisticData$ppn))/length(SMUlogisticData$ppn), 
            1-sum(ceiling(SMUlogisticData$ppn))/length(SMUlogisticData$ppn))
  }  
  
  
  minFreq <- min(Freq)
  minSampleSize <- 10/min(Freq)
  
  sampleSize <- length(SMUlogisticData$ppn)
  
  if (minSampleSize > sampleSize) 
    print("Sample is below minimum suggested levels") else 
      print("Sample size is sufficient")
  
  #-------------------------------------------------------------------------------
  # step 6. 
  # Evaluate statistical significance of model coefficients using Wald's test,
  # where Wald's statistic = coefficient/SE(coefficient))
  # Statistical significant is shown in p-value from model output.
  # Agresti et al. (2006); useful for large sample size, but often not adequate
  # when sample sizes are small
  #-------------------------------------------------------------------------------
  
  signTable <-  All_Ests %>% filter(Param %in% c("B_0", "B_1")) %>% 
    rename(P.value=Pr...z.2..) %>% dplyr::select(Param, Estimate, 
                                                 Std..Error, z.value, 
                                                 P.value)
  
  signTable
  
  #-----------------------------------------------------------------------------
  # Step 7. 
  # Evaluate goodness-of-fit based on ratio of Deviance to the null model 
  #   Is there an improvement in fit with the inclusion of the aggregate 
  #   abundance variable? 
  #-----------------------------------------------------------------------------
  
  # Evaluate goodness of fit by comparing the residual deviance to null 
  # deviance, and evaluating this ratio relative to  a Chi-square distribution 
  # (df = 1, the difference in the number of parameters) 
  # Agresti et al. 2007 (LRT); Ahmad et al. 2011 (definition of deviance as
  # -2LL)
  
  
  Deviance <- sum(DevResid^2)
  
  NullDev <- deviance(glm( SMUlogisticData$ppn ~ 1 , family = 
                             quasibinomial))
  
  # Am confused by this example, that doesnt show 1 - pchisq
  # https://stats.stackexchange.com/questions/166585/pearson-vs-deviance-residuals-in-logistic-regression
  # Deviance = - 2 * LL, and chisqr test statistic = 2 (LL-full-LL-null), p<0.05 sign diff in models presumably?
  # pDRT <- signif(pchisq( q = 2* ( (NullDev/-2) - (Deviance/-2) ), df=1, lower.tail=F), digits=2)
  # Aha! Set lower.tail = F (see here, https://api.rpubs.com/tomanderson_34/lrt)
  
  
  # secr pacakge, LR.test function (type, LR.test to see underlying code)
  # chisq test statistic = 2 * abs (LL1 - LL2), 1 - pchisq(statistic, df=1)
  
  # Now, make sure my use of deviances is correct, 
  # Agresit et al. 2007 Ch4 describe chi2 test statistic = -2(LL0-LL1), 
  # which is Dev0-Dev1 
  
  
  pDRT <- signif(1-pchisq(q= NullDev - Deviance, df=1), digits=2)
  #https://stats.stackexchange.com/questions/6505/likelihood-ratio-test-in-r
  
  names(pDRT) <- c("pDRT")
  pDRT
  # P-value <0.05 a indicates significant improvement in fit with addition of 
  # variable
  
  ## However, some caveats with this approach, and further work required to 
  # evaluate it and alternatives
  # https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-020-01055-2
  
  
  # Note, Roback and Legler 2021 suggest evaluating overall model fit based on 
  # chi-square distribution of the deviance itself, with df=n-p
  
  # values < 0.05 indicate statistically significant evidence for lack of fit
  # See section 6.5.6: https://bookdown.org/roback/bookdown-bysh/ch-logreg.html
  
  
  # p.DevChiSq <- 1-pchisq(q=Deviance, df=length(DevResid)-2)
  # https://stats.stackexchange.com/questions/108995/interpreting-residual-and-null-deviance-in-glm-r
  # Maybe don't use this version: https://online.stat.psu.edu/stat501/lesson/15/15.4
  # p.DevChiSq <- pchisq(q=Deviance, df=length(DevResid)-2)# unless I add this, lower.tail=FALSE) 
  # names(p.DevChiSq) <- c("p.DevChiSq")
  # p.DevChiSq
  
  
  
  # Or equivalently, based on Pearson residuals, where sum of squared 
  # Pearson's residuals is,
  # Pearson <- sum(PearResid^2)
  # p.PearChiSq <- 1 - pchisq(q=Pearson, df=length(PearResid)-2)
  ## p.PearChiSq <- pchisq(q=Pearson, df=length(PearResid)-2)# unless I add this, lower.tail=FALSE) 
  # names(p.PearChiSq) <- c("PearChiSq")
  # p.PearChiSq
  # Section 1.4 of https://www.flutterbys.com.au/stats/tut/tut10.5a.html
  
  # Note, the deviance is defined as the difference of likelihoods 
  # between the fitted model and the saturated model:
  # D = − 2 loglik(^β) + 2 loglik(saturated model), 
  # where loglik(saturated model) = 1 and so is dropped from the equation
  # Portugués et al. 2020 (online resource only)
  # https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html (section 4.7)
  
  
  #-------------------------------------------------------------------------------
  # Step 8. 
  # Evaluate quasi-Rsquared.
  #   What is the ratio of the fit to the null model?
  #-------------------------------------------------------------------------------
  
  # What is ratio of fit to the null model? This is a measure of the
  # strength of the relationship
  quasiR2 <- 1 - Deviance/NullDev
  # This is not the percentage of variance explained by the logistic model, 
  # but rather a ratio indicating how close is the fit to being perfect 
  # or the worst. It is not related to any correlation coefficient
  # https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
  # However, when we forced the y-intercept near zero,  the model actually fits 
  # worse than the null model and quasiR2 can be negative
  names(quasiR2) <- c("quasiR2")
  
  quasiR2
  # For WCVI Chinook, the null model is a better fit to the data than the 
  # fitted model because of the penalty, so the quasiR2 is negative.
  
  
  #-----------------------------------------------------------------------------
  # Step 9. 
  # Evaluate hit ratio from a confusion matrix
  #   What is the classification accuracy of the LRP based on the logistic 
  #   regression?
  #-----------------------------------------------------------------------------
  
  
  # In which years did the model predict aggregate abundances >LRP?
  yHat <- predPpnAboveBM > p
  # In which years were observed aggregate abundances >LRP?
  y <- SMUlogisticData$ppn > p
  
  # Confusion Matrix
  confMat <- table(y, yHat)
  confMat
  
  # What is the accuracy in classifying observed aggregate abundances?
  # Hit ratio = ratio of correct classification
  hitRatio <- sum(diag(confMat))/sum(confMat)
  hitRatio <- round(hitRatio, digits=2)
  names(hitRatio) <- c("hitRatio")
  hitRatio
  
  
   
  
  # The hit matrix will be always biased towards unrealistically good
  # classification rates when computed with the same sample used for
  # fitting the logistic model. An approach based on data-splitting/
  # cross-validation is therefore needed to estimate the hit matrix 
  # unbiasedly  https://bookdown.org/egarpor/PM-UC3M/glm-modsel.html
  
  #-----------------------------------------------------------------------------
  # Output
  out <- list()
  out$pBoxTidwell <- pBoxTidwell
  out$DevResid <- DevResid
  out$PearResid <- DevResid
  out$p1 <- p1
  out$p2 <- p2
  out$p3 <- p3
  out$p4 <- p4
  out$AR1.dev<-AR1.dev
  out$AR1.pear<-AR1.pear
  out$minSampleSize <- minSampleSize
  out$sampleSize <- sampleSize
  out$signTable <- signTable
  out$pDRT <- pDRT
  out$quasiR2 <- quasiR2
  out$confMat <- confMat
  out$hitRatio <- hitRatio

  
  return(out)
  
} # End of Function 1: LRdiagnostics()
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Function 2: LOO_LRdiagnostic, performs leave-one-out cross-validation by 
# iteratively re-running LRP estimation of subsets of data, leaving one data 
# point out each time, and estimating classification accuracy on removed data 
# point
# Purpose: answer the questions, what is the out-of-sample classification 
# accuracy?
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Steps:
# 1. Get time-series length
# 2. Run logistic regression iteratively removing one year in time-series 
# each time
# 3. Get observed data (1 or 0 indicating if all CUs were  above their lower 
# benchmark) for the year that was held out, iteratively (also accepts ppn of 
# CUs above lower benchmark)
# 4. Get predicted probability of all CUs above their lower benchmark (or 
# proportion of all CUs above their lower benchmark) for the year that was held 
# out (iteratively)
# 5. Evaluate confusion matrix hit rate on the resulting time-series of 
# predicted and observed probabilities (or ppns).


# Note, this function will have to be adapted to case studies based on the 
# logistic regression model used for those stocks. I revised the Get.LRP() 
# function for WCVI Chinook by adding an argument LOO, where LOO = numeric for 
# leave-one-out cross validation of the logistic regression. This number is the 
# index of the time-series of ppn of CUs and aggregate abundances that are 
# removed prior to implementing the logistic regression in TMB. It is set to NA 
# as a default (no values removed). Note, the outputted time-series ('out') 
# contain all the data, but parameter estimates are derived from 
# time-series without LOO index value. Within the Get.LRP function, I added 
# the following code when setting up data list for TMB:
# if(!is.na(LOO)) {
#  data$LM_Agg_Abund <- data$LM_Agg_Abund[-LOO]
#  data$N_Above_BM <- data$N_Above_BM[-LOO]
# }



#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Arguments:
# remove.EnhStocks = logical specifying if enhanced stocks are removed 
  # (PNI<0.5)



# Note, LOO_LRdiagnists() calls the function Get.LRP() which generates the 
# following outputs required for calculation of this diagnostic:
# out$All_Ests = Dataframe containing parameters of logistic regression, B_0 and B_1
  # with p values (All_Ests <- data.frame(summary(sdreport(obj), p.value=TRUE)))
# SMU_Esc = Vector of scaled observed aggregate abundances, with NAs removed
# out$Logistic_Data$yy = Vectors of observed ppn of CUs above their lower benchmarks
  # with NAs removed
# LRPppn = the proportion of CUs that must be > their benchmark when defining LRP


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Returns:
# - Hit Ratio from confusion matrix, based on leave-one-out predictions
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Function 2:

LOO_LRdiagnostics <- function(remove.EnhStocks=TRUE){
  
  
  # Step 1: Get length of time-series used in logistic regression, n
  n <- length(which(!is.na(Get.LRP(remove.EnhStocks = TRUE)$SMU_ppn)))[1]
  
  # Step 2: estimate logistic regression iteratively, removing a single year 
  # each time
  
  
  predPpnAboveBM <- NA
  
  for (i in 1:n){
    # Estimate logistic regression using function Get.LRP
    zz <- Get.LRP(remove.EnhStocks = TRUE, LOO=i)
    All_Ests <- zz$out$All_Ests

    if(i==1){ # These remain constant over iterations
      # Step 3: Get observed time-series of aggregate raw abundances that includes all
      # data and then scale to units near 1-10 
      AggAbundRaw <- zz$out$Logistic_Data$xx
      digits <- count.dig(AggAbundRaw)
      ScaleSMU <- min(10^(digits -1 ), na.rm=T)
      AggAbund <- AggAbundRaw/ScaleSMU
      # Get time-series of observed ppns of CUs> benchmark, including all
      # data
      obsPpnAboveBM <- zz$out$Logistic_Data$yy
      # Get threshold p value (ppn of CUs>benchmark) used to estimate LRP
      p <- zz$LRPppn
    }
 
   # Step 4: Get predicted ppn of CUs above their lower benchmark for the year 
   # that was held out
    B_0 <- All_Ests %>% filter(Param=="B_0") %>% pull(Estimate)
    B_1 <- All_Ests %>% filter(Param=="B_1") %>% pull(Estimate)
    predPpnAboveBM[i] <- inv_logit(B_0 + B_1 * AggAbund[i])
  } # End of for i in 1:n
  
  # Step 5: Calculate Hit Ratio
  
  # In which years did the model predict aggregate abundances >LRP?
  yHat <- predPpnAboveBM > p
  # In which years were observed aggregate abundances >LRP?
  y <- obsPpnAboveBM > p
  
  # Confusion Matrix
  confMat <- table(y, yHat)
  
  # What is the accuracy in classifying observed aggregate abundances?
  # Hit ratio = ratio of correct classification
  hitRatio_LOO <- sum(diag(confMat))/sum(confMat)
  hitRatio_LOO <- round(hitRatio_LOO, digits=2)
  
  return(hitRatio=hitRatio_LOO)
}# End of Function 2: LOO_LRdiagnostics()



# 
# #----------------------------------------------------------------------------
# # Run code using WCVI CK input data from above
# # load("DataIn/Input_LRdiagnostics.rda")
# 
# LRdiagOut <- LRdiagnostics(SMUlogisticData = input$SMUlogisticData,
#                            nCU = input$nCU,
#                            All_Ests = input$All_Ests, p = input$p,
#                            Bern_logistic = input$Bern_logistic,
#                            dir = input$dir, plotname = input$plotname)
# 
# LOO_LRdiagOut <- LOO_LRdiagnostics(remove.EnhStocks=TRUE)
# save(LRdiagOut, LOO_LRdiagOut, file="DataOut/LRdiagout.rda")

#-------------------------------------------------------------------------------

# 
# #----------------------------------------------------------------------------
# # Run code using ISC chum input data from above

#CW commented this out in order to run cofo example
# 
#LRdiagOut <- LRdiagnostics(SMUlogisticData = SMUlogisticData,
#                           nCU = nCU,
#                           All_Ests = All_Ests, p = p,
#                           Bern_logistic = Bern_logistic,
#                           dir = dir, plotname = plotname,
#                           caseStudy = "SCChum")
#
#save(LRdiagOut, file="DataOut/logisticFit_2018Output.rda")

#-------------------------------------------------------------------------------

LOO_LRdiagnostics_cohoModel <- function(year, p, useBern_Logistic,
                                        RunName,outputDir, TMB_Inputs) {
  
  # Step 1: Extract logistic regression data from .rda file
  load(paste(outputDir,"/logisticFit_",year,".rda",sep=""))
  print(LRP_Mod$Logistic_Data)
  
  # Save "observed" data from logistic model fit with all data
  # how Carrie scaled aggregate abundance data:
  AggAbundRaw<-LRP_Mod$Logistic_Data$xx
  digits <- count.dig(AggAbundRaw)
  ScaleSMU <- min(10^(digits -1 ), na.rm=T)
  obsAggAbund <- AggAbundRaw/ScaleSMU
  # an alternative that seems to give the same thing ...
  #obsAggAbund<-LRP_Mod$Logistic_Data$xx/TMB_Inputs$Scale
  
  obsPpnAboveBM <- LRP_Mod$Logistic_Data$yy
  
  LRPmodel<-"LRP_Logistic_Only_Bernoulli"
  
  # Loop over nyears and leave-one-out of data set
  predPpnAboveBM <- NA
  
  for (i in 1:nrow(LRP_Mod$Logistic_Data)) {
    
    # Step 2: Get observed time-series of aggregate raw abundances that includes all
    # data and then scale to units near 1-10
    
    # ---- 2a) Set-up data and parameter inputs ---------------------
    Logistic_Data<-LRP_Mod$Logistic_Data[-i,]
    
    # Test: Used following to confirm that logistic only model fit with all data gave the same LRP estimate
    # as the original integrated model (K.Holt, March 2021).Estimates and S.Error match.
    #Logistic_Data<-LRP_Mod$Logistic_Data
    
    Logistic_Data$yr_num <- group_by(as.data.frame(Logistic_Data), yr) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
    
    
    data<-list()
    data$LM_Agg_Status<-Logistic_Data$yy
    data$LM_Agg_Abund <- Logistic_Data$xx / TMB_Inputs$Scale
    data$LM_yr <- Logistic_Data$yr_num
    data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
    data$p <- p
    
    param<-list()
    param$B_0 <- 2
    param$B_1 <- 0.1
    
    
    #  ---- 2b) Run TMB code to fit logistic model ---------------------
    obj <- MakeADFun(data, param, DLL=LRPmodel, silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    
    #  ---- 2c) Create table of outputs ----------------------
    All_Ests <- data.frame(summary(sdreport(obj)))
    All_Ests$Param <- row.names(All_Ests)  
    
    # Put together readable data frame of values
    All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
    All_Ests[All_Ests$Param == "Agg_LRP", ] <-  All_Ests %>% filter(Param == "Agg_LRP") %>% 
      mutate(Estimate = Estimate*TMB_Inputs$Scale) %>% mutate(Std..Error = Std..Error*TMB_Inputs$Scale)
    Preds <- All_Ests %>% filter(Param == "Logit_Preds")
    All_Ests <- All_Ests %>% filter(!(Param == "Logit_Preds")) 
    
    
    # Step 3: Get predicted ppn of CUs above their lower benchmark for the year
    # that was held out
    B_0 <- All_Ests %>% filter(Param=="B_0") %>% pull(Estimate)
    B_1 <- All_Ests %>% filter(Param=="B_1") %>% pull(Estimate)
    # Step 3a: Calculate predicted ppn of CUs above lower benchmark for the year that was held out
    predPpnAboveBM[i] <- inv_logit(B_0 + B_1*obsAggAbund[i])
    
    # Testing plots
    #Preds.plot<-Preds %>% mutate(PredProp=inv_logit(Estimate))
    #plot(data$Pred_Abund,Preds.plot$PredProp)
    
    
  } # End of for i in 1:nyears
  
  # Step 4: Calculate Hit Ratio
  # In which years did the model predict aggregate abundances >LRP?
  yHat <- predPpnAboveBM > p
  
  # In which years were observed aggregate abundances >LRP?
  y <- obsPpnAboveBM > p
  #
  # Confusion Matrix
  confMat <- table(y, yHat)
  
  # What is the accuracy in classifying observed aggregate abundances?
  # Hit ratio = ratio of correct classification
  hitRatio <- sum(diag(confMat))/sum(confMat)
  hitRatio <- round(hitRatio, digits=2)
  
  return(hitRatio)
  

}