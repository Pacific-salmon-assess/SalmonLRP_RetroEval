## April 20234: For WCVI CK, need to change CUpars and mcmc files that are 
# read it, and make sure they are consistent with scenario being run

run_ScenarioProj <- function(SRDat, BMmodel, scenarioName, useGenMean, genYrs,
                             TMB_Inputs, outDir, runMCMC, nMCMC, nProj,
                             ERScalar=NULL, recCorScalar, cvER,
                             gammaSigScalar=NULL, cvERSMU=NULL, agePpnConst=NULL,
                             annualcvERCU=NULL,
                             corMat=NULL, alphaScalar=NULL, aNarrow=NULL, SREPScalar=NULL,
                             biasCorrectEst=NULL, biasCorrectProj=NULL, ER=NULL, evenPars=NULL){


  scenInputDir <- paste(outDir, "SamSimInputs", scenarioName, sep="/")
  scenOutputDir <- paste(outDir, "SamSimOutputs", sep="/")

    if (file.exists(scenInputDir) == FALSE){
      dir.create(scenInputDir)
    }
    if (file.exists(scenOutputDir) == FALSE){
      dir.create(scenOutputDir)
    }

  if(!is.null(SRDat)){
    if (all(is.na(SRDat$Recruits)) == FALSE){
      # Run MPD fit to parameterize samSim projections ==========
      mpdOut<-get_MPD_Fit(SRDat, BMmodel, TMB_Inputs, outDir, biasCorrectEst=biasCorrectEst)
      
      # save correlation matrix
      corMatrix<-mpdOut$corMatrix
      corMatrix<-corMatrix * recCorScalar
      corMatrix[col(corMatrix)==row(corMatrix)] <- 1
      write.table(corMatrix, paste(scenInputDir,"corrMat.csv",sep="/"),row.names=F, col.names=F, sep=",")

      # Run MCMC fit to parameterize samSim projections =======
      if (runMCMC == T) {
        mcmcOut<-get_MCMC_Fit(scenarioName, obj= mpdOut$obj, init=mpdOut$par,
                              upper=mpdOut$upper, lower=mpdOut$lower, nMCMC=nMCMC, Scale = TMB_Inputs$Scale)

        # Need to convert capacity parameters to beta if prior on capacity is used
        if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
          # -- calculate muLSurv (needed to get beta)
          #muSurv <- SRDat %>% group_by(CU_ID) %>%
           # summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
          #muLSurv <- log(muSurv$muSurv)
          # -- expand muLSurv for number of MCMC samples
          #muLSurv<-rep(muLSurv,length=nrow(mcmcOut))
          
          # Alternative method that does not weight by observed returns at age
          muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
          muLSurv <- muLSurv$muLSurv
          
          # -- update mcmcOut
          mcmcOut <- mcmcOut %>% add_column(beta = (mcmcOut$alpha + mcmcOut$gamma * muLSurv) / mcmcOut$cap, .after="alpha") %>% select(-cap)

        }

        mcmcOut<-as.data.frame(mcmcOut)
        
        # Save to high-level outDir folder so that it can be used for future runs with runMCMC=F
        write.csv(mcmcOut, paste(outDir,"/SamSimInputs/", BMmodel,"_mcmc.csv",sep=""),row.names=F)
        # Also save to scenario folder so it can be easily referenced for a run
        write.csv(mcmcOut, paste(scenInputDir,"/", BMmodel,"_mcmc.csv",sep=""),row.names=F)
        
      }

      if (runMCMC == F) {
        # Read-in previoulsy saved mcmcOut
        mcmcOut<-read.csv(paste(outDir,"/SamSimInputs/", BMmodel,"_mcmc.csv", sep=""))
        # Save to scenario folder so it can be easily referenced for a run
        write.csv(mcmcOut, paste(scenInputDir,"/", BMmodel,"_mcmc.csv",sep=""),row.names=F)
      }
    }# End of if all recruitment=NA

    # If there are recruitment data for some ages, but not all, fill in the
    # remaining ages with 0s
    if (all(is.na(SRDat$Recruits)) == FALSE){
      if(is.null(SRDat$Age_2_Recruits)){ rec2<-0} else {rec2 <-
        SRDat$Age_2_Recruits}
      if(is.null(SRDat$Age_3_Recruits)){ rec3<-0} else {rec3 <-
        SRDat$Age_3_Recruits}
      if(is.null(SRDat$Age_4_Recruits)){ rec4<-0} else {rec4 <-
        SRDat$Age_4_Recruits}
      if(is.null(SRDat$Age_5_Recruits)){ rec5<-0} else {rec5 <-
        SRDat$Age_5_Recruits}
      if(is.null(SRDat$Age_6_Recruits)){ rec6<-0} else {rec6 <-
        SRDat$Age_6_Recruits}
    } #End of if (all(is.na(SRDat$Recruits)) == FALSE){


  }# End of if(!is.null(SRDat)){



  #If there are NO recruitment data for any ages but SRDat dataset exists,
  # fill in  with NAs
  if(!is.null(SRDat)){
    if (all(is.na(SRDat$Recruits)) == TRUE) {
      rec2 <- rec3 <- rec4 <- rec5 <- rec6 <- NA
    }
  }

  # If SRDat exists, create recDatTrim for input into projections (for priming)
  if(!is.null(SRDat)){
    recDatTrim<-data.frame(stk=(SRDat$CU_ID + 1),
                         yr=SRDat$BroodYear,
                         ets=SRDat$Spawners,
                         totalSpwn=SRDat$Spawners,
                         rec2=rec2,#rep(0,length(SRDat$BroodYear)),
                         rec3=rec3,#SRDat$Age_3_Recruits,
                         rec4=rec4,#SRDat$Age_4_Recruits,
                         rec5=rec5,#rep(0,length(SRDat$BroodYear)),
                         rec6=rec6)#rep(0,length(SRDat$BroodYear)))
    write.csv(recDatTrim, paste(scenInputDir,"recDatTrim.csv", sep="/"), row.names=F)
  } else {recDatTrim <- NULL}

  
  # If there are no recruitment data, then pull correlation matrix from inputs
  if(is.null(SRDat) || all(is.na(SRDat$Recruits)) ){
    corMatrix <- corMat
    corMatrix <- corMatrix * recCorScalar
    corMatrix[col(corMatrix)==row(corMatrix)] <- 1
    write.table(corMatrix, paste(scenInputDir,"corrMat.csv",sep="/"),
                row.names=F, col.names=F, sep=",")

    if (runMCMC) {
      if(is.null(alphaScalar) & is.null(SREPScalar)){
        if(is.null(aNarrow)) {
          #mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc.csv", sep="/"))
          
          if(is.null(evenPars)) {
            mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc.csv",
            # mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_wEnh_2024.csv",
            # mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_CoreInd.csv",
            # mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_AllExMH.csv",
                                                                sep="/"))
          }        
          if(!is.null(evenPars)){# This mcmc SR set assume even pars across CUs
            if(evenPars==TRUE){
              mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Even_mcmc.csv",
                                        sep="/"))
            }
            if(evenPars=="Same"){
              mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Same_mcmc.csv",
                                        sep="/"))
            }
            if(evenPars=="SameProd"){
              mcmcOut <- read.csv(paste(outDir,"SamSimInputs/SameProd_mcmc.csv",
                                        sep="/"))
            }
            if(evenPars=="SameBeta"){
              mcmcOut <- read.csv(paste(outDir,"SamSimInputs/SameBeta_mcmc.csv",
                                        sep="/"))
            }
            if(evenPars=="SameSREP"){
              mcmcOut <- read.csv(paste(outDir,"SamSimInputs/SameSREP_mcmc.csv",
                                        sep="/"))
            }
            
          }
  
        }# End of  if(is.null(aNarrow) {
        if(!is.null(aNarrow)) {
          if(aNarrow){
            mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_narrow.csv",
                                      sep="/"))

          }
        }
      }#End of if(is.null(alphaScalar) & is.null(SREPScalar)){
      if(!is.null(alphaScalar)&!is.null(SREPScalar)){

        if(alphaScalar==1 & SREPScalar==1){
            mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc.csv",
                                      sep="/"))
        }
        if(alphaScalar==1.5 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1.5_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar==0.5 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar0.5_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar==0.75 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar0.75_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar=="lifeStageModel" & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_lifeStageModel.csv",
                                    sep="/"))
        }
        
        if(alphaScalar==1 & SREPScalar==1.5){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1_SREPScalar1.5.csv",
                                    sep="/"))
        }
        if(alphaScalar==1 & SREPScalar==0.5){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1_SREPScalar0.5.csv",
                                    sep="/"))
        }

      }# End of if(!is.null(alphaScalar)&!is.null(SREPscalar)){
    }# End of if (runMCMC) {

    if (!runMCMC) mcmcOut <- NULL
  }# End of  if(is.null(SRDat) || all(is.na(SRDat$Recruits)) ){


  # Read-in CU pars file and re-write with updated scenario pars =====================
  CUpars<-read.csv(paste(outDir, "SamSimInputs/CUPars.csv",sep="/"))
  # CUpars<-read.csv(paste(outDir, "SamSimInputs/CUPars_wEnh.csv",sep="/"))
  # CUpars<-read.csv(paste(outDir, "SamSimInputs/CUPars_CoreInd.csv",sep="/"))
  # CUpars<-read.csv(paste(outDir, "SamSimInputs/CUPars_AllExMH.csv",sep="/"))
  
  # -- specify ER scenario
  CUpars$cvER <- rep(cvER,length(unique(CUpars$stk)))

  # For Sensitivity analysis assuming even SR pars among CU (WCVI CK)
  if(!is.null(evenPars)){
    if(evenPars==TRUE|evenPars=="Same"){
      CUpars$Sinit <- rep(3384,5)
    }
  }
  
  # -- fill-in MPD fits, only for stocks with SR data
  if(!is.null(SRDat)){
    if (all(is.na(SRDat$Recruits)) == FALSE){
      CUpars$alpha <- mpdOut$All_Ests[grepl("logA", mpdOut$All_Ests$Param), "Estimate" ]

      if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
        CUpars$beta0 <- mpdOut$All_Ests[grepl("B", mpdOut$All_Ests$Param), "Estimate" ]/TMB_Inputs$Scale
      }
      else {
        CUpars$beta0 <- exp(mpdOut$All_Ests[grepl("logB", mpdOut$All_Ests$Param), "Estimate" ])/TMB_Inputs$Scale
      }

      # eliminate logASigma using dum before getting logSigma:
      dum<-mpdOut$All_Ests[mpdOut$All_Ests$Param != "logSigmaA",]
      CUpars$sigma <- exp(dum[grepl("logSigma",dum$Param),"Estimate"])

      CUpars$gamma <- rep(mpdOut$All_Ests[grepl("gamma", mpdOut$All_Ests$Param), "Estimate" ],length(unique(SRDat$CU_ID)))
      # -- add initial marine survival covariate based on recent average, scaled by stock-specific return at age
    
      
      # -- specify distribution for marine survival for annual sampling
      # Do not take out pre-1999 years (like Arbeider et al. 2020 did) so that estimation and projection are consistent
      mean_logMS<-SRDat %>% group_by(CU_ID) %>% summarise(mean_logMS=mean(log(STAS_Age_3)))
      CUpars$mu_logCovar1 <- mean_logMS$mean_logMS
      sig_logMS<-SRDat %>% group_by(CU_ID) %>% summarise(sig_logMS=sd(log(STAS_Age_3)))
      CUpars$sig_logCovar1 <- sig_logMS$sig_logMS
      
      min_logMS<-SRDat %>% group_by(CU_ID) %>% summarise(min_logMS=min(log(STAS_Age_3)))
      CUpars$min_logCovar <- min_logMS$min_logMS # lowest observed mSurv in updated 2020 dataset
      
      max_logMS<-SRDat %>% group_by(CU_ID) %>% summarise(max_logMS=max(log(STAS_Age_3)))
      CUpars$max_logCovar <- max_logMS$max_logMS # lowest observed mSurv in updated 2020 dataset
      
      
      # -- add mean recruitment, by age
      CUpars$meanRec2 <- rep(0,length(unique(SRDat$CU_ID)))
      rec3<-SRDat %>% group_by(CU_ID) %>% summarize(rec3=mean(Age_3_Recruits/Recruits))
      CUpars$meanRec3 <- rec3$rec3
      rec4<-SRDat %>% group_by(CU_ID) %>% summarize(rec4=mean(Age_4_Recruits/Recruits))
      CUpars$meanRec4 <- rec4$rec4
      # -- add median annual recruitment and quantiles
      quantRec<- SRDat %>% group_by(CU_ID) %>% summarize(med=median(Recruits), lowQ=quantile(Recruits,0.25),highQ=quantile(Recruits,0.75))
      CUpars$medianRec <- quantRec$med
      CUpars$lowQRec <- quantRec$lowQ
      CUpars$highQRec <- quantRec$highQ

    }# End of if (all(is.na(SRDat$Recruits)) == FALSE)
  }#  End of  if(!is.null(SRDat)){
  
  write.csv(CUpars, paste(scenInputDir,"CUPars.csv", sep="/"), row.names=F)

  # Read-in sim par file and re-write with updated scenario pars =====================
  simPars<-read.csv(paste(outDir, "SamSimInputs/SimPars.csv",sep="/"))
  simPars$nameOM<-rep(scenarioName,nrow(simPars))
  simPars$scenario<-paste(simPars$nameOM,simPars$nameMP,sep="_")
  #simPars$seed <- 1:nrow(simPars)#rep(1,nrow(simPars))#
  if(!is.null(ER)) simPars$canER <- ER

  if(!is.null(cvERSMU)){
    simPars$cvERSMU <- cvERSMU
  }
  if(!is.null(annualcvERCU)){
    simPars$annualcvERCU <- annualcvERCU
  }
  if(!is.null(agePpnConst)){
    simPars$agePpnConst <- agePpnConst
  }
  # If gammaSigScalar is specified in function call, add to simPars file
  if (is.null(gammaSigScalar)==FALSE) {
  
    if (is.null(mcmcOut) == TRUE) {
      # Use mpd fit standard error if no mcmc outputs available
      gammaSig<-mpdOut$All_Ests[mpdOut$All_Ests$Param=="gamma","Std..Error"]
    } else {
      # Use standard deviation of gamma posterior is mcmc output is available
      gammaSig<-as.numeric(mcmcOut %>% filter(stk==1) %>% summarize(sd(gamma)))
    }
    simPars$sampCU_gamma<-TRUE
    simPars$sigCU_gamma<-gammaSig*gammaSigScalar
  }
  
  # Add BiasCorr to simPars
  if (is.null(biasCorrectProj)) {
    simPars$biasCor <- rep(FALSE,nrow(simPars))
  } else {
    if (biasCorrectProj == TRUE) simPars$biasCor <- rep(TRUE,nrow(simPars))
    if (biasCorrectProj == FALSE) simPars$biasCor <- rep(FALSE,nrow(simPars))
  }
  
  write.csv(simPars, paste(scenInputDir,"SimPars.csv", sep="/"), row.names=F)

  ## Run projections =================================================================================

  ## Check if necessary packages are available and install if necessary
  listOfPackages <- c("here", "parallel", "doParallel", "foreach",
                      "tidyverse", "tictoc", "samSim")
  newPackages <- listOfPackages[!(listOfPackages %in%
                                    installed.packages()[ , "Package"])]
  if (length(newPackages)) {
    install.packages(newPackages)
  }
  lapply(listOfPackages, require, character.only = TRUE)

  dimnames(corMatrix)=NULL

  dirNames <- simPars$nameOM

  simsToRun <- split(simPars, seq(nrow(simPars)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 1) #save one core
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(samSim)))

  clusterExport(cl, c("simsToRun", "CUpars", "nProj",
                       "recDatTrim", "corMatrix", "mcmcOut"),
                envir=environment())

  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    genericRecoverySim(x, cuPar=CUpars, srDat=recDatTrim,
                variableCU=FALSE, ricPars=mcmcOut,
                cuCustomCorrMat = corMatrix,
                nTrials=nProj, makeSubDirs=FALSE,
                random=FALSE, outDir=outDir)
   })
   stopCluster(cl) #end cluster
  toc()

  # for (i in 1:nrow(simPars)) {
  #
  # genericRecoverySim(simPars[i, ], cuPar=CUpars, srDat=recDatTrim,
  #         variableCU=FALSE, ricPars=mcmcOut, cuCustomCorrMat = corMatrix,
  #          nTrials=nProj, makeSubDirs=FALSE, random=FALSE, outDir=outDir)
  #
  # }


  # Read-in projection outputs
  for (i in 1:nrow(simPars)) {

    filename<-paste(simPars[i, "nameOM"],simPars[i, "nameMP"],"CU_SRDat.csv",sep="_" )
    filename2<-paste(simPars[i, "nameOM"],simPars[i, "nameMP"],"lrpDat.csv",sep="_" )
    datCUSp.i<-read.csv(here(outDir,"SamSimOutputs", "simData", dirNames[[i]], filename))
    datCUSp.i$expRate<-simPars[i, "canER"] + simPars[i, "usER"]
    datLRP.i<-read.csv(here(outDir,"SamSimOutputs", "simData", dirNames[[i]], filename2))
    datLRP.i$expRate<-simPars[i, "canER"] + simPars[i, "usER"]

    if (i == 1) {
      projSpwnDat<-datCUSp.i
      projLRPDat<-datLRP.i
    }

    if (i > 1) {
      projSpwnDat<-rbind(projSpwnDat,datCUSp.i)
      projLRPDat<-rbind(projLRPDat,datLRP.i)
    }

  }


  write.csv(projSpwnDat,paste(here(outDir,"SamSimOutputs", "simData"), paste("projSpwnDat_",simPars$nameOM[1],".csv",sep=""),sep="/"))
  write.csv(projLRPDat,paste(here(outDir,"SamSimOutputs", "simData"), paste("projLRPDat_",simPars$nameOM[1],".csv",sep=""),sep="/"))
    
}

#______________________________________________________________________________

# Run a set of scenarios specified in SimPars in parallel. The set is labeled 
# "scenarioName", but the specific OM scenarios are in simPars. Each of the 
# individual scenarios in simPars must use the same cuPars (i.e., same BMmodel 
# and/or mcmc estimates). If testing different BMmodels, then create different 
# scencario sets.

# To run this code, first create SimPars with a list of OMs to test. labeling 
# the names in column B "namesOM"

# !!This won't work unless I can put covariance matrix created with recCorScalar 
# within simPars, so that it is iterated over parallel processing, and move cvER to simPars


run_ScenarioProjParallel <- function(SRDat, BMmodel, scenarioName, useGenMean, 
                                     genYrs, TMB_Inputs, outDir, runMCMC, nMCMC, 
                                     nProj, ERScalar=NULL,
                                     gammaSigScalar=NULL, 
                                     recCorScalar=NULL, 
                                     corMat=NULL, alphaScalar=NULL, 
                                     aNarrow=NULL, SREPScalar=NULL, 
                                     biasCorrectEst=NULL, biasCorrectProj=NULL){
  
  
  scenInputDir <- paste(outDir, "SamSimInputs", scenarioName, sep="/")
  scenOutputDir <- paste(outDir, "SamSimOutputs", sep="/")
  
  if (file.exists(scenInputDir) == FALSE){
    dir.create(scenInputDir)
  }
  if (file.exists(scenOutputDir) == FALSE){
    dir.create(scenOutputDir)
  }
  
  if(!is.null(SRDat)){
    if (all(is.na(SRDat$Recruits)) == FALSE){
      # Run MPD fit to parameterize samSim projections ==========
      mpdOut<-get_MPD_Fit(SRDat, BMmodel, TMB_Inputs, outDir, biasCorrectEst=biasCorrectEst)
      # save correlation matrix
      corMatrix<-mpdOut$corMatrix
      corMatrix<-corMatrix * recCorScalar
      corMatrix[col(corMatrix)==row(corMatrix)] <- 1
      write.table(corMatrix, paste(scenInputDir,"corrMat.csv",sep="/"),row.names=F, col.names=F, sep=",")
      
      # Run MCMC fit to parameterize samSim projections =======
      if (runMCMC == T) {
        mcmcOut<-get_MCMC_Fit(scenarioName, obj= mpdOut$obj, init=mpdOut$par,
                              upper=mpdOut$upper, lower=mpdOut$lower, nMCMC=nMCMC, Scale = TMB_Inputs$Scale)
        
        # Need to convert capacity parameters to beta if prior on capacity is used
        if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
          # -- calculate muLSurv (needed to get beta)
          #muSurv <- SRDat %>% group_by(CU_ID) %>%
          # summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
          #muLSurv <- log(muSurv$muSurv)
          # -- expand muLSurv for number of MCMC samples
          #muLSurv<-rep(muLSurv,length=nrow(mcmcOut))
          
          # Alternative method that does not weight by observed returns at age
          muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
          muLSurv <- muLSurv$muLSurv
          
          # -- update mcmcOut
          mcmcOut <- mcmcOut %>% add_column(beta = (mcmcOut$alpha + mcmcOut$gamma * muLSurv) / mcmcOut$cap, .after="alpha") %>% select(-cap)
          
        }
        
        mcmcOut<-as.data.frame(mcmcOut)
        
        # Save to high-level outDir folder so that it can be used for future runs with runMCMC=F
        write.csv(mcmcOut, paste(outDir,"/SamSimInputs/", BMmodel,"_mcmc.csv",sep=""),row.names=F)
        # Also save to scenario folder so it can be easily referenced for a run
        write.csv(mcmcOut, paste(scenInputDir,"/", BMmodel,"_mcmc.csv",sep=""),row.names=F)
        
      }
      
      if (runMCMC == F) {
        # Read-in previoulsy saved mcmcOut
        mcmcOut<-read.csv(paste(outDir,"/SamSimInputs/", BMmodel,"_mcmc.csv", sep=""))
        # Save to scenario folder so it can be easily referenced for a run
        write.csv(mcmcOut, paste(scenInputDir,"/", BMmodel,"_mcmc.csv",sep=""),row.names=F)
      }
    }# End of if all recruitment=NA
    
    # If there are recruitment data for some ages, but not all, fill in the
    # remaining ages with 0s
    if (all(is.na(SRDat$Recruits)) == FALSE){
      if(is.null(SRDat$Age_2_Recruits)){ rec2<-0} else {rec2 <-
        SRDat$Age_2_Recruits}
      if(is.null(SRDat$Age_3_Recruits)){ rec3<-0} else {rec3 <-
        SRDat$Age_3_Recruits}
      if(is.null(SRDat$Age_4_Recruits)){ rec4<-0} else {rec4 <-
        SRDat$Age_4_Recruits}
      if(is.null(SRDat$Age_5_Recruits)){ rec5<-0} else {rec5 <-
        SRDat$Age_5_Recruits}
      if(is.null(SRDat$Age_6_Recruits)){ rec6<-0} else {rec6 <-
        SRDat$Age_6_Recruits}
    } #End of if (all(is.na(SRDat$Recruits)) == FALSE){
    
    
  }# End of if(!is.null(SRDat)){
  
  
  
  #If there are NO recruitment data for any ages but SRDat dataset exists,
  # fill in  with NAs
  if(!is.null(SRDat)){
    if (all(is.na(SRDat$Recruits)) == TRUE) {
      rec2 <- rec3 <- rec4 <- rec5 <- rec6 <- NA
    }
  }
  
  # If SRDat exists, create recDatTrim for input into projections (for priming)
  if(!is.null(SRDat)){
    recDatTrim<-data.frame(stk=(SRDat$CU_ID + 1),
                           yr=SRDat$BroodYear,
                           ets=SRDat$Spawners,
                           totalSpwn=SRDat$Spawners,
                           rec2=rec2,#rep(0,length(SRDat$BroodYear)),
                           rec3=rec3,#SRDat$Age_3_Recruits,
                           rec4=rec4,#SRDat$Age_4_Recruits,
                           rec5=rec5,#rep(0,length(SRDat$BroodYear)),
                           rec6=rec6)#rep(0,length(SRDat$BroodYear)))
    write.csv(recDatTrim, paste(scenInputDir,"recDatTrim.csv", sep="/"), row.names=F)
  } else {recDatTrim <- NULL}
  
  # If there are no recruitment data, then pull correlation matrix from inputs
  if(is.null(SRDat) || all(is.na(SRDat$Recruits)) ){
    corMatrix <- corMat
    corMatrix <- corMatrix * recCorScalar
    corMatrix[col(corMatrix)==row(corMatrix)] <- 1
    write.table(corMatrix, paste(scenInputDir,"corrMat.csv",sep="/"),
                row.names=F, col.names=F, sep=",")
    
    if (runMCMC) {
      if(is.null(alphaScalar) & is.null(SREPScalar)){
        if(is.null(aNarrow)) {
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc.csv",
                                    sep="/"))
        }# End of  if(is.null(aNarrow) {
        if(!is.null(aNarrow)) {
          if(aNarrow){
            mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_narrow.csv",
                                      sep="/"))
            
          }
        }
      }#End of if(is.null(alphaScalar) & is.null(SREPScalar)){
      if(!is.null(alphaScalar)&!is.null(SREPScalar)){
        
        if(alphaScalar==1 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc.csv",
                                    sep="/"))
        }
        if(alphaScalar==1.5 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1.5_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar==1.25 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1.25_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar==0.5 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar0.5_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar==0.75 & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar0.75_SREPScalar1.csv",
                                    sep="/"))
        }
        if(alphaScalar=="lifeStageModel" & SREPScalar==1){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_lifeStageModel.csv",
                                    sep="/"))
        }
        
        if(alphaScalar==1 & SREPScalar==1.5){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1_SREPScalar1.5.csv",
                                    sep="/"))
        }
        if(alphaScalar==1 & SREPScalar==0.5){
          mcmcOut <- read.csv(paste(outDir,"SamSimInputs/Ricker_mcmc_alphaScalar1_SREPScalar0.5.csv",
                                    sep="/"))
        }
        
      }# End of if(!is.null(alphaScalar)&!is.null(SREPscalar)){
    }# End of if (runMCMC) {
    
    if (!runMCMC) mcmcOut <- NULL
  }
  
  
  # Read-in CU pars file and re-write with updated scenario pars =====================
  CUpars<-read.csv(paste(outDir, "SamSimInputs/CUPars.csv",sep="/"))
  
  #_____________________________________________________________________
  # Can I move this to simPars???
  # CUpars$cvER <- rep(cvER,length(unique(CUpars$stk)))
  #_____________________________________________________________________
  
  # -- fill-in MPD fits, only for stocks with SR data
  if(!is.null(SRDat)){
    if (all(is.na(SRDat$Recruits)) == FALSE){
      CUpars$alpha <- mpdOut$All_Ests[grepl("logA", mpdOut$All_Ests$Param), "Estimate" ]
      
      if (BMmodel %in% c("SR_IndivRicker_SurvCap","SR_HierRicker_SurvCap")) {
        CUpars$beta0 <- mpdOut$All_Ests[grepl("B", mpdOut$All_Ests$Param), "Estimate" ]/TMB_Inputs$Scale
      }
      else {
        CUpars$beta0 <- exp(mpdOut$All_Ests[grepl("logB", mpdOut$All_Ests$Param), "Estimate" ])/TMB_Inputs$Scale
      }
      
      # eliminate logASigma using dum before getting logSigma:
      dum<-mpdOut$All_Ests[mpdOut$All_Ests$Param != "logSigmaA",]
      CUpars$sigma <- exp(dum[grepl("logSigma",dum$Param),"Estimate"])
      
      CUpars$gamma <- rep(mpdOut$All_Ests[grepl("gamma", mpdOut$All_Ests$Param), "Estimate" ],length(unique(SRDat$CU_ID)))
      # -- add initial marine survival covariate based on recent average, scaled by stock-specific return at age
     
      
      # -- specify distribution for marine survival for annual sampling
      # Do not take out pre-1999 years (like Arbeider et al. 2020 did) so that estimation and projection are consistent
      mean_logMS<-SRDat %>% group_by(CU_ID) %>% summarise(mean_logMS=mean(log(STAS_Age_3)))
      CUpars$mu_logCovar1 <- mean_logMS$mean_logMS
      sig_logMS<-SRDat %>% group_by(CU_ID) %>% summarise(sig_logMS=sd(log(STAS_Age_3)))
      CUpars$sig_logCovar1 <- sig_logMS$sig_logMS
      CUpars$min_logCovar <- rep(log(0.000128),length(unique(SRDat$CU_ID))) # lowest observed mSurv in updated 2020 dataset
      CUpars$max_logCovar <- rep(log(0.036),length(unique(SRDat$CU_ID))) # highest observed mSurv in updated 2020 dataset
      
      # -- add mean recruitment, by age
      CUpars$meanRec2 <- rep(0,length(unique(SRDat$CU_ID)))
      rec3<-SRDat %>% group_by(CU_ID) %>% summarize(rec3=mean(Age_3_Recruits/Recruits))
      CUpars$meanRec3 <- rec3$rec3
      rec4<-SRDat %>% group_by(CU_ID) %>% summarize(rec4=mean(Age_4_Recruits/Recruits))
      CUpars$meanRec4 <- rec4$rec4
      # -- add median annual recruitment and quantiles
      quantRec<- SRDat %>% group_by(CU_ID) %>% summarize(med=median(Recruits), lowQ=quantile(Recruits,0.25),highQ=quantile(Recruits,0.75))
      CUpars$medianRec <- quantRec$med
      CUpars$lowQRec <- quantRec$lowQ
      CUpars$highQRec <- quantRec$highQ
      
    }# End of if (all(is.na(SRDat$Recruits)) == FALSE)
  }#  End of  if(!is.null(SRDat)){
  
  write.csv(CUpars, paste(scenInputDir,"CUPars.csv", sep="/"), row.names=F)
  
  # Read-in sim par file and re-write with updated scenario pars =====================
  simPars<-read.csv(paste(outDir, "SamSimInputs/SimPars.csv",sep="/"))
  # Keep the nameOM as in simPars
  #simPars$nameOM<-rep(scenarioName,nrow(simPars))
  simPars$scenario<-paste(simPars$nameOM,simPars$nameMP,sep="_")
  simPars$seed <- rep(1,nrow(simPars))#1:nrow(simPars)
  
  # #_____________________________________________________________________
  # # Specify this in SimPars,not here...
  # if(!is.null(cvERSMU)){
  #   simPars$cvERSMU <- cvERSMU
  # }
  # if(!is.null(annualcvERCU)){
  #   simPars$annualcvERCU <- annualcvERCU
  # }
  # if(!is.null(agePpnConst)){
  #   simPars$agePpnConst <- agePpnConst
  # }
  # # Can I move this to simPars??? But then I have to change samSim to get it from simPar instead of cuPar (with an option to look in cupar for KH?)
  # CUpars$cvER <- rep(cvER,length(unique(CUpars$stk)))
  #_____________________________________________________________________
  # This won't work unless I can put covariance matrix created with recCorScalar within simPars, so that it is iterated over parallel processing
  
  
  # If gammaSigScalar is specified in function call, add to simPars file
  if (is.null(gammaSigScalar)==FALSE) {
    if (is.null(mcmcOut) == TRUE) {
      # Use mpd fit standard error if no mcmc outputs available
      gammaSig<-mpdOut$All_Ests[mpdOut$All_Ests$Param=="gamma","Std..Error"]
    } else {
      # Use standard deviation of gamma posterior is mcmc output is available
      gammaSig<-as.numeric(mcmcOut %>% filter(stk==1) %>% summarize(sd(gamma)))
    }
    simPars$sampCU_gamma<-TRUE
    simPars$sigCU_gamma<-gammaSig*gammaSigScalar
  }
  
  # Add BiasCorr to simPars
  if (is.null(biasCorrectProj)) {
    simPars$biasCor <- rep(FALSE,nrow(simPars))
  } else {
    if (biasCorrectProj == TRUE) simPars$biasCor <- rep(TRUE,nrow(simPars))
    if (biasCorrectProj == FALSE) simPars$biasCor <- rep(FALSE,nrow(simPars))
  }
  
  write.csv(simPars, paste(scenInputDir,"SimPars.csv", sep="/"), row.names=F)
  
  ## Run projections =================================================================================
  
  ## Check if necessary packages are available and install if necessary
  listOfPackages <- c("here", "parallel", "doParallel", "foreach",
                      "tidyverse", "tictoc", "samSim")
  newPackages <- listOfPackages[!(listOfPackages %in%
                                    installed.packages()[ , "Package"])]
  if (length(newPackages)) {
    install.packages(newPackages)
  }
  lapply(listOfPackages, require, character.only = TRUE)
  
  dimnames(corMatrix)=NULL
  
  dirNames <- simPars$nameOM
  
  simsToRun <- split(simPars, seq(nrow(simPars)))
  Ncores <- detectCores()
  cl <- makeCluster(Ncores - 1) #save one core
  registerDoParallel(cl)
  clusterEvalQ(cl, c(library(samSim)))
  
  clusterExport(cl, c("simsToRun", "CUpars", "nProj",
                      "recDatTrim", "corMatrix", "mcmcOut"),
                envir=environment())
  
  tic("run in parallel")
  parLapply(cl, simsToRun, function(x) {
    genericRecoverySim(x, cuPar=CUpars, srDat=recDatTrim,
                       variableCU=FALSE, ricPars=mcmcOut,
                       cuCustomCorrMat = corMatrix,
                       nTrials=nProj, makeSubDirs=FALSE,
                       random=TRUE, outDir=outdir)
  })
  stopCluster(cl) #end cluster
  toc()
  
  # for (i in 1:nrow(simPars)) {
  #
  # genericRecoverySim(simPars[i, ], cuPar=CUpars, srDat=recDatTrim,
  #         variableCU=FALSE, ricPars=mcmcOut, cuCustomCorrMat = corMatrix,
  #          nTrials=nProj, makeSubDirs=FALSE, random=FALSE, outDir=outDir)
  #
  # }
  
  
  # Read-in projection outputs
  for (i in 1:nrow(simPars)) {
    
    filename<-paste(simPars[i, "nameOM"],simPars[i, "nameMP"],"CU_SRDat.csv",sep="_" )
    filename2<-paste(simPars[i, "nameOM"],simPars[i, "nameMP"],"lrpDat.csv",sep="_" )
    datCUSp.i<-read.csv(here(outDir,"SamSimOutputs", "simData", dirNames[[i]], filename))
    datCUSp.i$expRate<-simPars[i, "canER"] + simPars[i, "usER"]
    datLRP.i<-read.csv(here(outDir,"SamSimOutputs", "simData", dirNames[[i]], filename2))
    datLRP.i$expRate<-simPars[i, "canER"] + simPars[i, "usER"]
    
    if (i == 1) {
      projSpwnDat<-datCUSp.i
      projLRPDat<-datLRP.i
    }
    
    if (i > 1) {
      projSpwnDat<-rbind(projSpwnDat,datCUSp.i)
      projLRPDat<-rbind(projLRPDat,datLRP.i)
    }
    
  }
  
  
  write.csv(projSpwnDat,paste(here(outDir,"SamSimOutputs", "simData"), paste("projSpwnDat_",simPars$nameOM[1],".csv",sep=""),sep="/"))
  write.csv(projLRPDat,paste(here(outDir,"SamSimOutputs", "simData"), paste("projLRPDat_",simPars$nameOM[1],".csv",sep=""),sep="/"))
  
}

#______________________________________________________________________________
get_MPD_Fit<-function (SRDat, BMmodel, TMB_Inputs, outDir, biasCorrectEst) {

  # Run MCMC SR analysis to parameterize projections using available data

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
  data$Pred_Spwn <- rep(seq(0,max(data$S)*1.1,length=100), N_Stocks) # vectors of spawner abundance to use for predicted recruits, one vector of length 100 for each stock


  for (i in 1:N_Stocks) { # make a vector of same length as Pred_Spwn with CU ids
    if (i == 1) {
      data$stk_predS <- rep(unique(SRDat$CU_ID)[1],100)
    } else {
      data$stk_predS<-c(data$stk_predS,rep(unique(SRDat$CU_ID)[i],100))
    }
  }

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

  # Calculate correlation matrix in MPD recruitment residuals ========================
  resids<-as_tibble(data.frame(stock=data$stk,year=data$yr, resid=obj$report()$R_Resid))

  cor_matrix <- resids %>% pivot_wider(names_from = stock, names_prefix="stock",values_from = resid)  %>%
    select(-year) %>% cor()

  mpdOut<-list()
  mpdOut$obj<-obj
  mpdOut$par<-opt$par
  mpdOut$upper<-upper
  mpdOut$lower<-lower
  mpdOut$All_Ests<-All_Ests
  mpdOut$corMatrix<-cor_matrix

 mpdOut

}


get_MCMC_Fit<-function (scenarioName, obj, init, upper, lower, nMCMC, Scale) {


# # Fit mcmc with STAN to get parameter estimates for projections ===============================
fitmcmc <- tmbstan(obj, chains=3, iter=nMCMC, init=init, thin=5,
                   control = list(adapt_delta = 0.99),upper=upper, lower=lower)


## Can get ESS and Rhat from rstan::monitor
mon <- monitor(fitmcmc)
max(mon$Rhat)
min(mon$Tail_ESS)

# pull out posterior vals
post<-as.matrix(fitmcmc)
post<-as_tibble(post)
post<-post %>% add_column(iteration=as.numeric(row.names(post)))

# Extract alpha marginal posterior
post_long_alpha<-post %>% select(starts_with("logA"), iteration) %>% pivot_longer(starts_with("logA"),names_to="stock", values_to="logA")
post_long_alpha$stock<-rep(1:5,length=nrow(post_long_alpha))


## Extract Sgen marginal posterior
post_long_Sgen<-post %>% select(starts_with("logSgen"), iteration) %>% pivot_longer(starts_with("logSgen"),names_to="stock", values_to="logSgen")
post_long_Sgen$stock<-rep(1:5,length=nrow(post_long_Sgen))

# Extract beta marginal posterior (or cap parameter if Prior Cap model formulation)
if ("logB" %in% names(obj$par)) {
  post_long_beta<-post %>% select(starts_with("logB"), iteration) %>% pivot_longer(starts_with("logB"),names_to="stock", values_to="logB")
  post_long_beta$stock<-rep(1:5,length=nrow(post_long_beta))
}
if ("cap" %in% names(obj$par)) {
  post_long_cap<-post %>% select(starts_with("cap"), iteration) %>% pivot_longer(starts_with("cap"),names_to="stock", values_to="cap")
  post_long_cap$stock<-rep(1:5,length=nrow(post_long_cap))
}

# Extract sigma marginal posterior
post_long_sigma<-post %>% select(starts_with("logSigma") & !starts_with("logSigmaA"), iteration, gamma) %>%
  pivot_longer(starts_with("logSigma"),names_to="stock", values_to="logSigma")
post_long_sigma$stock<-rep(1:5,length=nrow(post_long_sigma))

# Extract gamma marginal posterior; I have based on an extraction of alpha to replicate gammas over all stocks
dum<-post %>% select(starts_with("logA"), iteration, gamma) %>% pivot_longer(starts_with("logA"),names_to="stock", values_to="logA")
dum <- dum %>% select(-"logA")
post_long_gamma<-dum
post_long_gamma$stock<-rep(1:5,length=nrow(post_long_gamma))

# Compile marginal posteriors to get joint posterior for export
if ("cap" %in% names(obj$par)) {
  post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(cap = post_long_cap$cap*Scale, sigma=exp(post_long_sigma$logSigma), gamma = post_long_gamma$gamma, Sgen=exp(post_long_Sgen$logSgen)*Scale)
 # post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(cap = post_long_cap$cap*Scale, sigma=exp(post_long_sigma$logSigma), gamma = post_long_gamma$gamma)
}

if ("logB" %in% names(obj$par)) {
  post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(beta = exp(post_long_beta$logB)/Scale, sigma=exp(post_long_sigma$logSigma), gamma = post_long_gamma$gamma, Sgen=exp(post_long_Sgen$logSgen)*Scale)
#  post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(beta = exp(post_long_beta$logB)/Scale, sigma=exp(post_long_sigma$logSigma), gamma = post_long_gamma$gamma)
}

post_long

}


