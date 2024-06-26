# Read in raw Chum data, infill and prepare for LRP Analysis, save as csv files

# Reads in raw wild escapement data and infills to get spawners, reads in wild 
# return data (from run reconstruction done by Pieter Van Will) and makes brood
# table based on returns and age composition data from fishery.

# Note that the run reconstruction to get wild returns by CU and area is done
# in an xlsx file separately by Pieter Van Will, not in an R code.

# Also plot data to examine infilling and saves figures.

library(ggplot2)
library(dplyr)
library(tidyr)
options(scipen = 1000000)
# ===========================================================================================
# Read-in Chum Data & Run Data Prep Functions:  
# ===========================================================================================

# Note: The following infilling & brood table reconstruction procedures are taken from: 
# Holt, C.A., Davis, B., Dobson, D., Godbout, L., Luedke, W., Tadey, J., and Van Will, P. 2018.
# Evaluating Benchmarks of Biological Status for Data-limited Conservation Units of Pacific
# Salmon, Focusing on Chum Salmon in Southern BC. DFO Can. Sci. Advis. Sec. Res. Doc.
# 2018/011. ix + 77 p. Available at: https://waves-vagues.dfo-mpo.gc.ca/Library/40759386.pdf

# All code in this section (including code in chumDataFunctions.r) was written by B. Davis (DFO) for the  
# above paper, and provided to us by Carrie Holt in October 2020 as part of "Retrospective Analysis BD" folder.

# With changes by Luke Warkentin

# Read in chum functions
source("R/chumDataFunctions.r")

# Read in raw escapement data from Pieter Van Will. 
rawdat <- readxl::read_excel("DataIn/Chum Escapement Data With Areas(CleanedFeb152021).xlsx", sheet="Updated 2018", trim_ws=TRUE) # strip.white for leading and trailing white spaces in Source and SummerRun columns
# save data summary table for sharing with Island Marine Aquatic Working Group
share <- rawdat %>% group_by(CU_Name, NME, Area, SummerRun, Rabcode) %>% summarise(n=n()) 
write.csv(share[ , -which(names(share)=="n")], "DataOut/streams_counted.csv", row.names = FALSE)

# ----------------------------------------------------#
# Notes on data (From Pieter Van Will): 
# ----------------------------------------------------#
# I. Data variables
#
#   Source:
#     -Wild are wild spawners
#     -Rack is the # of fish that are harvested at a facility and do not contribute to the system 
#           (part of the total return but not contributing as spawners or enhanced)
#     -Enhanced was an assignment we gave given the typically high magnitude of the enhancement, if it was 
#           a large facility on the system with a large chum production we tended to call the system as 
#           a whole enhanced (like Puntledge or the systems with large spawning channels like Qualicum)
#     -Brood are the fish from the system that are taken to support the supplementation of that stock 
#           (do contribute to the subsequent returns).
#
#   SummerRun: We removed the summer run fish as all the data in regards the reconstruction work is 
#       associated with Fall timed stocks.
#
#   NME: Individual streams. Note for Qualicum River, Little Qualicum River, Puntledge River - Historically we assume 
#       these three stocks are 100% enhanced at least since enhancement began at those locations.  
#       We have little data in the enhanced contribution found in the returns but for the purposes of 
#       pulling out wild we make the assumption they were 100% enhanced and not included.
#
# 
# II. Other notes
#
#   Note that NAs are unobserved creeks.
#
#  	Regarding blank escapement count values vs. 0 values, specifically where wild was NA but brood was 0:
#     Pieter Van Will: I would consider any of those as Blanks and would require infill 
#     (perhaps if the Escapement is 0, convert that to a blank and I doubt there any real 0 
#     escapement numbers, even if there are some fish or a 0 in brood stock.)
#
#   Little Qualicum, Qualicum, and Puntledge Rivers: assume these are 100% enhanced, and 
#     remove entirely for infilling for wild escapement

# ----------------------------------------------------#
# Step 1: Infill wild escapement data
# ----------------------------------------------------#
# Infill wild spawners only, and remove Little Qualicum, 
#               Qualicum, and Puntledge Rivers entirely

# Create look-up table for CU names
CU_raw <- unique(rawdat$CU_Name)
CU_short <- c("SCS", "NEVI", "UK", "LB", "BI", "GS", "HSBI")
CU_names<-c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight",
            "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound to Burrard Inlet" )
CUdf <- data.frame(CU_short, "CU_raw"=CU_raw[1:7], CU_names)

# Get number of summer vs fall run streams for each CU (for reference)
table(rawdat$CU_Name, rawdat$SummerRun)
CUsum <- merge(rawdat, CUdf, by.x="CU_Name", by.y="CU_raw")
CUsum1 <- CUsum %>% group_by(CU_Name, CU_short, SummerRun) %>% summarise("Streams" = n_distinct(NME)) %>% 
  pivot_wider(names_from=SummerRun, values_from=Streams)
names(CUsum1) <- c("CU Name", "Abbr.", "Fall run streams", "Summer run streams") # summarise number of streams
CUsum1$`CU Name` <- substr(CUsum1$`CU Name`, 5, 40) # remove leading number, spaces and dash
CUsum1$`Summer run streams`[is.na(CUsum1$`Summer run streams`)] <- 0
write.csv(CUsum1, "DataOut/CU_summary_table_report.csv", row.names=FALSE)

# Process and filter data
rawdat_f <- rawdat[rawdat$SummerRun==FALSE, ] # Remove summer run fish - earlier run timing means they are not intercepted in same fisheries, can't do run reconstruction with them
# Remove non wild fish
rawdat_w <- rawdat_f[rawdat_f$Source=="Wild", ]
# Remove Little Qualicum, Qualicum, and Puntledge Rivers, because they are historically highly enhanced
rawdat_w2 <- rawdat_w[which(rawdat_w$NME %notin% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER")), ]

# wide to long format. maintain NA values (uncounted streams)
ldat <- rawdat_w2 %>% pivot_longer(cols=grep("[[:digit:]]{4}", names(rawdat_f)), names_to="Year", values_to="Escape")

# summarise by stream, to collapse any brood/rack/enhanced/rack etc categories
ldat_s <- ldat %>% group_by(CU_Name, GroupName, GU_Name, NME, SummerRun, Rabcode, Area, Year, Source) %>%
  summarise(Escape=sum(Escape, na.rm=TRUE))
# for sum escapement values that are 0, make back into NA (unobserved)
# since na.rm was TRUE, these became 0 values; summing NA values with na.rm=TRUE returns 0
ldat_s$Escape[ ldat_s$Escape == 0 ] <- NA

#Now Infill
#Now Infill. Gives a list of 2 data frames. First is infilling by stream, second is a summary that is summarized by CU. 
wild_infill_by_stream_list <- Infill(data = ldat_s, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")
#AllData
#write.csv(NoQPDat[[1]], "DataOut/InfilledNoQP_all.csv")
#Sumamarised data
#write.csv(NoQPDat[[2]], "DataOut/InfilledNoQP.csv")

# Remove Fraser, make data frame to feed into infill again
wild_infill_by_stream_sum_no_fraser <- as.data.frame(wild_infill_by_stream_list[[2]][which(wild_infill_by_stream_list[[2]]$CU_Name %in% CUdf$CU_raw),])
wild_infill_by_stream <- wild_infill_by_stream_list[[1]][which(wild_infill_by_stream_list[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]

# Now infill missing years for entire CU

#Infill by CU for years + CU combinations where there are no observations (Knight and Bute)
wild_infill_by_CU <- Infill(data=wild_infill_by_stream_sum_no_fraser, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")

## Also want all infilled by site
wild_infill_join <- left_join(wild_infill_by_stream, data.frame(CU_Name=wild_infill_by_CU[[1]]$CU_Name, Year=wild_infill_by_CU[[1]]$Year, CUEsc = wild_infill_by_CU[[1]]$SiteEsc))
wild_infill_join$Escape <- ifelse(is.nan(wild_infill_join$SiteEsc), wild_infill_join$CUEsc*wild_infill_join$Props, wild_infill_join$SiteEsc)
wild_infill_join$CU <- CUdf$CU_short[match(wild_infill_join$CU_Name, CUdf$CU_raw)]
write.csv(wild_infill_join, "DataOut/wild_spawners_stream_infilled_by_site_CU.csv")

# Save data by CU
write.csv(wild_infill_by_CU[[1]], "DataOut/wild_spawners_CU_infilled_by_site_CU.csv")

# Compare to infilling for run reconstruction
# # Read in what was done for run reconstruction
# zz <- read.csv( "DataOut/infill_escapement_for_external_run_reconstruction/wild_spawners_infilled_by_site.csv")
# CUs <- unique(wild_infill_by_stream$CU_Name)
# compare_stream_infilling_fig <- function(CU, dat1, dat2) {
#   d1 <- dat1[dat1$CU_Name==CU, ]
#   d2 <- dat2[dat2$CU_Name==CU, ]
#   ggplot(d1, aes(y=SiteEsc, x=Year)) +
#     geom_point() +
#     geom_point(data=d2, aes(y=SiteEsc, x=as.integer(Year)), colour="dodgerblue", shape=1, stroke=1.2) +
#     facet_wrap(~NME, scales="free_y") +
#     scale_x_discrete(breaks=seq(1960,2020,10)) +
#     ggtitle(CU) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle=90, vjust=0.5))
# }
# # for total spawners
# fig_list <- as.list(CUs)
# fig_list <- purrr::map(fig_list, compare_stream_infilling_fig, dat1=zz, dat2=wild_infill_by_stream)
# fig_list[[1]] # check
# pdf("Figures/fig_compare_infill_run_reconstruction_brood_table.pdf", width=18, height=12, pointsize=8)
# for (i in 1:length(fig_list)){
#   print(fig_list[[i]])
# }
# dev.off()

# Looks identical 


# ----------------------------------------------------#
# Step 2: Construct Spawner Recruit Brood Tables 
# ----------------------------------------------------#

# Prep infilled escapement data
WildEsc <- wild_infill_by_CU[[1]]
WildEsc$CUinfill <- ifelse(is.na(WildEsc$Escape), TRUE, FALSE) # Flag 'Escape = NA' sites

# Read in estimated wild return data (wild spawners + wild harvest, by year) from Pieter Van Will received 2021-02-24
WildRetWide <- readxl::read_excel("DataIn/wild_ISC_chum_recruitment_PieterVanWill.xlsx", range="A84:BP101", trim_ws = TRUE)
names(WildRetWide)[2] <- "Area"

# Need to get PVW return data into long form:
# First need to "collapse" areas to CU level
WildByCU <- WildRetWide  %>% group_by(CU_Name)   %>%  summarise_at(vars("2018":"1953"),sum,na.rm = T)
WildByCU$CU_Name <- as.character(WildByCU$CU_Name)
# Now into long form
WildRetLong <- WildByCU %>% pivot_longer(cols=grep("[[:digit:]]{4}", names(WildByCU)), names_to="Year", values_to="Return")
# Read in age comp data
ACdat <- read.csv("DataIn/AgeComp_2018.csv")

# Construct Brood table

#Will need to merge EscDat and ECdat to get catch info, make three Brood tables
# Year needs to be integer to join
WildRetLong$Year <- as.integer(WildRetLong$Year)
Btable1 <- left_join(data.frame(Year=as.integer(WildEsc$Year), CU=WildEsc$CU_Name, Escape=WildEsc$SiteEsc, CUinfill=WildEsc$CUinfill), 
                     data.frame(Year=WildRetLong$Year, CU=WildRetLong$CU_Name, Return=WildRetLong$Return), by=c("Year", "CU"))

Btable <- left_join(Btable1, ACdat, by=c("Year"))

years <- sort(unique(Btable$Year), decreasing=F)
nyears <- length(years)
sites <- unique(Btable$CU)
nsites <- length(sites)

#need to enter age comp ages
ages<-c(3,4,5,6)
nages<-length(ages)

# go through years and calculate recruits by brood year
# Same age comp used for each CU
# will not get recruit estimate for first two years due to missing age comp data (NAs)
# will not get recruit estimate for last 
Btable$Recruit <- rep(NA, dim(Btable)[1])
for( i in 1:nsites){
  Sdat<- Btable[which(Btable$CU==sites[i]),]
  #cannot estimate returns for last 5 years
  for( j in 3:nyears ){
    # can only calculate up to nyears-6 -- only have age comps up to 2018
    if(j<=nyears-6){
      Rsum <- 0
      for( k in 1:nages ){
        #add up recruits from brood years
        Rsum <- Rsum + Sdat$Return[which(Sdat$Year==(years[j]+ages[k]))] * Sdat[which(Sdat$Year==years[j]+ages[k]), paste("Age", ages[k], sep="")]
      }
    } else {
      Rsum<-NA
    }
    Btable$Recruit[which(Btable$Year==years[j] & Btable$CU==sites[i])] <- Rsum
  }
}

write.csv(Btable, "DataOut/SRdatWild.csv", row.names = F)

# Check exploitation rate (total vs. wild)
# Read in total harvest (wild + hatchery)
harvest_total <- readxl::read_excel("DataIn/wild_ISC_chum_recruitment_PieterVanWill.xlsx", range="A22:BP39", trim_ws = TRUE)
names(harvest_total)[1] <- "CU" # change name

# Read in escapement (wild + hatchery)
escape_total <- readxl::read_excel("DataIn/wild_ISC_chum_recruitment_PieterVanWill.xlsx", range="A2:BP19", trim_ws = TRUE)

# function to format harvest and totals
format_totals <- function(x, new_name) {
  col <- ncol(x) -1 
  x1 <- x %>% group_by(CU) %>% summarise(across("2018":"1954", sum, na.rm=T)) # summarise by CU, note variables are indexed by number, sensitive to numner of years
  xl <- x1 %>% pivot_longer(cols=grep("[[:digit:]]{4}", names(x1)), names_to="Year", values_to=new_name) # wide to long format
  xl
}

h1 <- format_totals(harvest_total2, new_name="harvest")
e1 <- format_totals(escape_total, new_name="escape")
str(h1)
str(e1)
eh <- merge(h1, e1, by=c("CU", "Year"), all=TRUE) # merge harvest and returns
eh$exploit_rate <- eh$harvest / (eh$escape + eh$harvest) # calculate exploitation rate
# get mean exploitation rate by CU
ehs <- eh %>% group_by(CU) %>% summarise(mean_exploit_rate = mean(exploit_rate, na.rm=TRUE))
write.csv(ehs, "DataOut/mean_exploitation_rate_by_CU_missing_1953.csv", row.names = FALSE)

# plot exploitation rate
# ggplot(eh, aes(y=exploit_rate, x=Year, group=CU)) +
#   geom_point() +
#   geom_path() +
#   facet_wrap(~CU) +
#   scale_x_discrete(breaks=seq(1960,2020,10)) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5))

# Double check that wild exploitation rates (estimated) are matching total
# There are some years where estimated wild returns is lower than escapement.
#   This is for years*CUs with CU-level infilling.  
check_ER <- Btable[!Btable$Year==1953,] # remove 1953 without harvest data
check_ER <- check_ER[check_ER$CUinfill==FALSE, ] # keep only years*CUs without CU level infilling
check_ER$exploit_rate <- (check_ER$Return - check_ER$Escape) / check_ER$Return

ers <- check_ER %>% group_by(CU) %>% summarise(mean_ER = mean(exploit_rate, na.rm=TRUE))
# confirm yes, exploitation rates are nearly the same for total and wild fish
# (based on same data, but wanted to check)
# ggplot(check_ER, aes(y=exploit_rate, x=Year, group=CU)) +
#   geom_point() +
#   geom_path() +
#   facet_wrap(~CU) +
#   scale_x_discrete(breaks=seq(1960,2020,10)) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5))




# ----------------------------------------------------#
# Explore data with figures - LW ------------
# ----------------------------------------------------#

# format data
yrcols <- grep("[[:digit:]]{4}", names(rawdat)) # get position of yr columns
rdl <- rawdat %>% pivot_longer( cols= yrcols, names_to="year", values_to="escapement") # wide to long format for plotting
#check_names <- rdl %>% group_by(CU_Name, GroupName, GU_Name) %>% summarise(n=n()) # check correspondence of CU and GroupName, GU_Name

rdl <- rdl[rdl$Source=="Wild", ] # remove non-wild fish
rdl <- rdl[!(rdl$NME %in% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER")), ] # FLAG - why? remove three rivers removed from analysis below
rdl <- rdl[ rdl$SummerRun==FALSE, ] # remove summer run fish

#-----------#
# Look at raw escapement data - LW
#-----------#

# Appears that GroupName and GU_Name are not nested within CUs, but different categories that overlap irregularly

# Check "QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER"
rdl[rdl$NME %in% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER"),] %>%
  ggplot( aes(y=escapement, x=year, colour=Source, colour=NME, group=Source)) +
  geom_point() +
  geom_path() +
  scale_y_log10() +
  facet_wrap(~NME) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

#------------#
# Look at by-stream infilled data - LW
#------------#
# Plot observed and infilled escapement for each CU, one panel per stream
infill_by_stream <- wild_infill_by_stream_list[[1]] # get by-stream infilled data
CUs <- unique(infill_by_stream$CU_Name)
# add column with data type
infill_by_stream$data_type <- if_else(!is.na(infill_by_stream$Escape), "observed",
                                            ifelse (infill_by_stream$GroupEsc=="NaN" , "infilled by CU",
                                                    "infilled by stream"))
infill_by_stream1 <- infill_by_stream[!infill_by_stream$data_type=="infilled by CU", ]
make_stream_figs <- function(CU) {
  infill_by_stream1[infill_by_stream1$CU_Name==CU,] %>% 
    ggplot() +
    geom_path(aes(x=Year, y=SiteEsc, group=NME), colour="gray")+
    geom_point(aes(x=Year, y=SiteEsc, colour=data_type)) +
    facet_wrap(~NME, scales="free_y") +
    scale_x_discrete(breaks=seq(1960,2010,10)) +
    ggtitle(CU) +
    geom_hline(aes(yintercept=0))+
    coord_cartesian(expand=FALSE, clip="off")+
    scale_colour_manual(values=c("observed"= "black", "infilled by stream" = "dodgerblue"))+
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5),
          strip.background = element_blank(),
          strip.text = element_text(hjust=0, size=10))
}
fig_list <- as.list(CUs)
fig_list <- lapply(fig_list, make_stream_figs)
fig_list[[1]] # check
pdf("Figures/fig_by-stream_infill_escapement_streams.pdf", width=18, height=12, pointsize=8)
for (i in 1:length(fig_list)){
  print(fig_list[[i]])
}
dev.off()

# Plot correlations between actual stream escapements by CU. 
# Infilling by stream assumes escapements in streams within a CU have correlation. 
# change campbell river duplicate
rdl$NME[which(rdl$NME=="CAMPBELL RIVER" & rdl$GU_Name=="5 - Fraser")] <- "LITTLE CAMPBELL RIVER"

corCU <- function(CU) {
  df1 <- rdl %>% filter(CU_Name==CU) %>% select(NME, year, escapement) %>%
    pivot_wider(names_from=NME, values_from=escapement)
  cormat <- cor(df1[,-1], use="pairwise.complete.obs")
  cormat1 <- cormat[upper.tri(cormat, diag=FALSE)]
  hist(cormat1, main=CU, xlim=c(-1,1), xlab="Correlations between observed stream escapements")
  #PerformanceAnalytics::chart.Correlation(cor(df1[,-1], use="pairwise.complete.obs") )
}
png("Figures/fig_observed_spawners_by_stream_correlations.png", width=10, height=6, units="in", res=300)
layout(matrix(seq_along(CUs), ncol=3, byrow = TRUE))
for ( i in seq_along(CUs)) {
  corCU(CUs[i])
}
dev.off()

# ----------#
# Plot infilling by CU, to compare actual, infilled by stream, and infilling by CU - LW ------
# ----------#
infill_by_stream_no_fraser <- wild_infill_by_stream_list[[2]] %>% filter(!(CU_Name %in% c("8 - Lower Fraser", "9 - Fraser Canyon"))) 
png("Figures/fig_compare_actual_infill_by_stream_and_CU.png", height=7, width=7, res=300, units="in")
ggplot(data=wild_infill_by_CU[[1]], aes(y= SiteEsc, x=Year )) + # infilled by CU
  geom_point( colour="red") + # infill by CU is in red
  geom_point(data = infill_by_stream_no_fraser, aes(y=GroupEsc, x=Year), colour="dodgerblue") + # infill by stream is blue
  geom_point(data = infill_by_stream_no_fraser, aes(y=SumRawEsc, x=Year)) + # raw escapement data
  geom_path(data = infill_by_stream_no_fraser, aes(y=SumRawEsc, x=Year, group=CU_Name)) +
  geom_hline(aes(yintercept=0))+
  ylab("Escapement") +
  facet_wrap(~CU_Name, scales="free_y", ncol=2) +
  coord_cartesian(expand=FALSE, clip="off") +
  #scale_y_log10() +
  theme_classic() +
  scale_x_discrete(breaks=seq(1960,2010,10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        axis.line.x = element_line(colour=NULL, size=0),
        strip.background = element_blank())
dev.off()

# Plot all infilled escapement with 25% benchmark
benchmarks <- cdat%>% group_by(CU_Name) %>% summarise(benchmark_25 = quantile(SiteEsc, 0.25, na.rm=TRUE))
#png("Figures/fig_escapement_infilled_w_25_benchmark.png", height=5, width=10, res=300, units="in")
ggplot(data=wild_infill_by_CU[[1]], aes(y= SiteEsc, x=Year )) + # infilled by CU
  geom_point() + # infill by CU is in red
  geom_path(aes(group=CU_Name)) +
  geom_hline(aes(yintercept=0))+
  geom_hline(data=benchmarks, aes(yintercept=benchmark_25),alpha=0.7, colour="dodgerblue", linetype=1, size=1.1) +
  ylab("Escapement") +
  facet_wrap(~CU_Name, scales="free_y") +
  coord_cartesian(expand=FALSE, clip="off") +
  #scale_y_log10() +
  theme_classic() +
  scale_x_discrete(breaks=seq(1960,2010,10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        axis.line.x = element_line(colour=NULL, size=0),
        strip.background = element_blank())
#dev.off()

# Plot escapement and R/S time series on same x axis for each CU --------
# Merge escapement and recruitment data by CU and year
cdat <- merge(wild_infill_by_CU[[1]], Btable, by.x=c("CU_Name", "Year"), by.y=c("CU", "Year"), all=TRUE)
str(cdat)
cdat$Year <- as.numeric(cdat$Year)
cdat$RS <- cdat$Recruit/cdat$SiteEsc

# cdat$RS_pos <- ifelse(cdat$RS >1, 1, 2) # binary colours
# cols <- c("dodgerblue", "firebrick")
ncols <- 20
# get colours
colors <- paletteer::paletteer_c( palette = "scico::roma", n = ncols)
#colors <- paletteer::paletteer_c( palette = "pals::warmcool", n = ncols)

# Transform the numeric variable in bins
rank <- as.factor( as.numeric( cut(log(cdat$RS), ncols)))

options(scipen = 100000)

png("Figures/fig_escapement_RS.png", height=10, width=18, res=200, pointsize=20, units="in")
layout(mat=matrix(1:16, nrow=4, byrow=FALSE))
  yrs <- range(cdat$Year)
  CUs <- unique(cdat$CU_Name)
for(i in 1:length(unique(CUs))) {
  dat <- cdat[cdat$CU_Name==CUs[i], ] 
  par(mar=c(1,4,1,1)+0.1, bty="l", las=1)
  plot(y=dat$SiteEsc/1000, x=dat$Year, type="l", xlab="year", ylab="Escapement (thousands)", 
       main=CUs[i], xlim =c(yrs[1], yrs[2]), ylim=c(0,max(dat$SiteEsc, na.rm=TRUE)/1000*1.04), yaxs="i")
  grid(ny=0)
  # FLAG: upper Knight has returns greater than spawners in some years, impossible. Must be error.
  # segments(y0=dat$SiteEsc[-1]/1000, y1=dat$Return[-1]/1000, x0=dat$Year[-1],x1=dat$Year[-1], col=adjustcolor("red", alpha=0.6), lwd=4)   # add harvest, remove first year, as there is a problem with returns
  points(y=dat$SiteEsc/1000, x=dat$Year, xlab="year", pch=ifelse(dat$Year %% 2 == 0, 16,1))
  points(y=dat$SiteEsc[dat$CUinfill==TRUE]/1000, x=dat$Year[dat$CUinfill==TRUE], col="red", xlab="year", 
         pch=ifelse(dat$Year[dat$CUinfill==TRUE] %% 2 == 0, 16,1 ))
  par(mar=c(2,4,1,1)+0.1, bty="l", las=1)
  plot(y=dat$RS, x=dat$Year, log="y", type="l", ylab="Recruits/Spawner", xlim =c(yrs[1], yrs[2]), 
      xlab="Year", ylim=c(min(cdat$RS, na.rm=TRUE), max(cdat$RS, na.rm=TRUE)), yaxt="n")
  axis(side=2, at=c(0.01, 0.1, 1, 10, 100), labels=c(0.01, 0.1, 1, 10, 100))
  points(y=dat$RS, x=dat$Year, pch=16, col= colors[ rank[cdat$CU_Name==CUs[i]] ])
    grid(ny=0)
    abline(h=1, lty=3)
    abline(h=median(dat$RS, na.rm=TRUE),lty=2, col="orange")
    text(x=2020, y= median(dat$RS, na.rm=TRUE), col="orange", label=round(median(dat$RS, na.rm=TRUE), 1), adj=c(1,-0.5), )
    }
dev.off()

# Plot abundance on y axis and productivity on x axis
ggplot(cdat[!is.na(cdat$Recruit), ], aes(y=SiteEsc, x=RS, fill=Year)) +
  geom_vline(aes(xintercept=1), colour="gray", linetype=2) +
  geom_hline(aes(yintercept=0)) +
  geom_segment(aes(x=RS, xend=RS, y=SiteEsc, yend=Return), colour="coral", alpha=0.4, size=1.5) +
  geom_path(aes(colour=Year)) +
  geom_point(shape=21, size=4, alpha=0.8) +
  scale_x_log10() +
  scale_fill_viridis_c() +
  scale_colour_viridis_c() +
  coord_cartesian(expand=FALSE, clip="off") +
  ylab("Spawners") +
  xlab("Recruits/Spawner") +
  facet_wrap(~CU_Name, scales="free_y") +
  theme_classic()

ggplot(cdat, aes(x=SiteEsc, y=RS, fill=Year)) +
  geom_hline(aes(yintercept=1), colour="gray", linetype=2) +
  #geom_segment(aes(y=RS, yend=RS, x=SiteEsc, xend=Return), colour="coral", alpha=0.4, size=1.5) +
  geom_path(aes(colour=Year)) +
  geom_point(shape=21, size=4, alpha=0.8) +
  scale_y_log10() +
  scale_fill_viridis_c() +
  scale_colour_viridis_c() +
  xlab("Spawners") +
  ylab("Recruits/Spawner") +
  facet_wrap(~CU_Name, scales="free") +
  theme_classic()


# Look at density of spawners
png(filename="Figures/fig_spawner_distribution.png", width=8, height=4,units="in", res=300)
ggplot(cdat, aes(x=SiteEsc, colour=CU_Name, fill=CU_Name)) +
  geom_point(aes(y=0, x=SiteEsc), shape=108, colour="black", size=2) +
  geom_density(alpha=0.5) + 
  scale_x_log10( breaks= c(10^(1:10))) +
  xlab(bquote("log"[10]~"(spawners)")) +
  ylab("Density") +
  theme_classic()
dev.off()

# -----------#
# Compare sum of max escapement summed across streams to infilling data - LW 
# -----------#

head(Btable)

head(rdl)
CUs <- unique(rdl$CU_Name)

sum_escp <- rdl %>% group_by(NME, CU_Name) %>% summarise(max_escape = max(escapement, na.rm=TRUE), min_escape = min(escapement, na.rm=TRUE)) %>% 
  group_by(CU_Name) %>% summarise(sum_max_escape= sum(max_escape, na.rm=TRUE), sum_min_escape=sum(min_escape, na.rm=TRUE))

Btable_max_escp <- merge(Btable, sum_escp, by.x="CU", by.y="CU_Name", all.x=TRUE)
head(Btable_max_escp)

Btable_max_escp$RS <- Btable_max_escp$Recruit/ Btable_max_escp$Escape
Btable_max_escp$logRS <-log(Btable_max_escp$RS)

# Plot recruits vs. spawners using infilled data, add vertical lines at 
# sum(max observed escapement for each stream) and sum(min observed escapement for each stream)
ggplot(Btable_max_escp, aes(x=Escape, y=Recruit,  fill=logRS)) + 
  geom_point(aes(size=logRS), shape=21, alpha=0.5) + 
  scale_fill_viridis_c() +
  geom_vline( aes(xintercept=sum_max_escape), colour="dodgerblue") +
  geom_vline(aes(xintercept=sum_min_escape), colour="firebrick") +
  geom_text(aes(x =sum_min_escape, y=0, label="Min obs spawners", angle=90, hjust=0)) + 
  geom_text(aes(x =sum_max_escape, y=0, label="Max obs spawners", angle=90, hjust=0)) + 
  geom_abline(aes(intercept=0, slope=1), colour="gray", lty=2) +
  facet_wrap(~CU, scales="free")

ggplot(Btable_max_escp, aes(x=Year, y=RS)) + 
#ggplot(Btable_max_escp, aes(x=Year, y=RS, fill=logRS)) + 
  geom_path() +
  geom_point(shape=21) + 
  #geom_point(shape=21, alpha=0.5, aes( size=logRS)) + 
  scale_fill_viridis_c() +  geom_hline(aes(yintercept=1)) +
  coord_cartesian(ylim=c(0,5)) +
  #scale_y_log10() +
  facet_wrap(~CU ) 

med_RS <- Btable_max_escp %>% group_by(CU) %>% summarise(median_RS = median(RS, na.rm=TRUE))
ggplot(Btable_max_escp, aes(x=RS)) +
  geom_density(colour="dodgerblue") + 
  geom_vline(aes(xintercept=1), lty=2) +
  scale_x_log10() +
  geom_vline(data=med_RS, aes(xintercept=median_RS), colour="dodgerblue") +
  facet_wrap(~CU) 


# Look at difference between recruits and escapement (Fisheries harvest)
png("Figures/fig_escapement_returns_harvest.png", height=8, width=12, units="in", res=300)
ggplot(Btable, aes(x=Year, y=Return)) +
  geom_point() +
  geom_path() +
  geom_linerange(aes(x=Year, ymin=Escape, ymax=Return), colour="red", lwd=1.3) +
  geom_ribbon(aes(x=Year, ymin=0, ymax=Escape), fill="dodgerblue", colour="dodgerblue") +
  #geom_text(aes(x=Year, y=Return+5000, label=round((Return-Escape)/Return*100,0)), size=4) +
  facet_wrap(~CU, scales="free_y")
dev.off()

# How would an alternative to Sgen that is spawners required to get to 
# a benchmark like 50% of max observed spawners compare to Sgen?
# incorporates alpha. 
# read in alpha and beta estimates from LRP integrated model
ests <- read.csv("DataIn/ricker_est_to_compare_with_beta.csv")
head(ests)
sum_escp <- sum_escp[1:7, ]

CUs <- unique(ests$CU_name)


curve( ests$A[1] * x * exp( - ests$B[1] * x), xlim=c(0,100000))
abline(h=sum_escp$sum_max_escape[1]*0.2 )
abline(h=sum_escp$sum_max_escape[1]*0.2 )
abline(v= ests$Sgen[1], col="grey")


function()

layout(mat= matrix(1:9, byrow=TRUE))
for (i in seq_along(CUs)) {
  plot(x=0,y=0,xlim=c(0, sum_escp$sum_max_escape[sum_escp$CU_Name==i]*1.1), ylim=c(0,sum_escp$sum_max_escape[sum_escp$CU_Name==i]*6))
  abline(v=sum_escp$sum_max_escape[sum_escp$CU_Name==i], col="dodgerblue")
  abline(v=sum_escp$sum_min_escape[sum_escp$CU_Name==i], col="firebrick")
  abline(v= ests$Sgen[sum_escp$CU_Name==i], col="grey")
  curve( ests$A[sum_escp$CU_Name==i] * x * exp( - ests$B[sum_escp$CU_Name==i] * x))
  main(i) 
}

# ----------#
# Look at age composition - LW
# ----------#
head(ACdat)

ACdat %>% pivot_longer(cols=grep("Age", names(ACdat)), names_to="Age", values_to="proportion") %>%
  ggplot(., aes(y=proportion, x=Year, fill=Age)) +
  geom_col()


# OBS


# Plot raw escapement, normal vs. summer run, 1 page per CU
# options(scipen = 100000)
# CUs <- unique(rdl$CU_Name)
# make_CU_figs <- function(CU) {
#   rdl[rdl$CU_Name==CU,] %>% ggplot( aes(y=escapement, x=year, colour=SummerRun, group=NME)) +
#     geom_point() +
#     geom_path() +
#     scale_y_log10() +
#     #facet_grid(NME~.) +
#     scale_colour_manual(values=c("black", "red")) +
#     ggtitle(CU) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle=90, vjust=0.5))
# }
# fig_list <- as.list(CUs)
# fig_list <- lapply(fig_list, make_CU_figs)
# fig_list[[1]] # check
# pdf("Figures/fig_raw_escapement_by_CU.pdf", width=7, height=4, pointsize=8)
# for (i in 1:length(fig_list)){
#   print(fig_list[[i]])
# }
# dev.off()

# Plot actual vs. infilled data by stream, one page per CU

# CUs <- unique(infill_by_stream$CU_Name)
# make_CU_figs <- function(CU) {
#   infill_by_stream[infill_by_stream$CU_Name==CU,] %>% ggplot( aes(y=SiteEsc, x=Year, group=NME)) +
#     geom_point(colour="dodgerblue", shape=1) +
#     #geom_path() +
#     #geom_point(aes(y=ContrEsc, x=Year), colour="dodgerblue") +
#     geom_point(aes(y=Escape, x=Year), colour="black", shape=1) +
#     ggtitle(CU) +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle=90, vjust=0.5))
# }
# fig_list <- as.list(CUs)
# fig_list <- lapply(fig_list, make_CU_figs)
# fig_list[[1]] # check
# pdf("Figures/fig_by-CU_infill_escapement.pdf", width=7, height=4, pointsize=8)
# for (i in 1:length(fig_list)){
#   print(fig_list[[i]])
# }
# dev.off()




