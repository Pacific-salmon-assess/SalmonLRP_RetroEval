# Inside South Coast Chum Case Study

Primary contacts: Luke Warkentin (Luke.Warkentin@dfo-mpo.gc.ca), Kendra Holt (Kendra.Holt@dfo-mpo.gc.ca)

## Overview

> The SCChumStudy folder contains the files that are specific to the Limit Reference Point retrospective analysis for the Inside South Coast Chum (non-Fraser) Stock Management Unit. This includes raw input data, code to prepare data and run analyses, case study-specific functions, and outputs (data and figures). 

![The seven Conservation Units that make up the Inside South Coast Chum Stock Management Unit](./Figures/fig_chum_CU_map.png)

## Contents

### DataIn  
_(input data - read only)_  
Pieter Van Will (Pieter.VanWill@dfo-mpo.gc.ca) supplied data

#### AgeComp_2018.csv  
Age composition data from Johnstone Strait fishery aggregate

#### Chum Escapement Data With Areas_2018(CleanedFeb152021).xlsx  
Raw escapement data (spawner counts) with updated 2018 tab
Visual counts, some Area Under the Curve, some fences

#### wild_ISC_chum_recruitment_PieterVanWill.xlsx   
Wild returns generated by Pieter Van Will, using infilled escapement, harvest, and percent hatchery data from infilling escapement data by total and wild.

#### LRP_compare_methods.csv
Data output from multi-dimensional CU status analysis, State of the Salmon program, SOLV-Code/TEMP-Chum-Synoptic repository (currently private repo)


### DataOut

#### wild_spawners_CU_infilled_by_site_CU.csv  
Infilled wild spawners by CU

#### wild_spawners_stream_infilled_by_site_CU.csv  
Infilled wild spawners by stream

#### SRdatWild.csv  
Time series of spawners and recruits, made by make_brood_table.R

#### streams_counted.csv  
Summary of the streams in the data set

#### summary_n_infilled_by_year.csv  
Summary of how many streams in each CU were infilled each year

#### AnnualRetrospective  
Folder containing output data .csv files from LRP retrospective analysis

#### infill_escapement_for_external_run_reconstruction  
Folder with infilling of total and wild spawners by stream, Area, and CU, used by Pieter Van Will for run reconstruction. Files made by infill_escapement_for_external_run_reconstruction.R

### Figures  
Folder containing output figures

### R
folder containg R scripts

#### runSouthCoastChum.R
Master script to run LRP Retrospective Analyses.  
Set working directory to SalmonLRP_RetroEval/SCChumStudy folder to run (or make an .Rproj file in this folder).

#### make_brood_table.R
Prepares data for retrospective analysis, including infilling of raw escapement data and creation of stock-recruit brood table. 

#### chumDataFunctions.r
Functions specifically for chum case study, including infilling

#### infilling_for_reconstruction.R
This code takes the raw escapement data and infills by stream and CU. The outputs of this code were sent to Pieter Van Will to do the run reconstruction, which uses percent wild spawners, and catch to give wild returns (**wild_ISC_chum_recruitment_PieterVanWill.xlsx**). Does not need to be run to do LRP retrospective analysis.

#### make_map.R
Make a map of the Conservation Unit areas.




