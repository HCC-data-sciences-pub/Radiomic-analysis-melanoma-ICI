

rm(list = ls())
gc()

## ---------------------------------------------------------------

## parameters
if(TRUE) {

  my.tissue ='SoftTissue'  ## "Adrenal"  "Brain" "LN"   "Liver_L"   "Lung_L"  "SoftTissue"
  seed = 8
  run = 1 
  clinicvars.set = 'clin04' 
  my.cohort = 'ipinivo' ## ipinivo pd1 
  
  ## ---------------------------------------------------------------
  
  print(
    paste0(
      'Global flags: ',
      'my.tissue = ', my.tissue, ';',
      'seed = ', seed, ';',
      'run = ', run, ';',
      'clinicvars.set = ', clinicvars.set,';',
      'my.cohort = ', my.cohort,';'
    )
  )
}

## ---------------------------------------------------------------

## libraries
if(TRUE) {
  
  library(ggplot2)
  library(RColorBrewer)
  library(broom)
  library(survival)
  library(survminer)
  library(circlize)
  library(ggsci)
  library(scales)
  library(ggpubr)
  library(plotrix)
  library(fBasics)
  library(data.table)
  library(caret)
  library(AppliedPredictiveModeling)
  library(pROC)
  library(plotROC)
  library(xgboost)
  library(ggfortify)
  library(cluster)
  
}

## ---------------------------------------------------------------

source("00_modules/0401_data_settings.R")

## ---------------------------------------------------------------

source("00_modules/0402_data_import_preproc.R")

## ---------------------------------------------------------------

source("00_modules/0403_data_qc_filter.R")

## ---------------------------------------------------------------

source("00_modules/0404_analyze_figure3and4_S7S8S9.R")
    
## ---------------------------------------------------------------

sessionInfo()
