
rm(list = ls())
gc()

## --------------------------------------------

## libs 
if(TRUE) {
  
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(circlize)
  library(RColorBrewer)
  library(patchwork)
  library(limma)
  library(broom)
  library(ComplexHeatmap)
  library(ggfortify)
  
}

## --------------------------------------------

## settings 
if(TRUE) {
  
  cohorts = c('ipinivo','pd1','braf')

  tissues = c('Adrenal',
              'Brain','LN',
              'Liver_L','Lung_L','SoftTissue',
              'all')
  
  output = 'Mel.20240206.radiomics'
  
  ## --------------------------------------------
  
  plot.colors = c('DC' = '#CCCC00','PD' = '#0000CC')
  
  ## --------------------------------------------
  
  fdr = 0.10
  rawp = 0.01 
  mean.diff.min = 0.001
  expr.min = 0.30 
  
  use.fdr = 1
  use.rawp = 0
  if(use.fdr == 1) { use.rawp = 0 }
  
  
}

## --------------------------------------------

path = 'results/20240206/20240303b_archive_DE/DE_analysis'
setwd(path)

## --------------------------------------------

source("00_modules/0201_data_import_preproc.R")

## --------------------------------------------

source("00_modules/0202_analyze_figureS6.R")

## --------------------------------------------

sessionInfo()


