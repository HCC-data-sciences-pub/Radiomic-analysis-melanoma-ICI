
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
              'Brain',
              'LN','Liver_L','Lung_L','SoftTissue',
              'all')
  
  output = 'Mel.20240206.radiomics'
  
  ## --------------------------------------------
  
  plot.colors = c('DC' = '#CCCC00','PD' = '#0000CC')
  
  ## --------------------------------------------
  
  fdr = 0.10
  rawp = 0.01 ## 0.01 0.05
  mean.diff.min = 0.001
  expr.min = 0.30 ## higher than mean in at least 10% samples overall ...
  
  use.fdr = 1
  use.rawp = 0
  if(use.fdr == 1) { use.rawp = 0 }
  
  ## --------------------------------------------
  
  folder = 'fdr0.10'
  
}

## --------------------------------------------

path = file.path('results/20240206/20240303b_archive_DE/DE_analysis',folder)
print(path)
setwd(path)

## --------------------------------------------

source("00_modules/0301_data_import_preproc.R")

## --------------------------------------------

source("00_modules/0302_analyze_figure2.R")


## --------------------------------------------

