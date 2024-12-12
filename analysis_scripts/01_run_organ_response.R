
rm(list = ls())
gc()

## --------------------------------------------

## libraries 
if(TRUE) {
  
  library(ggplot2)
  library(ggpubr)
  library(ggsci)
  library(patchwork)
  library(broom)
  library(circlize)
  library(RColorBrewer)
  
  library(psych)
  library(dplyr)
  library(webr)
  library(ComplexHeatmap)
  library(ggridges)
  
}

## --------------------------------------------

## settings
if(TRUE) {
  
  plot.colors = c('PD' = '#0000CC', 'SD' = '#66B2FF', 'PR' = '#CC00CC', 'CR' = '#CC0000',
                  'NA' = '#C0C0C0', 'NE' = '#C0C0C0','DC' = '#CCCC00',
                  'ici' = '#E0E0E0', 'pd1' = '#202020','ipinivo' = '#606060',
                  'braf' = '#A0A0A0',
                  'uniformProgression' = '#00CCCC', 'mixedResponse' = '#00FF00', 
                  'uniformDiseaseControl' = '#CC6600')
  
  ## for plotting response organ heterogeneity scores
  Response.organ.het.norm.colors = colorRamp2(c(-1,0, 1), 
                                              c('#00FFFF','#FFFFFF','#FF8000'))
  

  tissues = list(
    Adrenal = 'adrenal',
    Brain = 'brain',
    Liver = 'liver',
    LN = 'nodes',
    Lung = 'lung',
    Soft.tissue = 'soft.tissue'
  )
  
}

## --------------------------------------------

path = 'results/20240206'
setwd(path)

## --------------------------------------------

source("00_modules/0101_data_import_preproc.R")

## --------------------------------------------

source("00_modules/0102_analyze_table1.R")

## --------------------------------------------

source("00_modules/0103_analyze_figure1.R")

## --------------------------------------------

source("00_modules/0104_analyze_figureS3.R")

## --------------------------------------------

source("00_modules/0105_analyze_figureS4S5.R")

## --------------------------------------------

sessionInfo()



