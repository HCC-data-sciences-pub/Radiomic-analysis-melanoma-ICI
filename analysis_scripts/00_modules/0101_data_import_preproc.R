
## data files 
if(TRUE) {
  
  clinical.file = 'Mel.20240206.clinical.pt265.csv'
  ptlist.file = 'radiomics_paper.pts_final.20240206.xlsx'
  clinical.wOrganResponse.file = 'Mel.20240206.clinical.pt265.wOrganResponse.fixed.csv'
  
}

## --------------------------------------------

## import data
if(TRUE) {
  
  clinical = read.csv(clinical.file, row.names = 1)
  ptlist = list(
    
    ipinivo = data.frame(readxl::read_excel(ptlist.file, sheet = 'ipinivo', skip=2,
                                            na = c('N/A','NA','','na')),
                         stringsAsFactors = F),
    pd1 = data.frame(readxl::read_excel(ptlist.file, sheet = 'pd1', skip=2,
                                        na = c('N/A','NA','','n/a')),
                     stringsAsFactors = F),
    braf = data.frame(readxl::read_excel(ptlist.file, sheet = 'braf', skip=2,
                                         na = c('N/A','NA','','n/a')),
                      stringsAsFactors = F)
    
    
  ) 
  
  clinical.wOrganResponse = read.csv(clinical.wOrganResponse.file,
                                     row.names = 1)
  
}

## --------------------------------------------

## preprocess data
if(TRUE) {
  
  dim(clinical) 
  dim(ptlist)
  
  ptlist$pd1 = ptlist$pd1[!is.na(ptlist$pd1),,drop=F]
  
  ## --------------------------------------------
  
  colnames(clinical)
  table(clinical$Cohort)
  
  ## --------------------------------------------
  
  print(all.equal(sort(clinical$Patient.identifier[clinical$Cohort=='ipinivo']),
                  paste0('P',sprintf("%09d",sort(ptlist$ipinivo$Patient.identifier)))))
  ## TRUE
  print(all.equal(sort(clinical$Patient.identifier[clinical$Cohort=='pd1']),
                  paste0('P',sprintf("%09d",sort(ptlist$pd1$Patient.identifier)))))
  ## TRUE
  print(all.equal(sort(clinical$Patient.identifier[clinical$Cohort=='braf']),
                  paste0('P',sprintf("%09d",sort(ptlist$braf$Patient.identifier)))))
  ## TRUE
  
  
}
