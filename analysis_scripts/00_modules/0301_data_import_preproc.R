

## data files 
if(TRUE) {
  

  deg.filelist = list()
  
  for(my.cohort in cohorts) {
    print(my.cohort)
    
    for(my.tissue in tissues) {
      print(my.tissue)
      
      deg.filelist[[my.cohort]][[my.tissue]] = ''
      
      my.file  = NULL 
      my.file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.wilcox.csv')
      
      if(file.exists(my.file)) {
        deg.filelist[[my.cohort]][[my.tissue]] = my.file
      }
      
    }
  }
  
  print(deg.filelist)
  
  saveRDS(deg.filelist,
          file = paste0(output,'.deg.filelist.rds'))
  
  ## --------------------------------------------
  
  data.filelist = list()
  
  for(my.cohort in cohorts) {
    print(my.cohort)
    
    for(my.tissue in tissues) {
      print(my.tissue)
      
      data.filelist[[my.cohort]][[my.tissue]] = ''
      
      my.file  = NULL 
      my.file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.csv')
      
      if(file.exists(my.file)) {
        data.filelist[[my.cohort]][[my.tissue]] = my.file
      }
      
    }
  }
  
  print(data.filelist)
  
  saveRDS(data.filelist,
          file = paste0(output,'.data.filelist.rds'))
  
}

## --------------------------------------------

## import data 
if(TRUE) {
  
  deg = list()
  
  for(my.cohort in cohorts) {
    print(my.cohort)
    
    for(my.tissue in tissues) {
      print(my.tissue)
      
      deg[[my.cohort]][[my.tissue]] = ''
      
      my.file  = NULL 
      my.file = deg.filelist[[my.cohort]][[my.tissue]]
      
      if(file.exists(my.file)) {
        deg[[my.cohort]][[my.tissue]] = 
          read.csv(my.file, row.names = 1)
      }
      
    }
  }
  
  head(deg$ipinivo$all)
  
  saveRDS(deg,
          file = paste0(output,'.deg.rds'))
  
  ## --------------------------------------------
  
  data = list()
  
  for(my.cohort in cohorts) {
    print(my.cohort)
    
    for(my.tissue in tissues) {
      print(my.tissue)
      
      data[[my.cohort]][[my.tissue]] = ''
      
      my.file  = NULL 
      my.file = data.filelist[[my.cohort]][[my.tissue]]
      
      if(file.exists(my.file)) {
        data[[my.cohort]][[my.tissue]] = 
          read.csv(my.file, row.names = 1)
      }
      
    }
  }
  
  head(data$ipinivo$all)
  
  saveRDS(data,
          file = paste0(output,'.data.rds'))
}