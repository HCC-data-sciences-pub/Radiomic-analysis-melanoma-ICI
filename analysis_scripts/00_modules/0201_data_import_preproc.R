
## data files 
if(TRUE) {
  
  file_list.file = '../Mel.20240206.radiomics.file_list.20240206.xlsx'
  
}

## --------------------------------------------

## import data 
if(TRUE) {
  
  file_list = data.frame(readxl::read_excel(file_list.file, sheet = 1),
                         stringsAsFactors = F)
  
  ## --------------------------------------------
  
  file_list
  
  file_list = file_list[file_list$Run ==1,]
  
  file_list
  
  ## --------------------------------------------
  
  radiomics = list()
  
  for(my.cohort in cohorts) {
    print(my.cohort)
    
    radiomics[[my.cohort]] = list()
    
    for(my.tissue in tissues) {
      print(my.tissue)
      
      my.df = NULL 
      my.df = file_list$Filename[file_list$Cohort == my.cohort &
                                   file_list$Tissue == my.tissue]
      my.df = read.csv(file.path('../',my.df), row.names = 1)
      
      my.df = my.df[,!colnames(my.df) %in% c('BRR')]
      
      ## --------------------------------------------
      
      if(TRUE & my.tissue == 'Brain') {
        x = NULL 
        x = my.df[,-1]
        x = data.frame(t(x))
        x$feat = gsub('[.]\\S+$','',row.names(x))
        x = reshape2::melt(x)
        colnames(x)[2:3] = c('Sample','value')
        x = reshape2::dcast(Sample ~ feat,
                            value.var = 'value',
                            data = x,
                            fun.aggregate = mean,
                            na.rm=T)
        row.names(x) = gsub('[.]','-',x$Sample)
        x = x[,-1]
        
        y = NULL 
        y = my.df[,1,drop=F]
        
        print(all.equal(row.names(x),row.names(y)))
        # TRUE
        
        my.df = cbind(y,x)
        dim(my.df) ## 23 401
        
      }
      
      ## --------------------------------------------
      
      radiomics[[my.cohort]][[my.tissue]] = my.df
      
    }
    
  }
  
  ## --------------------------------------------
  
  saveRDS(radiomics,
          file = paste0(output,'.data.rds'))
  
}
