

if(TRUE) {
  
  print(run)
  
  ## ---------------------------------------------------------------
  
  mel.file = NULL 
  mel = NULL 
  
  ## import data
  if(TRUE) {
    
    
    mel.file =  radiomics.filelist[[my.cohort]][[my.tissue]]
    print(mel.file)
    
    if(length(mel.file) != 1){ print("ERROR! more than one mel.file") ; next }
    
    
    mel = read.csv(mel.file, stringsAsFactors = F, row.names = 1)
    colnames(mel)[which(colnames(mel) == 'response')] = 'Response'
    
    print(mel.file)
    print(table(mel$Response))
    
  } 
  
  
  ## ---------------------------------------------------------------
  
  mel = mel[,!colnames(mel) %in% c('Best.Response.by.RECIST',
                                   'BRR','response',
                                   'PID')]
  
  print(mel[1:3,1:9])
  
  ## ---------------------------------------------------------------
  
  rm(out.dir.2)
  out.dir.2 = file.path(out.dir,seed,paste0('run',run))
  if(!dir.exists(out.dir.2)) { dir.create(out.dir.2, recursive = T) }
  
  ## ---------------------------------------------------------------
  
  write.csv(mel,
            file = file.path(out.dir.2,
                             paste0(output,'.sm',
                                    nrow(mel),
                                    '.data.csv')))
  
  
}