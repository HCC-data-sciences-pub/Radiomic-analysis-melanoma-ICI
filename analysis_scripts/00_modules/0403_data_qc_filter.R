 
if(TRUE) {
  
  ## prepare data 
  if(TRUE) {
    
    qc_flag = 1 
    predictor_list_ex = c()
    
    rm(predictor_list_full, predictor_list)
    predictor_list_full = list(
      run1 = c(imgvars),
      run7 = c(imgvars, clinicvars, molvars) 
      
    )
    
    predictor_list = predictor_list_full[[paste0('run',run)]]
    
    predictor_list = predictor_list[!predictor_list %in% predictor_list_ex]
    
    outcome_col = "Class"
    
    ## ---------------------------------------------------------------
    
    rm(data)
    data = mel[,c('Response',predictor_list)] 
    
    ## ---------------------------------------------------------------
    
    data$Response[which(data$Response %in% c('DC'))] = 'bDC'
    data$Response[which(data$Response %in% c('PD'))] = 'aPD'
    
    data$Class = NA
    data$Class[which(data$Response=='bDC')] = 1
    data$Class[which(data$Response=='aPD')] = 0
    
    data = data[,c('Class',predictor_list)]
    
    ## remove NA rows (83 rows removed!)
    dim(data) 

    
    ## ---------------------------------------------------------------
    
    ## stats 
    data.frame(table(data$Class))
    w = NULL 
    w = data.frame(is.NA = 
                     apply(data, 2, function(x) sum(is.na(x))))
    print(w[w$is.NA>0,,drop=F])

    if(nrow(w[w$is.NA>0,,drop=F]) > 0) {
      data = na.omit(data)
    }
    
    dim(data) 
    
    ## change outcome col to factor
    data$Class = as.character(data$Class)
    data$Class = gsub("1", "bDC", data$Class)
    data$Class = gsub("0", "aPD", data$Class)

    
    ## ---------------------------------------------------------------
    
    head(data)
    dim(data)
    
    table(data$Class)
    
  }
  
  ## ------------------------------------------------------
  
  # seed = 1
  set.seed(seed)
  
  ## convert categorical var to dummy vars 
  if(TRUE) {
    
    dim(data)
    
    ## ------------------------------------------------------
    
    print(dim(data)) ## 82 401
    
    ## convert categorical factors to dummy variables 
    col.select = NULL 
    for(i in 1:ncol(data)) {
      
      if(!colnames(data)[i] %in% c(outcome_col,'PID')) {
        if(class(data[,i]) == 'character') {
          col.select = c(col.select, colnames(data)[i])
        }
      }
      
    }
    
    # print(col.select)
    ## "Liver"     "PDL1"      "Smoking"   "Sex"       "ECOG"      "TumorType"
    
    ## numeric cols ....
    col.select.2 = NULL 
    for(i in 1:ncol(data)) {
      
      if(!colnames(data)[i] %in% c(outcome_col,'PID')) {
        if(class(data[,i]) %in% c('numeric','integer')) {
          col.select.2 = c(col.select.2, colnames(data)[i])
        }
      }
      
    }
    
    # print(col.select.2)
    ##  "PriorTx"  "NLR"      "Age"      "TumorSize" "Radiomic" "IL8" 
    
    ## ------------------------------------------------------
    
    if(length(col.select) > 0) {
      dmy = dummyVars( ~ ., data = data[,col.select,drop=F], fullRank = T)
      data_transformed = data.frame(predict(dmy, newdata = data))
      
      data = cbind(data[,outcome_col,drop=F],
                   data_transformed,
                   data[,col.select.2,drop=F])
      sum(duplicated(colnames(data))) ## 0
      
    }
    
    grep('^Sample',colnames(data)) 
    
    dim(data)
    
    output = output
    write.csv(data,
              file = file.path(out.dir.2,
                               paste0(output,'.sm',
                                      nrow(data),'_feat',ncol(data),
                                      '.data.csv')))
    
    ## ------------------------------------------------------
    
    ## cannot have NA values ...
    data = na.omit(data) ##  133 405
    dim(data)
    
    write.csv(data,
              file = file.path(out.dir.2,
                               paste0(output,'.sm',
                                      nrow(data),'_feat',ncol(data),
                                      '.data.csv')))
  }
  
  ## ------------------------------------------------------
  
  set.seed(seed)
  
  ## split training/test set 
  if(TRUE & split.flag == 1) {

    trainIndices = createDataPartition(data[,outcome_col], 
                                       p = ratio, list = F) 
    
    ## training and test sets
    data_train = data[trainIndices, ] 
    data_test = data[-trainIndices, ] 
    
    ## stats 
    data.frame(table(data_train$Class))
    data.frame(table(data_test$Class))
    
    data = data_train
    
    
  }
  
  dim(data)
  
  
  ## ------------------------------------------------------
  
  set.seed(seed)
  
  ## QC
  ## feature filtering (QC): remove low vars, high cor, and high colinearity 
  if(TRUE & qc_flag == 1 & my.tissue %in% c('all',
                                            'panorgan',
                                            
                                            as.vector(unlist(ct.tissues)),
                                            'Brain')) {
    
    dim(data) ## 141
    values = data[,-1,drop=F]
    samples = data[,1,drop=F]
    
    ## for testing purpose 
    if(seed == 2) {
      
      write.csv(values,
                file = file.path(out.dir.2,
                                 paste0(output,'.sm',
                                        nrow(data),'_feat',ncol(data),
                                        '.data.values.csv')))
      
      write.csv(samples,
                file = file.path(out.dir.2,
                                 paste0(output,'.sm',
                                        nrow(data),'_feat',ncol(data),
                                        '.data.samples.csv')))
    }
    
    ## -------------------------------------
    
    x = data.frame(var = apply(scale(values), 2, var))
    x = x[order(x$var),,drop=F]
    
    dim(values)
    ## 2. Zero- and Near Zero-Variance Predictors ("unique" or "unbalanced" predictor): 
    nzv = nearZeroVar(values, saveMetrics= TRUE) 
    nzv
    
    nzv.ex = row.names(nzv)[nzv$nzv==T]
    values = values[,! colnames(values) %in% nzv.ex,drop=F]
    dim(values) ## 116
    
    ## -------------------------------------
    
    ## remove high cor feat 
    if(FALSE & ncol(values) >= 2) {
      
      uniq.feat = colnames(values )
      feat.rescue = colnames(values)
      
      while(length(uniq.feat) > 0 ) {
        
        
        descrCor = cor(values[,uniq.feat], method = 'spearman')
        
        ## for testing purpose 
        if(seed == 2) {
          
          write.csv(descrCor,
                    file = file.path(out.dir.2,
                                     paste0(output,'.sm',
                                            nrow(data),'_feat',ncol(data),
                                            '.data.values.descrCor.csv')))
        }
        
        ## -------------------------------------
        
        if(TRUE) {
          
          dim(descrCor)
          
          highlyCor = list() 
          for(i in 1:nrow(descrCor)) {
            # print(i)
            
            x = NULL 
            y = NULL 
            
            x = row.names(descrCor)[i]
            y = descrCor[i,i:ncol(descrCor)]
            y = names(y)[abs(y) > cor.abs.cutoff]
            y = y[!y %in% x]
            
            ## keep this one, drop its highly correlated mates 
            highlyCor[[row.names(descrCor)[i]]] = y
            
          }
          
          
          ## -------------------------------------
          
          ## unique feats
          lowlyCor = NULL 
          ## features highly cor with unique feats
          highlyCor.2 = NULL 
          
          for(i in 1:length(highlyCor)) {
            # print(i)
            
            ## removing cor feat sequentially ...
            ## keep this feat, if this feat has not been removed yet
            
            x = NULL 
            y = NULL 
            
            x = names(highlyCor)[i]
            y = highlyCor[[i]]
            
            highlyCor.2 = c(highlyCor.2, y)
            ## exclude the already recorded unique feats ...
            highlyCor.2 = highlyCor.2[! highlyCor.2 %in% lowlyCor]
            if(! x %in% highlyCor.2) {
              lowlyCor = c(lowlyCor, x)
            }
            
            
            
          }
          
          lowlyCor = unique(sort(lowlyCor))
          
          if(length(highlyCor.2) == 0) { break }
          
          ## -------------------------------------
          
          if(TRUE) {
            
            ## feat that I kept ...
            if(TRUE) {
              
              x = NULL 
              y = NULL 
              x = values[,lowlyCor,drop=F]
              y = cor(x, method = 'spearman')
              
              y[abs(y)<= cor.abs.cutoff] = NA
              diag(y) = NA
              y
              # print(sum(!is.na(y))==0) ## 0
              ## shouild be all NAs!!
              ## TRUE :)
              
            }
            
            ## -------------------------------------
            
            ## feat that I did not keep ...
            feat.rescue = NULL 
            if(TRUE) {
              
              for(my.feat in row.names(descrCor)) {
                
                if(my.feat %in% lowlyCor) { next }
                
                # print(my.feat)
                
                x = NULL 
                y = NULL 
                
                x = descrCor[row.names(descrCor) == my.feat,lowlyCor]
                if(sum(abs(x) > cor.abs.cutoff) == 0) {
                  # print(my.feat)
                  
                  feat.rescue = c(feat.rescue, my.feat)
                }
                
                
              }
              
            }
            
            # print(length(feat.rescue))
            
            if(length(feat.rescue) == 0) { break }
            
            ## -------------------------------------
            
            intersect(lowlyCor, feat.rescue) ## should be no overlap!
            
            uniq.feat = NULL 
            uniq.feat = c(lowlyCor, feat.rescue)
            
          }
          
        }
        
        if(length(feat.rescue) == 0) { break }
        if(length(highlyCor.2) == 0) { break }
      } 
      
      ## -------------------------------------
      
      values = values[,lowlyCor,drop=F]
      
      print(dim(values)) ## 84 12
      
    }
    
    ## -------------------------------------
    
    ## 4. Linear Dependencies
    ## identify the linear dependencies between predictors
    comboInfo = findLinearCombos(values)
    comboInfo
    if(length(comboInfo$remove) > 0) {
      values = values[, -comboInfo$remove]
    }
    dim(values) ## 107
    
    ## -------------------------------------
    
    data = cbind(samples, values)
    dim(data)
    
    # output = mel.file
    write.csv(data,
              file =file.path(out.dir.2,
                              paste0(output,'.sm',
                                     nrow(data),'_feat',ncol(data),
                                     '.data_qc_flt.csv')))
  }
  
  
  
  
}