
set.seed(seed)

print(dim(data))

if(TRUE & model.type == 'xgboost') {
  
  set.seed(seed)
  
  if(my.tissue %in% c('all','panorgan')) {
    
    cv_opts = trainControl(method="cv", 
                           # seed = as.list(rep(seed,cv.fold+1)),  
                           summaryFunction=twoClassSummary,
                           savePredictions = T,
                           number= cv.fold,
                           returnResamp = "all", # save losses across all models
                           classProbs = TRUE, # set to TRUE for AUC to be computed
                           allowParallel = TRUE,
                           verboseIter = TRUE ## print training log on screen
    )
  } else if(my.tissue %in% c(as.vector(unlist(ct.tissues)),
                             'Brain')) {
    
    set.seed(seed)
    
    loocv.number = cv.fold 
    
    cv_opts = trainControl(method="LOOCV", 
                           # seed = as.list(rep(seed,cv.fold+1)), 
                           summaryFunction=twoClassSummary,
                           savePredictions = T,
                           # number= cv.fold * 4,
                           number= loocv.number, ## too small N, customize it
                           returnResamp = "all", 
                           classProbs = TRUE, 
                           allowParallel = TRUE,
                           verboseIter = TRUE 
    )
  }
  
  ## ------------------------------------------------------
  
  set.seed(seed)
  
  xgb.grid.test <- expand.grid(
    nrounds = 100,
    eta = c(0.2),
    max_depth = c(6),
    gamma = c(1), 
    subsample = c(0.5, 0.75, 1),
    min_child_weight = c(1), 
    colsample_bytree = c(1)
  )
  
  xgb.grid <- expand.grid(
    nrounds = 1000,
    eta = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
    max_depth = c(2, 3, 4, 5, 6, 7, 8, 9, 10),
    gamma = c(1, 2, 3), 
    subsample = c(0.5, 0.75, 1),
    min_child_weight = c(1, 2, 3), 
    colsample_bytree = c(1)
    
  )
  
  set.seed(seed)
  xgb_tune <-train(Class ~ .,
                   data=data,
                   method="xgbTree",
                   trControl=cv_opts,
                   tuneGrid=xgb.grid,
                   verbose=T,
                   metric="ROC",
                   nthread = threads
  )
  
  
  print(xgb_tune)
  
  saveRDS(xgb_tune,
          file = file.path(out.dir.2,
                           paste0(output,'.sm',
                                  nrow(data),'_feat',ncol(data),
                                  '.xgb_tune.rds')))
  
  ## ------------------------------------------------------
  
  # scatter plot of the AUC against max_depth and eta
  p1 = NULL 
  p1 = ggplot(xgb_tune$results, 
              aes(x = as.factor(eta), 
                  y = max_depth, 
                  size = ROC, 
                  color = ROC)) + 
    geom_point() + 
    theme_bw() + 
    scale_size_continuous(guide = "none")
  
  pdf(file = file.path(out.dir.2,
                       paste0(output,'.sm',
                              nrow(data),'_feat',ncol(data),
                              '.xgb_tune.rocTune.pdf')), 
      width = 10, height = 10)
  print(p1)
  dev.off()
  
  ## ------------------------------------------------------
  
  # library(plotROC)
  
  # Select a parameter setting
  print(xgb_tune$bestTune)
  xgb_tune.bestTune = xgb_tune$bestTune
  
  selectedIndices = NULL 
  selectedIndices = xgb_tune$pred
  selectedIndices = 
    which(
      selectedIndices$nrounds==xgb_tune.bestTune$nrounds
      & selectedIndices$max_depth==xgb_tune.bestTune$max_depth
      & selectedIndices$eta==xgb_tune.bestTune$eta
      & selectedIndices$gamma==xgb_tune.bestTune$gamma
      & selectedIndices$colsample_bytree==xgb_tune.bestTune$colsample_bytree
      & selectedIndices$min_child_weight==xgb_tune.bestTune$min_child_weight
      & selectedIndices$subsample==xgb_tune.bestTune$subsample
      
    )
  
  g <- ggplot(xgb_tune$pred[selectedIndices, ], 
              # aes(m=bDC, d=factor(obs, levels = c( "aPD", "bDC")))) + 
              aes(m=bDC, 
                  d=as.numeric(factor(obs, levels = c( "aPD", "bDC"))) - 1
              )
  ) + 
    geom_roc(n.cuts=0) + 
    coord_equal() +
    style_roc()
  
  g = g + annotate("text", x=0.75, y=0.25, 
                   label=paste("AUC =", round((calc_auc(g))$AUC, 4)))
  # Warning message:
  # In verify_d(data$d) : D not labeled 0/1, assuming NR = 0 and R = 1!
  
  
  pdf(file.path(out.dir.2,
                paste0(output,'.sm',
                       nrow(data),'_feat',ncol(data),
                       '.xgb_tune.roc.pdf')), 
      width = 8, height = 6)
  print(g)
  dev.off()
  
  ## ------------------------------------------------------
  
  ## saving cross-validation pred auc 
  if(TRUE) {
    
    library(MLeval)
    # library(caret)
    
    # Open a connection to a file for writing logs
    log_file = NULL 
    log_file <- file(file.path(out.dir.2,
                               paste0(output,'.sm',
                                      nrow(data),'_feat',ncol(data),
                                      '.xgb_tune.evalm.log')), open = "w")
    
    # Redirect standard output and standard error to the file
    sink(log_file, type = "output")
    sink(log_file, type = "message")
    
    
    ## run MLeval
    res = NULL 
    res <- evalm(xgb_tune,
                 positive = 'bDC') ## this is cv pred 
    
    # Close the connection to the file
    sink(type = "output")
    sink(type = "message")
    
    # Close the file
    close(log_file)
    
    ## make sure you run this ... otherwise you will no longer see prints on the screen lol
    # Unsink output
    sink(type = "output", append = FALSE)
    
    # Unsink messages
    sink(type = "message", append = FALSE)
    
    
    ## List containing: 1) A ggplot2 ROC curve object for printing 2) A ggplot2 PROC object for printing 3) A ggplot2 PRG curve for printing 4) Optimised results according to defined metric 5) P cut-off of 0.5 standard results
    
    saveRDS(res, 
            file = file.path(out.dir.2,
                             paste0(output,'.sm',
                                    nrow(data),'_feat',ncol(data),
                                    '.xgb_tune.evalm.rds')))
    
    pdf(file.path(out.dir.2,
                  paste0(output,'.sm',
                         nrow(data),'_feat',ncol(data),
                         '.xgb_tune.evalm.roc.pdf')), 
        width = 8, height = 6)
    ## get ROC
    print(res$roc)
    dev.off()
    
    pdf(file.path(out.dir.2,
                  paste0(output,'.sm',
                         nrow(data),'_feat',ncol(data),
                         '.xgb_tune.evalm.cc.pdf')), 
        width = 8, height = 6)
    ## get calibration curve
    print(res$cc)
    dev.off()
    
    pdf(file.path(out.dir.2,
                  paste0(output,'.sm',
                         nrow(data),'_feat',ncol(data),
                         '.xgb_tune.evalm.prg.pdf')), 
        width = 8, height = 6)
    ## get precision recall gain curve
    print(res$prg)
    dev.off()
    
    ## save performance metrics ...
    rm(a,b)
    a = data.frame(t(data.frame(res$optres)))
    b = data.frame(t(data.frame(res$stdres)))
    
    write.csv(a,
              file = file.path(out.dir.2,
                               paste0(output,'.sm',
                                      nrow(data),'_feat',ncol(data),
                                      '.xgb_tune.evalm.optres.csv')))
    write.csv(b,
              file = file.path(out.dir.2,
                               paste0(output,'.sm',
                                      nrow(data),'_feat',ncol(data),
                                      '.xgb_tune.evalm.stdres.csv')))
    
    
    
  }
  
  ## ------------------------------------------------------
  
  ## saving cross validation training auc 
  if(TRUE 
     # & split.flag == 0 ## print for all training models ....
  ) {
    
    selectedIndices = NULL 
    selectedIndices = xgb_tune$results ## this is cv training  
    selectedIndices = 
      which(
        selectedIndices$nrounds==xgb_tune.bestTune$nrounds
        & selectedIndices$max_depth==xgb_tune.bestTune$max_depth
        & selectedIndices$eta==xgb_tune.bestTune$eta
        & selectedIndices$gamma==xgb_tune.bestTune$gamma
        & selectedIndices$colsample_bytree==xgb_tune.bestTune$colsample_bytree
        & selectedIndices$min_child_weight==xgb_tune.bestTune$min_child_weight
        & selectedIndices$subsample==xgb_tune.bestTune$subsample
        
      )
    
    xgb_tune$results[selectedIndices,]
    
    write.csv( xgb_tune$results[selectedIndices,],
               file = file.path(out.dir.2,
                                paste0(output,'.sm',
                                       nrow(data),'_feat',ncol(data),
                                       '.train_rf.auc.csv')))
    
  }
  
  ## ------------------------------------------------------
  
  xgb_tune_varimp = varImp(xgb_tune, scale=F)
  print(warnings()) ## examine the warning msgs!!!!
  ## In verify_d(data$d) : D not labeled 0/1, assuming NR = 0 and R = 1!
  
  p1 = plot(xgb_tune_varimp)
  
  xgb_tune_varimp = data.frame(xgb_tune_varimp$importance)
  
  write.csv(xgb_tune_varimp,
            file = file.path(out.dir.2,
                             paste0(output,'.sm',
                                    nrow(data),'_feat',ncol(data),
                                    '.xgb_tune.varimp.1.csv')))
  
  
  pdf(file.path(out.dir.2,
                paste0(output,'.sm',
                       nrow(data),'_feat',ncol(data),
                       '.xgb_tune.varImp.1.pdf')), 
      width = 5, height = 5)
  print(p1)
  dev.off()
  
  ## ------------------------------------------------------
  
  xgb_tune_varimp = varImp(xgb_tune, scale=T)
  print(warnings()) ## examine the warning msgs!!!!
  
  p2 = plot(xgb_tune_varimp)
  
  xgb_tune_varimp = data.frame(xgb_tune_varimp$importance)
  
  write.csv(xgb_tune_varimp,
            file =file.path(out.dir.2,
                            paste0(output,'.sm',
                                   nrow(data),'_feat',ncol(data),
                                   '.xgb_tune.varimp.2.csv')))
  
  
  pdf(file.path(out.dir.2,
                paste0(output,'.sm',
                       nrow(data),'_feat',ncol(data),
                       '.xgb_tune.varImp.2.pdf')), 
      width = 5, height = 5)
  print(p2)
  dev.off()
  
  xgb_tune_varimp = xgb_tune_varimp[rev(order(xgb_tune_varimp$Overall)),,drop=F]
  
  # print(xgb_tune_varimp)
  
  ## ---------------------------------------------------------------
  
  if(TRUE & split.flag ==1) {
    
    set.seed(seed)
    
    data_test$Class = factor(data_test$Class, levels = c('aPD','bDC'))
    data_test_st = data_test
    
    
    ## ------------------------------------------------------
    
    preds_rf = predict(xgb_tune, data_test_st[,-1])
    preds_rf2 = predict(xgb_tune, data_test_st[,-1],
                        "prob") ## print class probablit
    
    ## ------------------------------------------------------
    
    # library(pROC)

    preds_rf_roc = roc(data_test$Class,
                       predict(xgb_tune,
                               data_test_st[,-1],
                               type = "prob")[,1],
                       levels = levels(data_test_st$Class))
    
    saveRDS(preds_rf_roc,
            file.path(out.dir.2,
                      paste0(output,'.sm',
                             nrow(data),'_feat',ncol(data),
                             '.preds_rf.roc.rds')))
    
    pdf(file.path(out.dir.2,
                  paste0(output,'.sm',
                         nrow(data),'_feat',ncol(data),
                         '.preds_rf.roc.pdf')), 
        width = 8, height = 6)
    pROC::plot.roc(preds_rf_roc, 
                   # print.thres = c(.5),
                   type = "S",
                   # print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
                   print.thres.cex = .8,
                   print.auc = T,
                   legacy.axes = TRUE) 
    
    dev.off() 
    
    ## ------------------------------------------------------
    
    conf_rf2 = confusionMatrix(preds_rf, data_test[,outcome_col], 
                               positive='bDC')
    conf_rf2_table = data.frame(conf_rf2$table)
    rm(x)
    x = data.frame(conf_rf2.byClass = preds_rf_roc$auc[1])
    row.names(x) = 'AUC'
    conf_rf2_byclass = rbind(data.frame(conf_rf2$byClass),
                             x)
    sens = conf_rf2_byclass['Sensitivity',]
    spec = conf_rf2_byclass['Specificity',]
    
    # print(paste0(seed, ',', ratio, ', ', sens, ',', spec))
    
    write.csv(conf_rf2_byclass,
              file = file.path(out.dir.2,
                               paste0(output,'.sm',
                                      nrow(data),'_feat',ncol(data),
                                      '.preds_rf.auc.csv')))
    
  }
  
  
}