
## Table 1
if(TRUE) {
  
  ## prepare input
  if(TRUE) {
    
    data.plot = NULL 
    data.plot = clinical
    data.plot$Cohort2 = data.plot$Cohort
    data.plot$Cohort2[which(data.plot$Cohort %in% c('pd1','ipinivo'))] = 'ici'
    
    ## --------------------------------------------
    
    ## check values ....
    if(TRUE) {
      
      table(data.plot$Melanoma.subtype[data.plot$Cohort=='pd1'])
      
      sum(is.na(data.plot$Melanoma.subtype[data.plot$Cohort=='pd1']))
      
      table(data.plot$Melanoma.subtype)
      
      table(data.plot$BRAF.status[data.plot$Cohort=='pd1'])
      
      table(data.plot$BRAF.status[data.plot$Cohort=='ipinivo'])
      
      table(data.plot$BRAF.status[data.plot$Cohort=='braf'])
      
      table(data.plot$BRAF.status)
      
      table(data.plot$NRAS[data.plot$Cohort=='pd1'])
      
      table(data.plot$NRAS[data.plot$Cohort=='ipinivo'])
      
      table(data.plot$NRAS[data.plot$Cohort=='braf'])
      
      table(data.plot$NRAS)
      
      table(data.plot$NF1)
      
      table(data.plot$Line.of.immunotherapy)
      
      table(data.plot$Line.of.immunotherapy.clp)
      
      table(data.plot$Mets.per.AJCC)
      
      table(data.plot$Mets.per.AJCC.01)
      
      sum(as.numeric(as.character(data.plot$LDH.pre))[data.plot$Cohort=='pd1'] > 170,
          na.rm = T)
      sort(as.numeric(as.character(data.plot$LDH.pre))[data.plot$Cohort=='pd1'])
      sort(as.numeric(as.character(data.plot$LDH.pre))[data.plot$Cohort=='ipinivo'])
      sort(as.numeric(as.character(data.plot$LDH.pre))[data.plot$Cohort=='braf'])
      
      table(data.plot$LDH.pre)
      
      table(data.plot$LDH.pre.ULN)
      
      table(data.plot$LDH.pre.2xULN)
      
      sum(as.numeric(as.character(data.plot$NLR.recalc))[data.plot$Cohort=='ipinivo'] > 3,
          na.rm = T)
      sort(as.numeric(as.character(data.plot$NLR.recalc))[data.plot$Cohort=='ipinivo'])
      sum(is.na(as.numeric(as.character(data.plot$NLR.recalc))[data.plot$Cohort=='ipinivo'])) ## 1
      
      table(data.plot$NLR)
      
      table(data.plot$NLR.recalc.high)
      
    }
    
    ## --------------------------------------------
    
    data.plot.0 = NULL 
    data.plot.0 = data.plot
    
    col.select = NULL 
    col.select = c('Cohort', 'Cohort2',
                   'Age', 'Sex', 
                   'Melanoma.subtype',
                   'BRAF.status','NRAS','NF1',
                   'Line.of.immunotherapy.clp',
                   'Mets.per.AJCC.01',
                   'LDH.pre', 'LDH.pre.ULN', 'LDH.pre.2xULN', 
                   'NLR.recalc', 'NLR.recalc.high')
    
    data.plot = data.plot[,col.select]
    head(data.plot)
    
    ## --------------------------------------------
    
    write.csv(data.plot.0,
              file = paste0(clinical.file,'.characteristics.csv'))
    
    write.csv(data.plot,
              file = paste0(clinical.file,'.characteristics.slim.csv'))
    
  }
  
  ## --------------------------------------------
  
  ## draw tables
  if(TRUE) {
    
    library(gtsummary)
    
    table1 = NULL 
    table1 <-
      tbl_summary(
        data.plot,
        include = c(Age, Sex,
                    Melanoma.subtype,
                    BRAF.status, NRAS, NF1,
                    Line.of.immunotherapy.clp,
                    Mets.per.AJCC.01,
                    LDH.pre, LDH.pre.ULN, LDH.pre.2xULN, 
                    NLR.recalc, NLR.recalc.high),
        by = Cohort, # split table by group
        # missing = "no" # don't list missing data separately
      ) %>%
      add_n() %>% # add column with total number of non-missing observations
      # add_p() %>% # test for a difference between groups
      modify_header(label = "**Characteristic**") %>% # update the column header
      bold_labels()
    
    print(table1)
    
    table2 = NULL 
    table2 <-
      tbl_summary(
        data.plot[data.plot$Cohort2 == 'ici',],
        include = c(Age, Sex,
                    Melanoma.subtype,
                    BRAF.status, NRAS, NF1,
                    Line.of.immunotherapy.clp,
                    Mets.per.AJCC.01,
                    LDH.pre, LDH.pre.ULN, LDH.pre.2xULN, 
                    NLR.recalc, NLR.recalc.high),
        by = Cohort2, # split table by group
        # missing = "no" # don't list missing data separately
      ) %>%
      add_n() %>% # add column with total number of non-missing observations
      # add_p() %>% # test for a difference between groups
      modify_header(label = "**Characteristic**") %>% # update the column header
      bold_labels()
    
    print(table2)
    
    ## --------------------------------------------
    
    table1 %>%
      as_gt() %>%
      gt::gtsave(filename = paste0(clinical.file,'.characteristics.slim.01.pdf')) # use extensions .png, .html, .docx, .rtf, .tex, .ltx
    
    table2 %>%
      as_gt() %>%
      gt::gtsave(filename = paste0(clinical.file,'.characteristics.slim.02.pdf')) # use extensions .png, .html, .docx, .rtf, .tex, .ltx
    
  }
  
}
