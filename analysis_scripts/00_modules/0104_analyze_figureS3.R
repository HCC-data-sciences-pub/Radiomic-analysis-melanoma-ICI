
## Figure S2
if(TRUE) {
  
  dim(clinical.wOrganResponse) 
  
  ## prepare input
  if(TRUE) {
    
    data.plot = NULL 
    data.plot = clinical.wOrganResponse
    data.plot$Cohort2 = data.plot$Cohort
    data.plot$Cohort2[which(data.plot$Cohort %in% c('pd1','ipinivo'))] = 'ici'
    
    ## --------------------------------------------
    
    table(data.plot$Cohort)
    
    table(data.plot$Best.Response.by.RECIST)
    
    sum(is.na(data.plot$Best.Response.by.RECIST)) 
    
    table(data.plot$Best.Response)
    
    sum(is.na(data.plot$Best.Response)) 
    
    print(all.equal(data.plot$Best.Response.by.RECIST,
                    data.plot$Best.Response))
    ## TRUE
    
    ## --------------------------------------------
    
    print(names(tissues))
    
    ## --------------------------------------------
    
    ## label mixed response vs uniform progression vs uniform regression
    data.plot$Response.organ.mixed = NA
    for(i in 1:nrow(data.plot)) {
      print(i)
      
      if(is.na(data.plot$Best.Response[i])) { next } 
      
      my.response = NULL
      my.response = data.plot$Best.Response[i]
      my.df = NULL 
      my.df = data.plot[i,names(tissues)]
      my.response.clp = NULL 
      if(my.response %in% c('PD')) {
        my.response.clp = c('PD')
      } else if(my.response %in% c('SD','PR','CR')) {
        my.response.clp = c('SD','PR','CR')
      }
      
      my.stats = list(
        consistent = 0,
        notconsistent = 0,
        missing = 0
      )
      my.stats = list(
        consistent = sum(my.df %in% my.response.clp),
        notconsistent = sum(! my.df %in% my.response.clp & my.df != "NE"),
        missing = sum(my.df == "NE")
      )
      
      if(sum(my.stats$consistent,
             my.stats$notconsistent,
             my.stats$missing) != length(tissues)) {
        print("organ response counts not correct!")
        
      }
      
      if(my.stats$notconsistent > 0) {
        data.plot$Response.organ.mixed[i] = 'mixedResponse'
      } else if(my.stats$notconsistent == 0 & my.stats$consistent > 0 &
                my.response %in% c('PD')) {
        data.plot$Response.organ.mixed[i] = 'uniformProgression'
      } else if(my.stats$notconsistent == 0 & my.stats$consistent > 0 &
                my.response %in% c('SD','PR','CR')) {
        data.plot$Response.organ.mixed[i] = 'uniformDiseaseControl'
      } else {
        data.plot$Response.organ.mixed[i] = NA
      }
      
    }
    
    table(data.plot[,c('Cohort','Response.organ.mixed')])
    
  }
  
  ## --------------------------------------------
  
  ## draw figures
  if(TRUE) {
    
    ## mixed response 
    if(TRUE) {
      
      ## heatmap - response
      if(TRUE) {
        
        ## prep input; add organ response heterogeneity score 
        if(TRUE) {
          
          data.plot.2 = NULL 
          data.plot.2 = data.plot[,c('Cohort','Deidentified.name',
                                     names(tissues),'Best.Response',
                                     'Response.organ.mixed')]
          
          dim(data.plot.2)
          head(data.plot.2)
          
          ## --------------------------------------------
          
          data.plot.2$Response = paste0('zz',data.plot.2$Best.Response)
          
          ## --------------------------------------------
          
          data.plot.2[data.plot.2=='PD'] = -1
          data.plot.2[data.plot.2=='SD'] = 0
          data.plot.2[data.plot.2=='PR'] = 1
          data.plot.2[data.plot.2=='CR'] = 2
          data.plot.2[data.plot.2=='NE'] = NA
          head(data.plot.2)
          
          data.plot.2[,names(tissues)] = apply(
            data.plot.2[,names(tissues)] , 2, function(x) as.numeric(as.character(x))
          )
          head(data.plot.2)
          
          ## --------------------------------------------
          
          data.plot.2$Response = gsub('^zz','',data.plot.2$Response)
          
          ## --------------------------------------------
          
          ## sort by cohort first and then by overall response 
          data.plot.2 = data.plot.2[order(data.plot.2$Cohort,
                                          data.plot.2$Best.Response),]
          data.plot.2 = rbind(
            data.plot.2[data.plot.2$Cohort=='pd1',],
            data.plot.2[data.plot.2$Cohort=='ipinivo',],
            data.plot.2[data.plot.2$Cohort=='braf',]
            
          )
          
          head(data.plot.2)
          
          write.csv(data.plot.2,
                    file = paste0(clinical.file,'.ptByResponseOrgan.heatmap.csv'))
          
        }
        
        ## --------------------------------------------
        
        ## plotting
        if(TRUE) {
          
          centered.data = NULL 
          centered.data = data.plot.2[,c('Deidentified.name',names(tissues))]
          row.names(centered.data) = centered.data$Deidentified.name
          centered.data = data.frame(t(centered.data[,-1]))
          centered.data[,1:10]
          
          write.csv(centered.data,
                    file = paste0(clinical.file,'.ptByResponseOrgan.heatmap.data.csv'))
          
          if(2==2) {
            
            ## add annotation 
            sample.anno = NULL 
            sample.anno = data.plot.2[,c('Cohort','Deidentified.name','Response',
                                         'Response.organ.mixed')]
            row.names(sample.anno) = gsub('-','.',sample.anno$Deidentified.name)
            
            sample.anno[is.na(sample.anno)] = NA
            
            ## sort sample anno same as expression matrix
            sample.anno = sample.anno[order(match(row.names(sample.anno),
                                                  colnames(centered.data))),]
            print(all.equal(row.names(sample.anno), colnames(centered.data)))
            # TRUE 
            
            sample.anno.colors = list(
              Cohort = plot.colors,
              Response.organ.mixed = plot.colors,
              Response =  plot.colors
            )
            
            plot.anno = HeatmapAnnotation(df = sample.anno[,c('Cohort',
                                                              'Response.organ.mixed',
                                                              'Response'),drop=F]
                                          ,col = sample.anno.colors
            )
            
            write.csv(sample.anno,
                      file = paste0(clinical.file,'.ptByResponseOrgan.heatmap.anno.csv'))
            
          }
          
          if(3==3) {
            
            my.heatmap.colors = colorRamp2(c(-1, 0, 1, 2), 
                                           c(plot.colors['PD'], plot.colors['SD'],
                                             plot.colors['PR'], plot.colors['CR']))
            
            col.title = paste0(ncol(centered.data),' samples')
            row.title = paste0(nrow(centered.data),' organ sites')
            
            p2 = NULL 
            p2 =Heatmap(centered.data,
                        col = my.heatmap.colors,
                        name = 'OrganResponse',
                        na_col = "#FFFFFF",
                        column_title = col.title, 
                        row_title = row.title,
                        # column_title_side = 'bottom',
                        column_dend_height = unit(2, "cm"), 
                        row_dend_width = unit(4, "cm"), 
                        # km = 2,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        clustering_distance_rows = "euclidean",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "ward.D2",
                        show_row_names = T,
                        show_column_names = F,
                        column_names_gp =  gpar(fontsize = 22),
                        top_annotation = plot.anno)
            
            
            pdf(file = paste0(clinical.file,'.ptByResponseOrgan.heatmap.pdf'), 
                width = 20, height = 4)
            print(p2)
            dev.off()
            
          }
          
        }
      }
      
      ## --------------------------------------------
      
      ## heatmap - progression/disease control
      if(TRUE) {
        
        ## prep input; add organ response heterogeneity score 
        if(TRUE) {
          
          data.plot.2 = NULL 
          data.plot.2 = data.plot[,c('Cohort','Deidentified.name',
                                     names(tissues),'Best.Response',
                                     'Response.organ.mixed')]
          
          dim(data.plot.2)
          head(data.plot.2)
          
          ## --------------------------------------------
          
          data.plot.2$Response = paste0('zz',data.plot.2$Best.Response)
          
          ## --------------------------------------------
          
          data.plot.2[data.plot.2=='PD'] = -1
          data.plot.2[data.plot.2=='SD'] = 1
          data.plot.2[data.plot.2=='PR'] = 1
          data.plot.2[data.plot.2=='CR'] = 1
          data.plot.2[data.plot.2=='NE'] = NA
          head(data.plot.2)
          
          data.plot.2[,names(tissues)] = apply(
            data.plot.2[,names(tissues)] , 2, function(x) as.numeric(as.character(x))
          )
          head(data.plot.2)
          
          
          ## --------------------------------------------
          
          ## organ response heterogeneity 
          data.plot.2$Response.organ.het = apply(data.plot.2[,names(tissues)],1,sum,na.rm=T)
          table(data.plot.2$Response.organ)
          
          data.plot.2$Response.organ.notNA = apply(data.plot.2[,names(tissues)],1,
                                                   function(x) sum(!is.na(x)))
          
          data.plot.2$Response.organ.het.norm = 
            data.plot.2$Response.organ.het / data.plot.2$Response.organ.notNA
          table(data.plot.2$Response.organ.het.norm)
          data.plot.2$Response.organ.het.norm[
            is.infinite(data.plot.2$Response.organ.het.norm)] = 
            NA
          data.plot.2$Response.organ.het.norm[
            is.nan(data.plot.2$Response.organ.het.norm)] = 
            NA
          
          data.plot.2[is.na(data.plot.2$Response.organ.het.norm),]
          
          ## --------------------------------------------
          
          data.plot.2$Response = gsub('^zz','',data.plot.2$Response)
          data.plot.2$Response[which(data.plot.2$Response %in% c('PD'))] = 'PD'
          data.plot.2$Response[which(data.plot.2$Response %in% c('SD','PR','CR'))] = 'DC'
          
          ## --------------------------------------------
          
          ## sort by cohort first and then by overall response 
          data.plot.2 = data.plot.2[order(data.plot.2$Cohort,
                                          data.plot.2$Best.Response),]
          data.plot.2 = rbind(
            data.plot.2[data.plot.2$Cohort=='pd1',],
            data.plot.2[data.plot.2$Cohort=='ipinivo',],
            data.plot.2[data.plot.2$Cohort=='braf',]
            
          )
          
          head(data.plot.2)
          
          write.csv(data.plot.2,
                    file = paste0(clinical.file,'.ptByResponseOrganDCR.heatmap.csv'))
          
        }
        
        ## --------------------------------------------
        
        ## plotting
        if(TRUE) {
          
          centered.data = NULL 
          centered.data = data.plot.2[,c('Deidentified.name',names(tissues))]
          row.names(centered.data) = centered.data$Deidentified.name
          centered.data = data.frame(t(centered.data[,-1]))
          centered.data[,1:10]
          
          write.csv(centered.data,
                    file = paste0(clinical.file,'.ptByResponseOrganDCR.heatmap.data.csv'))
          
          if(2==2) {
            
            ## add annotation 
            sample.anno = NULL 
            sample.anno = data.plot.2[,c('Cohort','Deidentified.name','Response',
                                         'Response.organ.het.norm',
                                         'Response.organ.mixed')]
            row.names(sample.anno) = gsub('-','.',sample.anno$Deidentified.name)
            
            sample.anno[is.na(sample.anno)] = NA
            
            ## sort sample anno same as expression matrix
            sample.anno = sample.anno[order(match(row.names(sample.anno),
                                                  colnames(centered.data))),]
            print(all.equal(row.names(sample.anno), colnames(centered.data)))
            # TRUE 
            
            min(sample.anno$Response.organ.het.norm, na.rm=T) 
            max(sample.anno$Response.organ.het.norm, na.rm=T) 
            
            
            sample.anno.colors = list(
              Cohort = plot.colors,
              Response.organ.het.norm = Response.organ.het.norm.colors,
              Response.organ.mixed = plot.colors,
              Response =  plot.colors
            )
            
            plot.anno = HeatmapAnnotation(df = sample.anno[,c('Cohort',
                                                              'Response.organ.het.norm',
                                                              'Response.organ.mixed',
                                                              'Response'),drop=F]
                                          ,col = sample.anno.colors
            )
            
            write.csv(sample.anno,
                      file = paste0(clinical.file,'.ptByResponseOrganDCR.heatmap.anno.csv'))
            
          }
          
          if(3==3) {
            
            my.heatmap.colors = colorRamp2(c(-1, 1), 
                                           c(plot.colors['PD'], plot.colors['DC']))
            
            col.title = paste0(ncol(centered.data),' samples')
            row.title = paste0(nrow(centered.data),' organ sites')
            
            p2 = NULL 
            p2 =Heatmap(centered.data,
                        col = my.heatmap.colors,
                        name = 'OrganResponse',
                        na_col = "#FFFFFF",
                        column_title = col.title, 
                        row_title = row.title,
                        # column_title_side = 'bottom',
                        column_dend_height = unit(2, "cm"), 
                        row_dend_width = unit(4, "cm"), 
                        # km = 2,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        clustering_distance_rows = "euclidean",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "ward.D2",
                        show_row_names = T,
                        show_column_names = F,
                        column_names_gp =  gpar(fontsize = 22),
                        top_annotation = plot.anno)
            
            
            pdf(file = paste0(clinical.file,'.ptByResponseOrganDCR.heatmap.pdf'), 
                width = 20, height = 4)
            print(p2)
            dev.off()
            
          }
          
        }
      }
      
      ## --------------------------------------------
      
      ## bargraph - response
      if(TRUE) {
        
        data.plot.2 = NULL 
        
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          ## --------------------------------------------
          
          my.df = NULL
          
          if(my.cohort == 'ici') {
            my.df = data.plot[
              data.plot$Cohort2 ==my.cohort, c('Cohort2','Best.Response.by.RECIST',
                                               'Response.organ.mixed')]
            colnames(my.df)[1] = 'Cohort'
          } else {
            my.df =data.plot[
              data.plot$Cohort==my.cohort,c('Cohort','Best.Response.by.RECIST',
                                            'Response.organ.mixed')]
          }
          
          print(table(my.df$Response.organ.mixed))
          
          print(sum(is.na(my.df$Response.organ.mixed)))
          
          ## convert to %
          my.df = rbind(
            data.frame(prop.table(table(
              my.df[my.df$Best.Response.by.RECIST == 'PD',
                    c('Cohort','Best.Response.by.RECIST','Response.organ.mixed')]))),
            data.frame(prop.table(table(
              my.df[my.df$Best.Response.by.RECIST == 'SD',
                    c('Cohort','Best.Response.by.RECIST','Response.organ.mixed')]))),
            data.frame(prop.table(table(
              my.df[my.df$Best.Response.by.RECIST == 'PR',
                    c('Cohort','Best.Response.by.RECIST','Response.organ.mixed')]))),
            data.frame(prop.table(table(
              my.df[my.df$Best.Response.by.RECIST == 'CR',
                    c('Cohort','Best.Response.by.RECIST','Response.organ.mixed')])))
            
          )
          
          ## --------------------------------------------
          
          data.plot.2 = rbind(data.plot.2,
                              my.df)
        }
        
        data.plot.2
        
        ## --------------------------------------------
        
        data.plot.2$Pct = round(data.plot.2$Freq * 100, 
                                digits = 0)
        
        data.plot.2$Best.Response.by.RECIST = factor(data.plot.2$Best.Response.by.RECIST,
                                                     levels = c('PD','SD','PR','CR',
                                                                'NA'))
        
        write.csv(data.plot.2,
                  file = paste0(clinical.file,'.ptByResponseOrganMixed.bar.csv'))
        
        
        p1 = NULL 
        p1 = ggplot(data.plot.2, aes(Best.Response.by.RECIST,
                                     Pct)) +
          geom_bar(aes(fill = Response.organ.mixed), stat = 'identity',
                   position = 'stack',
                   color = '#000000', width = 0.6) +
          scale_fill_manual(values = plot.colors) + 
          facet_grid(~ Cohort, scales = 'free_x') +
          theme_pubr(x.text.angle = 90) 
        
        p1
        
        pdf(file = paste0(clinical.file,'.ptByResponseOrganMixed.bar.pdf'),
            height = 4, width = 10)
        print(p1)
        dev.off()
        
        
        
      }
      
      ## --------------------------------------------
      
      ## bargraph - progression/disease control
      if(TRUE) {
        
        data.plot.2 = NULL 
        
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          ## --------------------------------------------
          
          my.df = NULL
          
          if(my.cohort == 'ici') {
            my.df = data.plot[
              data.plot$Cohort2 ==my.cohort, c('Cohort2','Best.Response.by.RECIST',
                                               'Response.organ.mixed')]
            colnames(my.df)[1] = 'Cohort'
          } else {
            my.df =data.plot[
              data.plot$Cohort==my.cohort,c('Cohort','Best.Response.by.RECIST',
                                            'Response.organ.mixed')]
          }
          
          my.df$Best.Response.by.RECIST.clp = NA
          my.df$Best.Response.by.RECIST.clp[which(
            my.df$Best.Response.by.RECIST %in% c('PD')
          )] = "PD"
          my.df$Best.Response.by.RECIST.clp[which(
            my.df$Best.Response.by.RECIST %in% c('SD','PR','CD')
          )] = "DC"
          
          print(table(my.df$Best.Response.by.RECIST.clp))
          
          print(table(my.df$Response.organ.mixed))
          
          print(sum(is.na(my.df$Response.organ.mixed)))
          
          ## convert to %
          my.df = rbind(
            data.frame(prop.table(table(
              my.df[my.df$Best.Response.by.RECIST.clp == 'PD',
                    c('Cohort','Best.Response.by.RECIST.clp','Response.organ.mixed')]))),
            data.frame(prop.table(table(
              my.df[my.df$Best.Response.by.RECIST.clp == 'DC',
                    c('Cohort','Best.Response.by.RECIST.clp','Response.organ.mixed')])))
            
          )
          
          ## --------------------------------------------
          
          data.plot.2 = rbind(data.plot.2,
                              my.df)
        }
        
        data.plot.2
        
        ## --------------------------------------------
        
        data.plot.2$Pct = round(data.plot.2$Freq * 100, 
                                digits = 0)
        
        data.plot.2$Best.Response.by.RECIST.clp = factor(
          data.plot.2$Best.Response.by.RECIST.clp,
          levels = c('PD','DC',
                     'NA'))
        
        write.csv(data.plot.2,
                  file = paste0(clinical.file,'.ptByResponseOrganMixedDCR.bar.csv'))
        
        
        p1 = NULL 
        p1 = ggplot(data.plot.2, aes(Best.Response.by.RECIST.clp,
                                     Pct)) +
          geom_bar(aes(fill = Response.organ.mixed), stat = 'identity',
                   position = 'stack',
                   color = '#000000', width = 0.6) +
          scale_fill_manual(values = plot.colors) + 
          facet_grid(~ Cohort, scales = 'free_x') +
          theme_pubr(x.text.angle = 90) 
        
        p1
        
        pdf(file = paste0(clinical.file,'.ptByResponseOrganMixedDCR.bar.pdf'),
            height = 4, width = 10)
        print(p1)
        dev.off()
        
        
        
      }
    }
    
  }
  
  ## --------------------------------------------
  
  ## compute stats
  if(TRUE) {
    
    ## mixed response 
    if(TRUE) {
      
      data.plot.2 = NULL 
      data.plot.2 = data.plot
      
      data.plot.2$Best.Response.by.RECIST.clp = NA
      data.plot.2$Best.Response.by.RECIST.clp[which(
        data.plot.2$Best.Response.by.RECIST %in% c('PD')
      )] = 'PD'
      data.plot.2$Best.Response.by.RECIST.clp[which(
        data.plot.2$Best.Response.by.RECIST %in% c('SD','PR','CR')
      )] = 'DC'
      
      data.plot.2$Response.organ.mixed[ data.plot.2$Response.organ.mixed %in% 
                                          c('uniformProgression','uniformDiseaseControl')] = 'uniformResponse'
      
      dim(data.plot.2) 
      
      data.stats = NULL 
      for(my.cohort in c('ici','pd1','ipinivo','braf')) {
        print(my.cohort)
        
        ## --------------------------------------------
        
        my.df = NULL 
        
        if(my.cohort == 'ici') {
          my.df = data.plot.2[data.plot.2$Cohort2==my.cohort,]
          colnames(my.df)[1] = 'Cohort'
        } else {
          my.df = data.plot.2[data.plot.2$Cohort==my.cohort,]
        }
        
        my.df = data.frame(table(my.df[,c('Best.Response.by.RECIST.clp',
                                          'Response.organ.mixed')]))
        my.df = reshape2::dcast(Best.Response.by.RECIST.clp ~ Response.organ.mixed,
                                value.var = 'Freq',
                                data = my.df,
                                fun.aggregate = sum
        )
        row.names(my.df) = my.df$Best.Response.by.RECIST.clp
        my.df = my.df[,-1]
        my.df = as.matrix(my.df)
        my.df = t(my.df)
        my.df = my.df[,c('PD','DC')]
        
        my.stats = NULL 
        my.stats = tidy(fisher.test(my.df))
        
        data.stats = rbind(data.stats,
                           data.frame(
                             Cohort = my.cohort,
                             PD.mixedResponse = my.df[1,1],
                             PD.uniformResponse = my.df[2,1],
                             DC.mixedResponse = my.df[1,2],
                             DC.uniformResponse = my.df[2,2],
                             my.stats,
                             stringsAsFactors = F
                           ))
        
      }
      
      write.csv(data.stats,
                file = paste0(clinical.file,'.ptByResponseOrganMixedDCR.fisher.csv'))
    }
    
  }
  
}