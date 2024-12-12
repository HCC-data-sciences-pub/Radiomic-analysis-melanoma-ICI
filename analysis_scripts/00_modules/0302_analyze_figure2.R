
## draw figures 
if(TRUE) {
  
  ## between cohorts: ipinivo vs pd1 vs braf 
  if(TRUE) {
    
    ## prep input
    if(TRUE) {
      data.plot = NULL 
      i = 0
      
      for(my.cohort in cohorts) {
        
        print(my.cohort)
        my.df = NULL 
        my.df = deg[[my.cohort]][['all']]
        
        if(is.null(my.df)) { next }
        if(length(nrow(my.df)) == 0) { next }
        
        row.names(my.df) = my.df$feat
        colnames(my.df) = paste0(my.cohort,'.all.',colnames(my.df))
        
        i = i + 1
        
        if(i == 1) {
          data.plot = my.df
        } else {
          data.plot = merge(data.plot, my.df,
                            by = 'row.names', all = T)
          row.names(data.plot) = data.plot$Row.names
          data.plot = data.plot[,-1]
        }
        
      }
      
      dim(data.plot) # 343  36
      
      data.plot = data.frame(feat = row.names(data.plot),
                             data.plot,
                             stringsAsFactors = F)
      
      write.csv(data.plot,
                file = paste0(output, '.combined.diffcohorts_all.csv'))
      
      ## --------------------------------------------
      
      data.plot.flt = NULL 
      deg.filename = NULL 
      
      if(use.fdr == 1) {
        data.plot.flt = data.plot[(!is.na(data.plot$ipinivo.all.p.adj) & 
                                     data.plot$ipinivo.all.p.adj < fdr) |
                                    (!is.na(data.plot$pd1.all.p.adj) & 
                                       data.plot$pd1.all.p.adj < fdr) | 
                                    (!is.na(data.plot$braf.all.p.adj) & 
                                       data.plot$braf.all.p.adj < fdr),,drop=F]
        
        dim(data.plot.flt) # 39  36
        
        deg.filename = paste0(output, '.combined.diffcohorts_all.fdr',fdr)
        
      }
      if(use.rawp == 1) {
        data.plot.flt = data.plot[(!is.na(data.plot$ipinivo.all.p.value) & 
                                     data.plot$ipinivo.all.p.value < rawp) |
                                    (!is.na(data.plot$pd1.all.p.value) & 
                                       data.plot$pd1.all.p.value < rawp) | 
                                    (!is.na(data.plot$braf.all.p.value) & 
                                       data.plot$braf.all.p.value < rawp),,drop=F]
        
        dim(data.plot.flt) # 29 36
        
        deg.filename = paste0(output, '.combined.diffcohorts_all.rawp',rawp)
        
      }
      
      dim(data.plot.flt) ## 39 36
      print(deg.filename)
      
      write.csv(data.plot.flt,
                file = paste0(deg.filename,'.csv'))
    }
    
    ## --------------------------------------------
    
    if(nrow(data.plot.flt) > 0) { 
      
      ## heatmaps 
      if(TRUE) {
        
        ## ipinivo 
        if(TRUE) {
          
          my.cohort = 'ipinivo'
          
          ## ------------------------------------------------------
          
          my.df = NULL 
          centered_data = NULL 
          p1 = NULL
          my.df = data[[my.cohort]][['all']]
          dim(my.df)
          
          centered_data = t(my.df[,data.plot.flt$feat,drop=F])
          centered_data = t(scale(t(centered_data)))
          dim(centered_data)
          
          ## ------------------------------------------------------
          
          if(TRUE) {
            ## feat order (row order) was derived from ipinivo heatmap clusters (see below)
            heatmap.order = NULL 
            heatmap.order = read.csv(paste0(deg.filename,'.heatmap.row_order.csv'))
            centered_data.score = NULL 
            centered_data.score = centered_data
            ## flip the value ...
            x = NULL 
            x = heatmap.order$feat[heatmap.order$direction %in% c('up in DC')]
            x = which(row.names(centered_data.score) %in% x)
            centered_data.score[x,] = centered_data.score[x,] * -1
            centered_data.score = data.frame(feat.score = 
                                               apply(centered_data.score,
                                                     2, mean, na.rm=T))
            centered_data.score = data.frame(
              Sample = row.names(centered_data.score),
              centered_data.score,
              stringsAsFactors = F
            )
            
            centered_data.score = 
              centered_data.score[rev(order(centered_data.score$feat.score)),,drop=F]
            x = NULL
            y = NULL 
            x = centered_data.score[centered_data.score$Sample %in%
                                      row.names(my.df)[my.df$Response=='PD'],,drop=F]
            dim(x)
            x = x[rev(order(x$feat.score)),,drop=F]
            y = centered_data.score[centered_data.score$Sample %in%
                                      row.names(my.df)[my.df$Response=='DC'],,drop=F]
            dim(y)
            y = y[rev(order(y$feat.score)),,drop=F]
            centered_data.score = rbind(x,y)
            centered_data.score$heatmap.col_order = 1:nrow(centered_data.score)
            
            write.csv(centered_data.score,
                      file = paste0(deg.filename,'.heatmap.col_order.csv'))
            
            
          }
          
          ## ------------------------------------------------------
          
          centered_data = centered_data[,centered_data.score$Sample]
          print(all.equal(colnames(centered_data),
                          centered_data.score$Sample))
          ## TRUE
          
          ## ------------------------------------------------------
          
          if(TRUE) {
            
            ## add annotation 
            sample.anno = NULL 
            sample.anno = my.df[,1,drop=F]
            
            ## sort sample anno same as expression matrix
            sample.anno = sample.anno[row.names(sample.anno) %in% colnames(centered_data),,
                                      drop=F]
            ## sot samples 
            sample.anno = sample.anno[order(match(row.names(sample.anno),
                                                  colnames(centered_data))),,drop=F]
            print(all.equal(row.names(sample.anno), colnames(centered_data)))
            #3 TRUE 
            
            sample.anno.colors = list(
              Response = plot.colors
            )
            
            plot.anno = HeatmapAnnotation(df = sample.anno[,!colnames(sample.anno)%in%
                                                             c('Sample','Subject.ID'),
                                                           drop=F]
                                          ,col = sample.anno.colors
            )
            
          }
          
          write.csv(sample.anno,
                    file = paste0(deg.filename, '.heatmap.anno.csv'))
          
          write.csv(centered_data,
                    file = paste0(deg.filename, '.heatmap.data.csv'))
          
          ## ------------------------------------------------------
          
          myheatcol = colorRamp2(c(-2, 0, 2), c("blue","white", "red"))
          col.title = paste0(my.cohort,' ', my.tissue, " ",ncol(centered_data), ' samples')
          row.title = paste0(nrow(centered_data), ' features')
          
          p1 = NULL 
          p1 = Heatmap(centered_data,
                       na_col = "#000000",
                       col = myheatcol,
                       rect_gp = gpar(col = NA),
                       show_heatmap_legend = T,
                       column_title = col.title,
                       row_title = row.title,
                       # column_title_side = 'bottom',
                       column_names_side = 'bottom',
                       row_dend_width = unit(5, "cm"),
                       column_dend_height = unit(5, "cm"),
                       # km = 2,
                       cluster_rows = T,
                       cluster_columns = T,
                       clustering_distance_rows = "euclidean",
                       clustering_method_rows = "ward.D2",
                       clustering_distance_columns = "euclidean",
                       clustering_method_columns = "ward.D2",
                       show_row_names = T,
                       show_column_names = T,
                       top_annotation = plot.anno,
                       heatmap_legend_param = list(title = 'Cell density', 
                                                   color_bar = "continuous")
          )
          
          p1
          
          p2 = NULL 
          p2 = Heatmap(centered_data,
                       na_col = "#000000",
                       col = myheatcol,
                       rect_gp = gpar(col = NA),
                       show_heatmap_legend = T,
                       column_title = col.title,
                       row_title = row.title,
                       # column_title_side = 'bottom',
                       column_names_side = 'bottom',
                       row_dend_width = unit(5, "cm"),
                       column_dend_height = unit(5, "cm"),
                       # km = 2,
                       cluster_rows = T,
                       cluster_columns = F,
                       clustering_distance_rows = "euclidean",
                       clustering_method_rows = "ward.D2",
                       clustering_distance_columns = "euclidean",
                       clustering_method_columns = "ward.D2",
                       show_row_names = T,
                       show_column_names = T,
                       top_annotation = plot.anno,
                       heatmap_legend_param = list(title = 'Cell density', 
                                                   color_bar = "continuous")
          )
          
          p2
          
          pdf(file = paste0(deg.filename, '.heatmap.pdf'), 
              width = 10, height = 10)
          print(p1)
          dev.off()
          
          pdf(file = paste0(deg.filename, '.heatmap.srtBySM.pdf'), 
              width = 10, height = 10)
          print(p2)
          dev.off()
          
        }
        
      }
      
      ## --------------------------------------------
      
      ## bubble plot...
      if(TRUE) {
        
        heatmap.order = read.csv(paste0(deg.filename,'.heatmap.row_order.csv'))
        heatmap.order = heatmap.order[order(heatmap.order$heatmap.row_order),,drop=F]
        
        dim(heatmap.order) ## 39  4
        dim(data.plot.flt) ## 39 37
        
        
        data.plot.2 = NULL 
        data.plot.2 = data.plot.flt
        
        x = NULL 
        y = NULL 
        x = data.plot.2[,c('feat', paste0(cohorts,'.all.p.value'))]
        y = data.plot.2[,c('feat',paste0(cohorts,'.all.mean.diff'))]
        x = reshape2::melt(x)
        y = reshape2::melt(y)
        x$variable = as.character(x$variable)
        y$variable = as.character(y$variable)
        
        data.plot.2 = data.frame(
          feat = x$feat,
          Cohort = gsub('.all.p.value$','',x$variable),
          Tissue = 'all',
          p.value = x$value,
          mean.diff = y$value,
          stringsAsFactors = F
          
        )
        
        data.plot.2$direction = NA
        data.plot.2$direction[which(data.plot.2$mean.diff>0)] = 'DC'
        data.plot.2$direction[which(data.plot.2$mean.diff<0)] = 'PD'
        
        data.plot.2$feat = factor(data.plot.2$feat,
                                  levels = as.character(rev(heatmap.order$feat)))
        
        data.plot.2$Cohort = factor(data.plot.2$Cohort,
                                    levels = c('ipinivo','pd1','braf'))
        
        p1 = NULL 
        p1 = ggplot(data.plot.2, aes(Cohort, feat)) +
          geom_point(aes(size = (-log10(p.value))^2,
                         fill = direction),
                     shape = 21) +
          scale_fill_manual(values = plot.colors) +
          theme_minimal()
        
        p1
        
        pdf(file = paste0(deg.filename,'.bubble.pdf'),
            width = 5, height = 6)
        print(p1)
        dev.off()
        
        data.plot.2$Key = paste0(data.plot.2$Tissue,'!',data.plot.2$Cohort)
        data.plot.2$Key = factor(data.plot.2$Key,
                                 levels = c('all!ipinivo','all!pd1','all!braf'))
        
        p2 = NULL 
        p2 = ggplot(data.plot.2, aes(-log10(p.value), feat)) +
          geom_segment( aes(y=feat, yend=feat, x=0, xend=-log10(p.value))) +
          geom_point(aes(fill = direction), shape=21,
                     size = 3) +
          scale_fill_manual(values = plot.colors) +
          theme_minimal() +
          facet_grid(~ Key) +
          geom_vline(xintercept = -log10(0.05), linetype = 'dashed',color = '#CC0000') +
          geom_vline(xintercept = -log10(1), linetype = 'solid',color = '#000000')
        
        
        p2
        
        pdf(file = paste0(deg.filename,'.lollipop.pdf'),
            width = 7, height = 5)
        print(p2)
        dev.off()
        
        
      }
      
    }
    
  }
  
  ## --------------------------------------------
  
  ## same cohort, all versus different organs 
  if(TRUE) {
    
    ## ipinivo
    if(TRUE) {
      
      my.cohort = 'ipinivo' ## ipinivo pd1
      print(my.cohort)
      
      ## --------------------------------------------
      
      ## prep input
      if(TRUE) {
        
        
        data.plot = NULL 
        i = 0
        
        for(my.tissue in tissues) {
          
          print(my.tissue)
          my.df = NULL 
          my.df = deg[[my.cohort]][[my.tissue]]
          
          if(is.null(my.df)) { next }
          if(length(nrow(my.df)) == 0) { next }
          
          row.names(my.df) = my.df$feat
          colnames(my.df) = paste0(my.cohort,'.',my.tissue,'.',colnames(my.df))
          
          i = i + 1
          
          if(i == 1) {
            data.plot = my.df
          } else {
            data.plot = merge(data.plot, my.df,
                              by = 'row.names', all = T)
            row.names(data.plot) = data.plot$Row.names
            data.plot = data.plot[,-1]
          }
          
        }
        
        dim(data.plot) # 329  72
        
        data.plot = data.frame(feat = row.names(data.plot),
                               data.plot,
                               stringsAsFactors = F)
        
        write.csv(data.plot,
                  file = paste0(output, '.combined.difforgans_',my.cohort,'.csv'))
        
        ## --------------------------------------------
        
        data.plot.flt = NULL 
        deg.filename = NULL 
        
        if(use.fdr == 1) {
          data.plot.flt = data.plot[(!is.na(data.plot$ipinivo.all.p.adj) & 
                                       data.plot$ipinivo.all.p.adj < fdr) 
                                    ,,drop=F]
          
          dim(data.plot.flt) # 39  72
          
          deg.filename = paste0(output, '.combined.difforgans_',
                                my.cohort,'.fdr',fdr)
          
        }
        if(use.rawp == 1) {
          data.plot.flt = data.plot[(!is.na(data.plot$ipinivo.all.p.value) & 
                                       data.plot$ipinivo.all.p.value < rawp) 
                                    ,,drop=F]
          
          dim(data.plot.flt) # 27 73
          
          deg.filename = paste0(output, '.combined.difforgans_',
                                my.cohort,'.rawp',rawp)
          
        }
        
        dim(data.plot.flt) ## 27 73
        print(deg.filename)
        
        write.csv(data.plot.flt,
                  file = paste0(deg.filename,'.csv'))
      }
      
      ## --------------------------------------------
      
      if(nrow(data.plot.flt) > 0) { 
        
        ## heatmaps 
        if(TRUE) {
          
          ## ipinivo 
          if(TRUE) {
            
            my.cohort = 'ipinivo'
            
            ## ------------------------------------------------------
            
            my.df = NULL 
            centered_data = NULL 
            p1 = NULL
            my.df = data[[my.cohort]][['all']]
            dim(my.df)
            
            centered_data = t(my.df[,data.plot.flt$feat,drop=F])
            centered_data = t(scale(t(centered_data)))
            dim(centered_data)
            
            ## ------------------------------------------------------
            
            if(TRUE) {
              ## feat order (row order) was derived from ipinivo heatmap clusters (see below)
              heatmap.order = NULL 
              heatmap.order = read.csv(paste0(deg.filename,'.heatmap.row_order.csv'))
              centered_data.score = NULL 
              centered_data.score = centered_data
              ## flip the value ...
              x = NULL 
              x = heatmap.order$feat[heatmap.order$direction %in% c('up in DC')]
              x = which(row.names(centered_data.score) %in% x)
              centered_data.score[x,] = centered_data.score[x,] * -1
              centered_data.score = data.frame(feat.score = 
                                                 apply(centered_data.score,
                                                       2, mean, na.rm=T))
              centered_data.score = data.frame(
                Sample = row.names(centered_data.score),
                centered_data.score,
                stringsAsFactors = F
              )
              
              centered_data.score = 
                centered_data.score[rev(order(centered_data.score$feat.score)),,drop=F]
              x = NULL
              y = NULL 
              x = centered_data.score[centered_data.score$Sample %in%
                                        row.names(my.df)[my.df$Response=='PD'],,drop=F]
              dim(x)
              x = x[rev(order(x$feat.score)),,drop=F]
              y = centered_data.score[centered_data.score$Sample %in%
                                        row.names(my.df)[my.df$Response=='DC'],,drop=F]
              dim(y)
              y = y[rev(order(y$feat.score)),,drop=F]
              centered_data.score = rbind(x,y)
              centered_data.score$heatmap.col_order = 1:nrow(centered_data.score)
              
              write.csv(centered_data.score,
                        file = paste0(deg.filename,'.heatmap.col_order.csv'))
              
              
            }
            
            ## ------------------------------------------------------
            
            centered_data = centered_data[,centered_data.score$Sample]
            print(all.equal(colnames(centered_data),
                            centered_data.score$Sample))
            ## TRUE
            
            ## ------------------------------------------------------
            
            if(TRUE) {
              
              ## add annotation 
              sample.anno = NULL 
              sample.anno = my.df[,1,drop=F]
              
              ## sort sample anno same as expression matrix
              sample.anno = sample.anno[row.names(sample.anno) %in% colnames(centered_data),,
                                        drop=F]
              ## sot samples 
              sample.anno = sample.anno[order(match(row.names(sample.anno),
                                                    colnames(centered_data))),,drop=F]
              print(all.equal(row.names(sample.anno), colnames(centered_data)))
              #3 TRUE 
              
              sample.anno.colors = list(
                Response = plot.colors
              )
              
              plot.anno = HeatmapAnnotation(df = sample.anno[,!colnames(sample.anno)%in%
                                                               c('Sample','Subject.ID'),
                                                             drop=F]
                                            ,col = sample.anno.colors
              )
              
            }
            
            write.csv(sample.anno,
                      file = paste0(deg.filename, '.heatmap.anno.csv'))
            
            write.csv(centered_data,
                      file = paste0(deg.filename, '.heatmap.data.csv'))
            
            ## ------------------------------------------------------
            
            myheatcol = colorRamp2(c(-2, 0, 2), c("blue","white", "red"))
            col.title = paste0(my.cohort,' ', my.tissue, " ",ncol(centered_data), ' samples')
            row.title = paste0(nrow(centered_data), ' features')
            
            p1 = NULL 
            p1 = Heatmap(centered_data,
                         na_col = "#000000",
                         col = myheatcol,
                         rect_gp = gpar(col = NA),
                         show_heatmap_legend = T,
                         column_title = col.title,
                         row_title = row.title,
                         # column_title_side = 'bottom',
                         column_names_side = 'bottom',
                         row_dend_width = unit(5, "cm"),
                         column_dend_height = unit(5, "cm"),
                         # km = 2,
                         cluster_rows = T,
                         cluster_columns = T,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "ward.D2",
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "ward.D2",
                         show_row_names = T,
                         show_column_names = T,
                         top_annotation = plot.anno,
                         heatmap_legend_param = list(title = 'Cell density', 
                                                     color_bar = "continuous")
            )
            
            p1
            
            p2 = NULL 
            p2 = Heatmap(centered_data,
                         na_col = "#000000",
                         col = myheatcol,
                         rect_gp = gpar(col = NA),
                         show_heatmap_legend = T,
                         column_title = col.title,
                         row_title = row.title,
                         # column_title_side = 'bottom',
                         column_names_side = 'bottom',
                         row_dend_width = unit(5, "cm"),
                         column_dend_height = unit(5, "cm"),
                         # km = 2,
                         cluster_rows = T,
                         cluster_columns = F,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "ward.D2",
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "ward.D2",
                         show_row_names = T,
                         show_column_names = T,
                         top_annotation = plot.anno,
                         heatmap_legend_param = list(title = 'Cell density', 
                                                     color_bar = "continuous")
            )
            
            p2
            
            pdf(file = paste0(deg.filename, '.heatmap.pdf'), 
                width = 10, height = 10)
            print(p1)
            dev.off()
            
            pdf(file = paste0(deg.filename, '.heatmap.srtBySM.pdf'), 
                width = 10, height = 10)
            print(p2)
            dev.off()
            
          }
          
        }
        
        ## --------------------------------------------
        
        ## bubble plot...
        if(TRUE) {
          
          heatmap.order = read.csv(paste0(deg.filename,'.heatmap.row_order.csv'))
          heatmap.order = heatmap.order[order(heatmap.order$heatmap.row_order),,drop=F]
          
          dim(heatmap.order) ## 27  4
          dim(data.plot.flt) ## 27 73
          
          
          data.plot.2 = NULL 
          data.plot.2 = data.plot.flt
          
          x = NULL 
          y = NULL 
          x = data.plot.2[,colnames(data.plot.2) %in%
                            c('feat',paste0(my.cohort,'.',tissues,'.p.value'))]
          y = data.plot.2[,colnames(data.plot.2) %in%
                            c('feat',paste0(my.cohort,'.',tissues,'.mean.diff'))]
          x = reshape2::melt(x)
          y = reshape2::melt(y)
          x$variable = as.character(x$variable)
          y$variable = as.character(y$variable)
          
          data.plot.2 = data.frame(
            feat = x$feat,
            Cohort = my.cohort,
            Tissue = gsub(paste0(my.cohort,'.'),'',
                          gsub('.p.value$','',x$variable)),
            p.value = x$value,
            mean.diff = y$value,
            stringsAsFactors = F
            
          )
          
          data.plot.2$direction = NA
          data.plot.2$direction[which(data.plot.2$mean.diff>0)] = 'DC'
          data.plot.2$direction[which(data.plot.2$mean.diff<0)] = 'PD'
          
          data.plot.2$feat = factor(data.plot.2$feat,
                                    levels = as.character(rev(heatmap.order$feat)))
          
          data.plot.2$Tissue = factor(data.plot.2$Tissue,
                                      levels = c('all','Lung_L','LN','Liver_L',
                                                 'Brain','SoftTissue'))
          
          p1 = NULL 
          p1 = ggplot(data.plot.2, aes(Tissue, feat)) +
            geom_point(aes(size = (-log10(p.value))^2,
                           fill = direction),
                       shape = 21) +
            scale_fill_manual(values = plot.colors) +
            theme_minimal()
          
          p1
          
          pdf(file = paste0(deg.filename,'.bubble.pdf'),
              width = 6, height = 6)
          print(p1)
          dev.off()
          
          data.plot.2$Key = paste0(data.plot.2$Tissue,'!',data.plot.2$Cohort)
          data.plot.2$Key = factor(data.plot.2$Key,
                                   levels = c('all!ipinivo','Lung_L!ipinivo',
                                              'LN!ipinivo','Liver_L!ipinivo',
                                              'Brain!ipinivo','SoftTissue!ipinivo'))
          p2 = NULL 
          p2 = ggplot(data.plot.2, aes(-log10(p.value), feat)) +
            geom_segment( aes(y=feat, yend=feat, x=0, xend=-log10(p.value))) +
            geom_point(aes(fill = direction), shape=21,
                       size = 3) +
            scale_fill_manual(values = plot.colors) +
            theme_minimal() +
            facet_grid(~ Key) +
            geom_vline(xintercept = -log10(0.05), linetype = 'dashed',color = '#CC0000') +
            geom_vline(xintercept = -log10(1), linetype = 'solid',color = '#000000')
          
          
          p2
          
          pdf(file = paste0(deg.filename,'.lollipop.pdf'),
              width = 10, height = 8)
          print(p2)
          dev.off()
          
          
        }
        
      }
      
    }
    
    ## --------------------------------------------
    
    ## pd1
    if(TRUE) {
      
      my.cohort = 'pd1' ## pd1 
      print(my.cohort)
      
      ## --------------------------------------------
      
      ## prep input
      if(TRUE) {
        
        
        data.plot = NULL 
        i = 0
        
        for(my.tissue in tissues) {
          
          print(my.tissue)
          my.df = NULL 
          my.df = deg[[my.cohort]][[my.tissue]]
          
          if(is.null(my.df)) { next }
          if(length(nrow(my.df)) == 0) { next }
          
          row.names(my.df) = my.df$feat
          colnames(my.df) = paste0(my.cohort,'.',my.tissue,'.',colnames(my.df))
          
          i = i + 1
          
          if(i == 1) {
            data.plot = my.df
          } else {
            data.plot = merge(data.plot, my.df,
                              by = 'row.names', all = T)
            row.names(data.plot) = data.plot$Row.names
            data.plot = data.plot[,-1]
          }
          
        }
        
        dim(data.plot) # 350  60
        
        data.plot = data.frame(feat = row.names(data.plot),
                               data.plot,
                               stringsAsFactors = F)
        
        write.csv(data.plot,
                  file = paste0(output, '.combined.difforgans_',my.cohort,'.csv'))
        
        ## --------------------------------------------
        
        data.plot.flt = NULL 
        deg.filename = NULL 
        
        if(use.fdr == 1) {
          data.plot.flt = data.plot[(!is.na(data.plot$pd1.all.p.adj) & 
                                       data.plot$pd1.all.p.adj < fdr) 
                                    ,,drop=F]
          
          dim(data.plot.flt) #  5  61
          
          deg.filename = paste0(output, '.combined.difforgans_',
                                my.cohort,'.fdr',fdr)
          
        }
        if(use.rawp == 1) {
          data.plot.flt = data.plot[(!is.na(data.plot$pd1.all.p.value) & 
                                       data.plot$pd1.all.p.value < rawp) 
                                    ,,drop=F]
          
          dim(data.plot.flt) # 10 61
          
          deg.filename = paste0(output, '.combined.difforgans_',
                                my.cohort,'.rawp',rawp)
          
        }
        
        dim(data.plot.flt) ## 5  61
        print(deg.filename)
        
        write.csv(data.plot.flt,
                  file = paste0(deg.filename,'.csv'))
      }
      
      ## --------------------------------------------
      
      if(nrow(data.plot.flt) > 0) { 
        
        ## heatmaps 
        if(TRUE) {
          
          ## pd1 
          if(TRUE) {
            
            my.cohort = 'pd1'
            
            ## ------------------------------------------------------
            
            my.df = NULL 
            centered_data = NULL 
            p1 = NULL
            my.df = data[[my.cohort]][['all']]
            dim(my.df)
            
            centered_data = t(my.df[,data.plot.flt$feat,drop=F])
            centered_data = t(scale(t(centered_data)))
            dim(centered_data)
            
            ## ------------------------------------------------------
            
            if(TRUE) {
              ## feat order (row order) was derived from pd1 heatmap clusters (see below)
              heatmap.order = NULL 
              heatmap.order = read.csv(paste0(deg.filename,'.heatmap.row_order.csv'))
              centered_data.score = NULL 
              centered_data.score = centered_data
              ## flip the value ...
              x = NULL 
              x = heatmap.order$feat[heatmap.order$direction %in% c('up in DC')]
              x = which(row.names(centered_data.score) %in% x)
              centered_data.score[x,] = centered_data.score[x,] * -1
              centered_data.score = data.frame(feat.score = 
                                                 apply(centered_data.score,
                                                       2, mean, na.rm=T))
              centered_data.score = data.frame(
                Sample = row.names(centered_data.score),
                centered_data.score,
                stringsAsFactors = F
              )
              
              centered_data.score = 
                centered_data.score[rev(order(centered_data.score$feat.score)),,drop=F]
              x = NULL
              y = NULL 
              x = centered_data.score[centered_data.score$Sample %in%
                                        row.names(my.df)[my.df$Response=='PD'],,drop=F]
              dim(x)
              x = x[rev(order(x$feat.score)),,drop=F]
              y = centered_data.score[centered_data.score$Sample %in%
                                        row.names(my.df)[my.df$Response=='DC'],,drop=F]
              dim(y)
              y = y[rev(order(y$feat.score)),,drop=F]
              centered_data.score = rbind(x,y)
              centered_data.score$heatmap.col_order = 1:nrow(centered_data.score)
              
              write.csv(centered_data.score,
                        file = paste0(deg.filename,'.heatmap.col_order.csv'))
              
              
            }
            
            ## ------------------------------------------------------
            
            centered_data = centered_data[,centered_data.score$Sample]
            print(all.equal(colnames(centered_data),
                            centered_data.score$Sample))
            ## TRUE
            
            ## ------------------------------------------------------
            
            if(TRUE) {
              
              ## add annotation 
              sample.anno = NULL 
              sample.anno = my.df[,1,drop=F]
              
              ## sort sample anno same as expression matrix
              sample.anno = sample.anno[row.names(sample.anno) %in% colnames(centered_data),,
                                        drop=F]
              ## sot samples 
              sample.anno = sample.anno[order(match(row.names(sample.anno),
                                                    colnames(centered_data))),,drop=F]
              print(all.equal(row.names(sample.anno), colnames(centered_data)))
              #3 TRUE 
              
              sample.anno.colors = list(
                Response = plot.colors
              )
              
              plot.anno = HeatmapAnnotation(df = sample.anno[,!colnames(sample.anno)%in%
                                                               c('Sample','Subject.ID'),
                                                             drop=F]
                                            ,col = sample.anno.colors
              )
              
            }
            
            write.csv(sample.anno,
                      file = paste0(deg.filename, '.heatmap.anno.csv'))
            
            write.csv(centered_data,
                      file = paste0(deg.filename, '.heatmap.data.csv'))
            
            ## ------------------------------------------------------
            
            myheatcol = colorRamp2(c(-2, 0, 2), c("blue","white", "red"))
            col.title = paste0(my.cohort,' ', my.tissue, " ",ncol(centered_data), ' samples')
            row.title = paste0(nrow(centered_data), ' features')
            
            p1 = NULL 
            p1 = Heatmap(centered_data,
                         na_col = "#000000",
                         col = myheatcol,
                         rect_gp = gpar(col = NA),
                         show_heatmap_legend = T,
                         column_title = col.title,
                         row_title = row.title,
                         # column_title_side = 'bottom',
                         column_names_side = 'bottom',
                         row_dend_width = unit(5, "cm"),
                         column_dend_height = unit(5, "cm"),
                         # km = 2,
                         cluster_rows = T,
                         cluster_columns = T,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "ward.D2",
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "ward.D2",
                         show_row_names = T,
                         show_column_names = T,
                         top_annotation = plot.anno,
                         heatmap_legend_param = list(title = 'Cell density', 
                                                     color_bar = "continuous")
            )
            
            p1
            
            p2 = NULL
            p2 = Heatmap(centered_data,
                         na_col = "#000000",
                         col = myheatcol,
                         rect_gp = gpar(col = NA),
                         show_heatmap_legend = T,
                         column_title = col.title,
                         row_title = row.title,
                         # column_title_side = 'bottom',
                         column_names_side = 'bottom',
                         row_dend_width = unit(5, "cm"),
                         column_dend_height = unit(5, "cm"),
                         # km = 2,
                         cluster_rows = T,
                         cluster_columns = F,
                         clustering_distance_rows = "euclidean",
                         clustering_method_rows = "ward.D2",
                         clustering_distance_columns = "euclidean",
                         clustering_method_columns = "ward.D2",
                         show_row_names = T,
                         show_column_names = T,
                         top_annotation = plot.anno,
                         heatmap_legend_param = list(title = 'Cell density',
                                                     color_bar = "continuous")
            )
            
            p2
            
            pdf(file = paste0(deg.filename, '.heatmap.pdf'), 
                width = 10, height = 10)
            print(p1)
            dev.off()
            
            pdf(file = paste0(deg.filename, '.heatmap.srtBySM.pdf'),
                width = 10, height = 10)
            print(p2)
            dev.off()
            
          }
          
        }
        
        ## --------------------------------------------
        
        ## bubble plot...
        if(TRUE) {
          
          heatmap.order = read.csv(paste0(deg.filename,'.heatmap.row_order.csv'))
          heatmap.order = heatmap.order[order(heatmap.order$heatmap.row_order),,drop=F]
          
          dim(heatmap.order) ## 39  4
          dim(data.plot.flt) ## 39 37
          
          
          data.plot.2 = NULL 
          data.plot.2 = data.plot.flt
          
          x = NULL 
          y = NULL 
          x = data.plot.2[,colnames(data.plot.2) %in%
                            c('feat',paste0(my.cohort,'.',tissues,'.p.value'))]
          y = data.plot.2[,colnames(data.plot.2) %in%
                            c('feat',paste0(my.cohort,'.',tissues,'.mean.diff'))]
          x = reshape2::melt(x)
          y = reshape2::melt(y)
          x$variable = as.character(x$variable)
          y$variable = as.character(y$variable)
          
          data.plot.2 = data.frame(
            feat = x$feat,
            Cohort = my.cohort,
            Tissue = gsub(paste0(my.cohort,'.'),'',
                          gsub('.p.value$','',x$variable)),
            p.value = x$value,
            mean.diff = y$value,
            stringsAsFactors = F
            
          )
          
          data.plot.2$direction = NA
          data.plot.2$direction[which(data.plot.2$mean.diff>0)] = 'DC'
          data.plot.2$direction[which(data.plot.2$mean.diff<0)] = 'PD'
          
          # data.plot.2$feat = factor(data.plot.2$feat,
          #                           levels = as.character(rev(heatmap.order$feat)))
          
          data.plot.2$Tissue = factor(data.plot.2$Tissue,
                                      levels = c('all','Liver_L','LN','Lung_L',
                                                 'SoftTissue'))
          
          p1 = NULL 
          p1 = ggplot(data.plot.2, aes(Tissue, feat)) +
            geom_point(aes(size = (-log10(p.value))^2,
                           fill = direction),
                       shape = 21) +
            scale_fill_manual(values = plot.colors) +
            theme_minimal()
          
          p1
          
          pdf(file = paste0(deg.filename,'.bubble.pdf'),
              width = 6, height = 6)
          print(p1)
          dev.off()
          
          data.plot.2$Key = paste0(data.plot.2$Tissue,'!',data.plot.2$Cohort)
          data.plot.2$Key = factor(data.plot.2$Key,
                                   levels = c('all!pd1','Liver_L!pd1','LN!pd1',
                                              'Lung_L!pd1',
                                              'SoftTissue!pd1'))
          
          p2 = NULL 
          p2 = ggplot(data.plot.2, aes(-log10(p.value), feat)) +
            geom_segment( aes(y=feat, yend=feat, x=0, xend=-log10(p.value))) +
            geom_point(aes(fill = direction), shape=21,
                       size = 3) +
            scale_fill_manual(values = plot.colors) +
            theme_minimal() +
            facet_grid(~ Key) +
            geom_vline(xintercept = -log10(0.05), linetype = 'dashed',color = '#CC0000') +
            geom_vline(xintercept = -log10(1), linetype = 'solid',color = '#000000')
          
          
          p2
          
          pdf(file = paste0(deg.filename,'.lollipop.pdf'),
              width = 10, height = 8)
          print(p2)
          dev.off()
          
          
        }
        
      }
      
    }
    
  }
  
  ## --------------------------------------------
  
  ## same organ, ipinivo vs pd1 
  if(TRUE) {
    
    ## --------------------------------------------
    
    ## prep input
    if(TRUE) {
      
      
      data.plot = NULL 
      i = 0
      
      for(my.cohort in c('ipinivo','pd1')) {
        print(my.cohort)
        
        for(my.tissue in tissues) {
          
          print(my.tissue)
          my.df = NULL 
          my.df = deg[[my.cohort]][[my.tissue]]
          
          if(is.null(my.df)) { next }
          if(length(nrow(my.df)) == 0) { next }
          
          row.names(my.df) = my.df$feat
          colnames(my.df) = paste0(my.cohort,'.',my.tissue,'.',colnames(my.df))
          
          i = i + 1
          
          if(i == 1) {
            data.plot = my.df
          } else {
            data.plot = merge(data.plot, my.df,
                              by = 'row.names', all = T)
            row.names(data.plot) = data.plot$Row.names
            data.plot = data.plot[,-1]
          }
          
        }
        
      }
      
      dim(data.plot) #  372 132
      
      data.plot = data.frame(feat = row.names(data.plot),
                             data.plot,
                             stringsAsFactors = F)
      
      write.csv(data.plot,
                file = paste0(output, '.combined.sameorgan_diffcohorts.csv'))
      
      ## --------------------------------------------
      
      data.plot.flt = NULL 
      deg.filename = NULL 
      
      
      data.plot.flt = data.plot[
        (!is.na(data.plot$ipinivo.Brain.p.value) & 
           data.plot$ipinivo.Brain.p.value < rawp) |
          (!is.na(data.plot$ipinivo.LN.p.value) & 
             data.plot$ipinivo.LN.p.value < rawp) |
          (!is.na(data.plot$ipinivo.Liver_L.p.value) & 
             data.plot$ipinivo.Liver_L.p.value < rawp) |
          (!is.na(data.plot$ipinivo.Lung_L.p.value) & 
             data.plot$ipinivo.Lung_L.p.value < rawp) |
          (!is.na(data.plot$ipinivo.SoftTissue.p.value) & 
             data.plot$ipinivo.SoftTissue.p.value < rawp) |
          
          (!is.na(data.plot$pd1.LN.p.value) & 
             data.plot$pd1.LN.p.value < rawp) |
          (!is.na(data.plot$pd1.Liver_L.p.value) & 
             data.plot$pd1.Liver_L.p.value < rawp) |
          (!is.na(data.plot$pd1.Lung_L.p.value) & 
             data.plot$pd1.Lung_L.p.value < rawp) |
          (!is.na(data.plot$pd1.SoftTissue.p.value) & 
             data.plot$pd1.SoftTissue.p.value < rawp)
        
        ,,drop=F]
      
      dim(data.plot.flt) # 14  133
      
      deg.filename = paste0(output, '.combined.sameorgan_diffcohorts.rawp',rawp)
      
      # }
      
      dim(data.plot.flt) ## 76 133
      print(deg.filename)
      
      write.csv(data.plot.flt,
                file = paste0(deg.filename,'.csv'))
    }
    
    ## --------------------------------------------
    
    if(nrow(data.plot.flt) > 0) { 
      
      ## bubble plot...
      if(TRUE) {
        
        
        dim(data.plot.flt) ## 76 133
        
        
        data.plot.2 = NULL 
        data.plot.2 = data.plot.flt
        
        x = NULL 
        y = NULL 
        x = data.plot.2[,colnames(data.plot.2) %in%
                          c('feat',paste0('ipinivo.',tissues,'.p.value'),
                            paste0('pd1.',tissues,'.p.value'))]
        y = data.plot.2[,colnames(data.plot.2) %in%
                          c('feat',paste0('ipinivo.',tissues,'.p.value'),
                            paste0('pd1.',tissues,'.mean.diff'))]
        x = reshape2::melt(x)
        y = reshape2::melt(y)
        x$variable = as.character(x$variable)
        y$variable = as.character(y$variable)
        
        data.plot.2 = data.frame(
          feat = x$feat,
          Cohort = gsub('[.]\\S+$','',x$variable),
          Tissue = gsub('^ipinivo[.]|^pd1[.]','',
                        gsub('.p.value$','',x$variable)),
          p.value = x$value,
          mean.diff = y$value,
          stringsAsFactors = F
          
        )
        
        data.plot.2$direction = NA
        data.plot.2$direction[which(data.plot.2$mean.diff>0)] = 'DC'
        data.plot.2$direction[which(data.plot.2$mean.diff<0)] = 'PD'
        
        data.plot.2 = data.plot.2[!data.plot.2$Tissue %in% c('all','Brain'),]
        
        # data.plot.2$feat = factor(data.plot.2$feat,
        #                           levels = as.character(rev(heatmap.order$feat)))
        
        data.plot.2$Tissue = factor(data.plot.2$Tissue,
                                    levels = c('all','Lung_L','LN','Liver_L',
                                               'Brain','SoftTissue'))
        
        
        p1 = NULL 
        p1 = ggplot(data.plot.2, aes(Cohort, feat)) +
          geom_point(aes(size = (-log10(p.value))^2,
                         fill = direction),
                     shape = 21) +
          scale_fill_manual(values = plot.colors) +
          theme_minimal() +
          facet_grid(~ Tissue)
        
        p1
        
        pdf(file = paste0(deg.filename,'.bubble.pdf'),
            width = 6, height = 6)
        print(p1)
        dev.off()
        
        data.plot.2$Key = paste0(data.plot.2$Tissue,'!',data.plot.2$Cohort)
        data.plot.2$Key = factor(data.plot.2$Key,
                                 levels = c('Lung_L!ipinivo','Lung_L!pd1',
                                            'LN!ipinivo','LN!pd1',
                                            'Liver_L!ipinivo','Liver_L!pd1',
                                            'SoftTissue!ipinivo','SoftTissue!pd1'
                                 ))
        
        p2 = NULL 
        p2 = ggplot(data.plot.2, aes(-log10(p.value), feat)) +
          geom_segment( aes(y=feat, yend=feat, x=0, xend=-log10(p.value))) +
          geom_point(aes(fill = direction), shape=21,
                     size = 3) +
          scale_fill_manual(values = plot.colors) +
          theme_minimal() +
          facet_grid(~ Key) +
          # geom_vline(xintercept = -log10(0.05), linetype = 'dashed',color = '#CC0000') +
          geom_vline(xintercept = -log10(0.01), linetype = 'dashed',color = '#CC0000') +
          geom_vline(xintercept = -log10(1), linetype = 'solid',color = '#000000')
        ## rawp 0.01 for same organ feat diff
        
        p2
        
        pdf(file = paste0(deg.filename,'.lollipop.pdf'),
            width = 12, height = 5)
        print(p2)
        dev.off()
        
      }
      
    }
    
  }
  
}
