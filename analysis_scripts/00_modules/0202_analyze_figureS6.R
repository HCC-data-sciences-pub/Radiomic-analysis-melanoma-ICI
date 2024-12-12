
## DIFF FEAT analysis 
if(TRUE) {
  
  my.cohort = 'braf'
  my.tissue = 'all'
  print(tissues)

  for(my.cohort in cohorts) {
    
    print(my.cohort)
    
    for(my.tissue in tissues) {
      print(my.tissue)
      
      ## --------------------------------------------
      
      my.df = NULL 
      my.df = radiomics[[my.cohort]][[my.tissue]]
      dim(my.df) 
      
      my.df =my.df[!is.na(my.df$Response),,drop=F]
      
      write.csv(my.df,
                file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.csv'))
      
      ## --------------------------------------------
      
      table(my.df$Response)

      if(sum(my.df$Response=='DC') < 5 | 
         sum(my.df$Response=='PD') < 5) {
        next 
      }
      
      ## --------------------------------------------
      
      ## log transform 
      min(my.df[,-1]) ## -2.360627
      max(my.df[,-1]) ## 73085.17
      
      if(min(my.df[,-1]) < 0) {
        my.df[,-1] = my.df[,-1] + -1 * min(my.df[,-1])
      }
      
      min(my.df[,-1]) ## 0
      max(my.df[,-1]) ## 73087.53
      
      my.df[,-1] = log10(my.df[,-1] + 1)
      
      min(my.df[,-1]) ## 0
      max(my.df[,-1]) ## 4.863849
      
      write.csv(my.df,
                file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.csv'))
      
      ## --------------------------------------------
      
      ## visualize ...
      if(TRUE) {
        
        ## distribution?
        if(TRUE) {
          
          data.plot = NULL 
          data.plot = data.frame(patient = row.names(my.df),
                                 my.df,
                                 stringsAsFactors = )
          
          data.plot = reshape2::melt(data.plot)
          colnames(data.plot) = c('patient','Response','radiomics.feat','radiomics.value')
          
          write.csv(data.plot,
                    file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.box.csv'))
          
          p1 = ggplot(data.plot, aes(radiomics.value)) +
            geom_density(aes(color = patient)) +
            theme_pubr(x.text.angle = 90) +
            theme(legend.position = 'none') +
            facet_grid(~ Response)
          
          p1
          
          p2 = ggplot(data.plot, aes(patient, radiomics.value)) +
            geom_boxplot() +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ Response, scales = "free_x", space = 'free_x')
          
          p2
          
          
          p3 = ggplot(data.plot, aes(patient, radiomics.value)) +
            geom_violin() +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ Response, scales = "free_x", space = 'free_x')
          
          p3
          
          pdf(file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.box.pdf'),
              width = 10, height = 10)
          print(p1 + p2 + p3 +
                  plot_layout(ncol = 1))
          dev.off()
          
        }
        
        ## --------------------------------------------
        
        ## PCA 
        if(TRUE) {
          
          rm(data.plot.2, sample.anno)
          
          data.plot.2 = t(my.df[,-1])
          feat.total = nrow(data.plot.2)
          
          
          sample.anno = my.df[,1,drop=F]
          table(sample.anno$Response)

          ## ------------------------------------------------------
          
          data.plot.2 = t(data.plot.2)
          data.plot.2[1:3,1:2]
          head(sample.anno)
          
          ## ------------------------------------------------------
          
          data.plot.2 = merge(data.plot.2,
                              sample.anno, 
                              by = 'row.names')
          row.names(data.plot.2) = data.plot.2[,1]
          data.plot.2 = data.plot.2[,-1]
          
          write.csv(data.plot.2,
                    file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.pca.csv'))
          
          ## ------------------------------------------------------
          
          data.pca = prcomp(data.plot.2[,1:feat.total])
          
          saveRDS(data.pca,
                  file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.pca.rds'))
          
          ## ------------------------------------------------------
          
          p5 = autoplot(data.pca, data = data.plot, colour = 'Response',
                        label = F, label.size = 3) +
            scale_color_manual(values = plot.colors) +
            theme_minimal() +
            ggtitle(paste0(my.cohort,' ',my.tissue))
          
          p5
          
          pdf(file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.pca.pdf'), 
              width = 6, height = 6)
          print(p5 + 
                  plot_layout(ncol = 1))
          dev.off()
          
          
          
        }
        
        ## --------------------------------------------

      }
      
      ## --------------------------------------------
      
      ## remove low expr variables....
      if(TRUE) {
        
        dim(my.df) 
        
        x = NULL 
        x = data.frame(t(my.df[,-1]))
        x$mean = apply(x,1,mean,na.rm=T)
        x$higher.than.mean = apply(x[,-ncol(x)], 1, function(x) sum(x>mean(x,na.rm=T)))
        x = x[x$higher.than.mean >= round(nrow(my.df)* expr.min),,drop=F]
        dim(x)
        
        my.df = my.df[,c('Response',row.names(x)),drop=F]
        dim(my.df)
        
        write.csv(my.df,
                  file = paste0(output,'.',my.cohort,'.',my.tissue,'.data.log10.exLow.csv'))
      }
      
      ## --------------------------------------------
      
      ## DIFF FEAT analysis 
      if(TRUE) {
        
        ## stats 
        if(TRUE) {
          data.stats = NULL 
          
          ## --------------------------------------------
          
          for(i in 2:ncol(my.df)) {
            print(i)
            my.feat = NULL 
            my.feat = colnames(my.df)[i]
            
            my.df.2 = NULL 
            
            my.df.2 = my.df[,c('Response',my.feat)]
            
            my.stats = NULL 
            rm(x,y)
            x = my.df.2[my.df.2$Response == 'DC',my.feat]
            y = my.df.2[my.df.2$Response == 'PD',my.feat]
            my.stats = tidy(wilcox.test(x,y))
            
            data.stats = rbind(data.stats,
                               data.frame(
                                 feat = my.feat,
                                 comp = 'DC_vs_PD',
                                 DC.total = length(x),
                                 PD.total = length(y),
                                 DC.mean = mean(x, na.rm=T),
                                 PD.mean = mean(y, na.rm=T),
                                 my.stats,
                                 stringsAsFactors = F
                               ))
            
          }
          
          data.stats$p.adj = p.adjust(data.stats$p.value, method = 'fdr')
          data.stats$mean.diff = (data.stats$DC.mean - data.stats$PD.mean) 
          
          
          ## --------------------------------------------
          
          
          
          data.stats.sig = NULL 
          
          if(use.fdr == 1) {
            data.stats.sig = data.stats[!is.na(data.stats$p.adj) &
                                          data.stats$p.adj < fdr,,drop=F]
          } else if (use.rawp == 1) {
            data.stats.sig = data.stats[!is.na(data.stats$p.value) &
                                          data.stats$p.value < rawp,,drop=F]
          }
          
          
          # data.stats.sig = data.stats.sig[!is.na(data.stats.sig$mean.diff) &
          #                                   abs(data.stats.sig$mean.diff) >= mean.diff.min,,drop=F]
          print(dim(data.stats.sig))
          
          write.csv(data.stats,
                    file = paste0(output,'.',my.cohort,'.',my.tissue,
                                  '.data.log10.wilcox.csv'))
          
          write.csv(data.stats.sig,
                    file = paste0(output,'.',my.cohort,'.',my.tissue,
                                  '.data.log10.wilcox.sig.csv'))
          
        }
        
        ## --------------------------------------------
        
        ## DIFF FEAT heatmap
        if(TRUE & nrow(data.stats.sig) > 0) {
          
          centered_data = NULL 
          p1 = NULL
          
          centered_data = t(my.df[,data.stats.sig$feat,drop=F])
          centered_data = t(scale(t(centered_data)))
          dim(centered_data)
          
          ## ------------------------------------------------------
          
          if(TRUE) {
            
            ## add annotation 
            sample.anno = NULL 
            sample.anno = my.df[,1,drop=F]
            
            ## sort sample anno same as expression matrix
            sample.anno = sample.anno[row.names(sample.anno) %in% colnames(centered_data),,
                                      drop=F]
            print(all.equal(row.names(sample.anno), colnames(centered_data)))
            
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
                    file = paste0(output,'.',my.cohort,'.',my.tissue,
                                  '.data.log10.wilcox.sig.heatmap.anno.csv'))
          
          write.csv(centered_data,
                    file = paste0(output,'.',my.cohort,'.',my.tissue,
                                  '.data.log10.wilcox.sig.heatmap.data.csv'))
          
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
          
          pdf(file = paste0(output,'.',my.cohort,'.',my.tissue,
                            '.data.log10.wilcox.sig.heatmap.pdf'), 
              width = 10, height = 10)
          print(p1)
          dev.off()
          
          
        }
        
        ## --------------------------------------------
      
      }
      
    }
  }
  
}
