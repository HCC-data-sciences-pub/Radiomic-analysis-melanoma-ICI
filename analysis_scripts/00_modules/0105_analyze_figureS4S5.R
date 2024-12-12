
## Figure S3-S5
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
    
    ## add disease control col
    data.plot$Best.Response.clp = NA
    data.plot$Best.Response.clp[which(data.plot$Best.Response %in% c('PD'))] = 'PD'
    data.plot$Best.Response.clp[which(data.plot$Best.Response %in% c('SD','PR','CR'))] = 'DC'
    
    table(data.plot[,c('Cohort','Best.Response.clp')])
    sum(is.na(data.plot$Best.Response.clp)) ## 0
    
    
    ## --------------------------------------------
    
    write.csv(data.plot,
              file = gsub('.csv','.fixed.csv',clinical.wOrganResponse.file))
    
  }
  
  ## --------------------------------------------
  
  ## draw figures
  if(TRUE) {
    
    ## prep input
    if(TRUE) {
      
      organ.response.tumorsize.cols = NULL 
      organ.response.tumorsize.cols = c(
        paste0(unlist(tissues),'..baseline.'),
        paste0(unlist(tissues),'.2..baseline.'),
        paste0(unlist(tissues),'..best.response.'),
        paste0(unlist(tissues),'.2..best.response.')
      )
      data.plot.2 = NULL 
      data.plot.2 = data.plot[,c('Cohort','Cohort2','Deidentified.name',
                                 'Best.Response', 'Best.Response.clp',
                                 data.plot.0.cols,
                                 names(tissues),
                                 organ.response.tumorsize.cols)]
      
      head(data.plot.2[,organ.response.tumorsize.cols])
      data.plot.2[,organ.response.tumorsize.cols] = 
        apply(data.plot.2[,organ.response.tumorsize.cols],
              2,
              function(x) gsub(' +','',x))
      data.plot.2[,organ.response.tumorsize.cols] = 
        apply(data.plot.2[,organ.response.tumorsize.cols],2,
              function(x) as.numeric(as.character(x)))
      
      head(data.plot.2[,organ.response.tumorsize.cols])
      
      ## --------------------------------------------
      
      ## 1/25/2024 update: swap lesion 01 and 02 if 02 baseline size is bigger
      if(TRUE) {
        
        for(i in 1:nrow(data.plot.2)) {
          # print(i)
          
          for(my.tissue in unlist(tissues)) {
            # print(my.tissue)
            
            x = NULL 
            y = NULL 
            x = data.plot.2[i,paste0(my.tissue,'..baseline.')]
            y = data.plot.2[i,paste0(my.tissue,'.2..baseline.')]
            
            if(is.na(x) & is.na(y)) {
              ## no lesion in this organ site for this patient
            } else if(is.na(x) | is.na(y)) {
              ## only one lesion
            } else if (x >= y) {
              ## then lesion 01 is equal to/bigger than 02, no need to swap
            } else if (x < y) {
              ## need to swap lesion 01 and 02 baseline cols
              data.plot.2[i,paste0(my.tissue,'..baseline.')] = y
              data.plot.2[i,paste0(my.tissue,'.2..baseline.')] = x
              
              ## swap corresponding best response cols too
              z = NULL 
              w = NULL 
              z = data.plot.2[i,paste0(my.tissue,'..best.response.')]
              w = data.plot.2[i,paste0(my.tissue,'.2..best.response.')]
              
              data.plot.2[i,paste0(my.tissue,'..best.response.')] = w
              data.plot.2[i,paste0(my.tissue,'.2..best.response.')] = z
              
              
            } else {
              
              print("other situations???")
              
            }
            
            
          }
          
        }
        
        write.csv(data.plot.2,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.csv'))
      }
      
      ## --------------------------------------------
      
      ## add tumor size change% for each lesion
      ## 2/24 update: add weighted lesion change % per Alex formula....
      for(my.tissue in names(tissues)) {
        print(my.tissue)
        
        my.df = NULL 
        my.df = data.frame(matrix(NA, ncol = 4, nrow = nrow(data.plot.2)))
        colnames(my.df) = c(
          paste0(tissues[[my.tissue]],'..recist.'), ## lesion 1
          paste0(tissues[[my.tissue]],'.2..recist.'), ## lesion 2
          paste0(tissues[[my.tissue]],'.mean..recist.'), ## lesion 1 and 2 avg
          paste0(tissues[[my.tissue]],'.weighted..recist.') ## lesion 1 and 2 weighted
        )
        
        ## compute weights 
        if(TRUE) {
          rm(w1,w2)
          
          w1 = c(NA, rep = nrow(data.plot.2))
          w2 = w1
          
          ## have to loop through every row ...
          for(j in 1:nrow(data.plot.2)) {
            print(j)
            
            my.df.2 = NULL 
            my.df.2 = data.plot.2[j,,drop=F]
            
            if(!is.na(my.df.2[,paste0(tissues[[my.tissue]],'..baseline.')]) &
               !is.na(my.df.2[,paste0(tissues[[my.tissue]],'.2..baseline.')])) {
              
              w1[j] = (my.df.2[,paste0(tissues[[my.tissue]],'..baseline.')]) / 
                (my.df.2[,paste0(tissues[[my.tissue]],'..baseline.')] + 
                   my.df.2[,paste0(tissues[[my.tissue]],'.2..baseline.')])
              w2[j] = (my.df.2[,paste0(tissues[[my.tissue]],'.2..baseline.')]) / 
                (my.df.2[,paste0(tissues[[my.tissue]],'..baseline.')] + 
                   my.df.2[,paste0(tissues[[my.tissue]],'.2..baseline.')])
              
            } else if(!is.na(my.df.2[,paste0(tissues[[my.tissue]],'..baseline.')]) &
                      is.na(my.df.2[,paste0(tissues[[my.tissue]],'.2..baseline.')])) {
              
              w1[j] = 1
              w2[j] = NA
              
            } else if(is.na(my.df.2[,paste0(tissues[[my.tissue]],'..baseline.')]) &
                      is.na(my.df.2[,paste0(tissues[[my.tissue]],'.2..baseline.')])) {
              
              w1[j] = NA
              w2[j] = NA
            } else {
              print("L1 baseline data missing but L2 baseline data present??")
            }
            
          }
          
          
          dim(cbind(w1,w2)) 
          head(cbind(w1,w2))
          
        }
        
        ## --------------------------------------------
        
        ## lesion 1
        rm(x,y,z)
        x = data.plot.2[,paste0(tissues[[my.tissue]],'..baseline.')]
        y = data.plot.2[,paste0(tissues[[my.tissue]],'..best.response.')]
        z = (y-x) / x * 100
        
        data.plot.2[,paste0(tissues[[my.tissue]],'..recist.')] = z
        
        ## lesion 2
        rm(x,y,z)
        x = data.plot.2[,paste0(tissues[[my.tissue]],'.2..baseline.')]
        y = data.plot.2[,paste0(tissues[[my.tissue]],'.2..best.response.')]
        z = (y-x) / x * 100
        
        data.plot.2[,paste0(tissues[[my.tissue]],'.2..recist.')] = z
        
        ## --------------------------------------------
        
        ## lesion 1 and 2 avg
        data.plot.2[,paste0(tissues[[my.tissue]],'.mean..recist.')] = 
          apply(data.plot.2[,c(
            paste0(tissues[[my.tissue]],'..recist.'),
            paste0(tissues[[my.tissue]],'.2..recist.')
          )], 1, mean, na.rm=T)
        
        ## --------------------------------------------
        
        ## lesion 1 and 2 weighted 
        if(TRUE) {
          
          rm(a,b,c)
          a = data.plot.2[,paste0(tissues[[my.tissue]],'..recist.')]
          b = data.plot.2[,paste0(tissues[[my.tissue]],'.2..recist.')]
          
          c = data.frame(
            Patient.identifier = data.plot.2$Patient.identifier,
            Cohort = data.plot.2$Cohort,
            Tissue = my.tissue,
            L1.bs = data.plot.2[,paste0(tissues[[my.tissue]],'..baseline.')],
            L1.br = data.plot.2[,paste0(tissues[[my.tissue]],'..best.response.')],
            a = a,
            w1 = w1,
            axw1 = a * w1,
            L2.bs = data.plot.2[,paste0(tissues[[my.tissue]],'.2..baseline.')],
            L2.br = data.plot.2[,paste0(tissues[[my.tissue]],'.2..best.response.')],
            b = b,
            w2 = w2,
            bxw2 = b * w2, 
            weighted = a * w1 + b * w2, 
            mean = (a+b)/2,
            stringsAsFactors = F)
          # row.names(c) = 1:nrow(c)
          head(c)
          
          
          ## fix it...
          d = NULL 
          d = which(!is.na(c$axw1) & !is.na(c$bxw2))
          if(length(d) > 0 ) { c$weighted[d] = c$axw1[d] + c$bxw2[d] }
          d = NULL 
          d = which(!is.na(c$axw1) & is.na(c$bxw2))
          if(length(d) > 0 ) { c$weighted[d] = c$axw1[d] }
          d = NULL 
          d = which(is.na(c$axw1) & is.na(c$bxw2))
          if(length(d) > 0 ) { c$weighted[d] = NA }
          d = NULL 
          d = which(is.na(c$axw1) & !is.na(c$bxw2))
          if(length(d) > 0 ) { 
            print('something is wrong here!! L2 has data but L1 data missing') 
            
            c$weighted[d] = NA
            
          }
          print(c[d,])
          
          
          
          c$mean = apply(c[,c('a','b')],1, mean, na.rm=T)
          head(c)
          
          ## --------------------------------------------
          
          ## sanity check .... since I computed mean using two approaches here 
          print(all.equal(data.plot.2[,paste0(tissues[[my.tissue]],'.mean..recist.')],
                          c$mean))
          # TRUE; yay, all consistent!!
          
          
          write.csv(c,
                    file = paste0(clinical.file,'.ptByResOrganRECIST',my.tissue,'.csv'))
          
          ## --------------------------------------------
          
          data.plot.2[,paste0(tissues[[my.tissue]],'.weighted..recist.')]  = c$weighted
          
        }
        
        ## --------------------------------------------
        
        rm(x)
        x = data.plot.2[,paste0(tissues[[my.tissue]],'.mean..recist.')]
        x[is.nan(x)] = NA
        
        data.plot.2[,paste0(tissues[[my.tissue]],'.mean..recist.')] = x
        
        ## --------------------------------------------
        
        rm(x)
        x = data.plot.2[,paste0(tissues[[my.tissue]],'.weighted..recist.')]
        x[is.nan(x)] = NA
        
        data.plot.2[,paste0(tissues[[my.tissue]],'.weighted..recist.')] = x
        
        head(data.plot.2[,c(paste0(tissues[[my.tissue]],'.mean..recist.'),
                            paste0(tissues[[my.tissue]],'.weighted..recist.'))])
        
      }
      
      head(data.plot.2)
      
      write.csv(data.plot.2,
                file = paste0(clinical.file,'.ptByResOrganRECIST.csv'))
      
      ## --------------------------------------------
      
      if(TRUE) {
        
        x = NULL 
        x = read.csv(file = paste0(clinical.file,'.ptByResOrganRECIST.csv'),
                     row.names = 1)
        
        dim(x)
        
        ## --------------------------------------------
        
        y = NULL 
        y = x[,c('Cohort','Cohort2','Deidentified.name',
                 'Best.Response', 'Best.Response.clp',
                 'Patient.identifier',
                 paste0(unlist(tissues),'..baseline.'),
                 paste0(unlist(tissues),'.2..baseline.'))]
        
        for(my.tissue in unlist(tissues)) {
          print(my.tissue)
          
          y$new = NA 
          y$new = y[,paste0(my.tissue,'..baseline.')] - 
            y[,paste0(my.tissue,'.2..baseline.')]
          colnames(y)[ncol(y)] = paste0(my.tissue,'.lesion01minus02')
          
        }
        
        write.csv(y,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.lesion0102.csv'))
        
        ## --------------------------------------------
        
        ## direction .....
        z = NULL 
        z =  reshape2::melt(y[,c('Cohort','Cohort2','Deidentified.name',
                                 'Best.Response', 'Best.Response.clp',
                                 'Patient.identifier',
                                 paste0(unlist(tissues),'.lesion01minus02'))])
        z$tissue = gsub('.lesion01minus02','',z$variable)
        z$direction = NA
        z$direction[which(z$value >=0)] = 'equalOrLarger'
        z$direction[which(z$value <0)] = 'smaller'
        z$Key = paste0(z$Cohort,'!',z$Deidentified.name,'!',
                       z$tissue)
        
        write.csv(z,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.lesion0102.dir.csv'))
        
        ## --------------------------------------------
        
        dim(z[z$direction %in% c('smaller'),]) # 48  11
        length(unique(z$Deidentified.name[z$direction %in% c('smaller')])) ## 46 pts 
        
        w = NULL 
        w = data.frame(table(z[,c('Cohort','tissue','direction')]))
        w = reshape2::dcast(Cohort ~ tissue + direction,
                            value.var = 'Freq',
                            data = w,
                            fun.aggregate = sum)
        
        write.csv(w,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.lesion0102.dir.stats.csv'))
        
        ## --------------------------------------------
        
        y = reshape2::melt(y[,c('Cohort','Cohort2','Deidentified.name',
                                'Best.Response', 'Best.Response.clp',
                                'Patient.identifier',
                                paste0(unlist(tissues),'..baseline.'),
                                paste0(unlist(tissues),'.2..baseline.'))])
        y$tissue = gsub('..baseline.|.2..baseline.','',y$variable)
        y$Key = paste0(y$Cohort,'!',y$Deidentified.name,'!',
                       y$tissue)
        
        dim(y) ## 3492    9
        dim(z) ## 1746    9
        sum(duplicated(z$Key)) ## 0
        y = merge(y,z[,c('Key','direction')], by = 'Key')
        dim(y) ## 3492   10
        
        y$lesion = NA
        y$lesion[grep('..baseline.',y$variable)] = 'lesion01'
        y$lesion[grep('.2..baseline.',y$variable)] = 'lesion02'
        y$baseline = y$value
        
        dim(y) ## 3492   12
        y = y[!is.na(y$baseline),]
        dim(y) ## 726  12
        
        write.csv(y,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.lesion0102.dir.mgd.csv'))
        
        ## --------------------------------------------
        
        p.multi = list()
        for(my.cohort in c('pd1','ipinivo','braf')) {
          print(my.cohort)
          
          my.df = NULL 
          my.df = y[y$Cohort==my.cohort,]
          
          p1 = NULL 
          p1 = ggplot(my.df, aes(lesion, baseline)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#C0C0C0', lwd=0.5, fatten=2) +
            geom_point(aes(fill = Best.Response.clp), 
                       size = 2, shape=21) +
            geom_line(aes(group = Deidentified.name)) +
            scale_fill_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(direction ~ tissue)  +
            ggtitle(my.cohort)
          p1
          
          p.multi[[my.cohort]] = p1
          
        }
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.lesion0102.dir.mgd.line.pdf'),
            height = 8, width = 15)
        print(p.multi[['pd1']] + p.multi[['ipinivo']] + p.multi[['braf']] +
                plot_layout(nrow = 1))
        dev.off()
        
        
        
      }
      
      
    }
    
    ## --------------------------------------------
    
    ## boxplot - organ tumor size change % PD vs DC
    if(TRUE) {
      
      ## plotting - Fig. S3
      if(TRUE) {
        
        table(data.plot.2[,c('Cohort','Best.Response.clp')])
        
        data.plot.3 = NULL 
        data.plot.3 = data.plot.2[,c('Cohort','Deidentified.name',
                                     'Best.Response','Best.Response.clp',
                                     paste0(unlist(tissues),'.weighted..recist.'))]
        data.plot.3 = reshape2::melt(data.plot.3)
        data.plot.3$tissue = gsub('.weighted..recist.','',data.plot.3$variable)
        data.plot.3$.weighted..recist. = data.plot.3$value
        
        data.plot.3$Best.Response.clp = factor(
          data.plot.3$Best.Response.clp,levels = c('PD','DC','NA')
        )
        
        data.plot.3$Best.Response = factor(
          data.plot.3$Best.Response,levels = c('PD','SD','PR','CR','NA')
        )
        
        write.csv(data.plot.3,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.01.csv'))
        
        ## --------------------------------------------
        
        ## add organ response 
        x = NULL
        x = data.plot.2[,c('Cohort','Deidentified.name',
                           names(tissues))]
        x = reshape2::melt(x, id.vars = c('Cohort','Deidentified.name'))
        x$variable = as.character(x$variable)
        x$tissue = NA
        for(my.tissue in names(tissues)) {
          print(my.tissue)
          x$tissue[which(x$variable == my.tissue)] = tissues[[my.tissue]]
        }
        
        x$Organ.response = x$value
        table(x$Organ.response)
        
        x$Organ.response[x$Organ.response=='NE'] = NA
        x$Key = paste0(x$Cohort,'!',x$Deidentified.name,'!',x$tissue)
        
        data.plot.3$Key = paste0(data.plot.3$Cohort,'!',
                                 data.plot.3$Deidentified.name,'!',
                                 data.plot.3$tissue)
        sum(duplicated(x$Key))
        dim(data.plot.3) # 1452   24
        data.plot.3 = merge(data.plot.3, x[,c('Key','Organ.response')],
                            by = 'Key')
        dim(data.plot.3) ## 1452   25
        
        data.plot.3$Organ.response.clp = NA
        data.plot.3$Organ.response.clp[which(data.plot.3$Organ.response %in%
                                               c('PD'))] = 'PD'
        data.plot.3$Organ.response.clp[which(data.plot.3$Organ.response %in%
                                               c('SD','PR','CR'))] = 'DC'
        table(data.plot.3$Organ.response.clp)
        sum(is.na(data.plot.3$Organ.response.clp))
        
        write.csv(data.plot.3,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.01.csv'))
        
        ## --------------------------------------------
        
        p.multi = list()
        p.multi.2 = list()
        p.multi.3 = list()
        p.multi.4 = list()
        y.min = min(data.plot.3$.weighted..recist., na.rm=T)
        y.max = max(data.plot.3$.weighted..recist., na.rm=T)
        y.min
        y.max
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          
          print(my.cohort)
          my.df = NULL 
          
          if(my.cohort == 'ici') {
            my.df = data.plot.3[data.plot.3$Cohort %in% c('pd1','ipinivo'),]
          } else {
            my.df = data.plot.3[data.plot.3$Cohort == my.cohort,]
          }
          
          print(dim(my.df))
          print(length(unique(my.df$Deidentified.name)))
          
          ## --------------------------------------------
          
          p1 = NULL 
          p1 = ggplot(my.df, aes(Best.Response.clp,
                                 .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(fill = Organ.response.clp), 
                        width = 0.2, height = 0,
                        size = 2, shape=21) +
            scale_fill_manual(values = plot.colors) +
            theme_pubr() +
            facet_grid(~ tissue) +
            geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed',
                       color = '#C0C0C0') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('PD','DC'))) +
            ggtitle(my.cohort) +
            ylim(y.min, y.max +50)
          p1
          
          p.multi[[my.cohort]] = p1
          
          ## --------------------------------------------
          
          p1.1 = NULL 
          p1.1 = ggplot(my.df, aes(tissue,
                                   .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(fill = Organ.response.clp), 
                        width = 0.2, height = 0,
                        size = 2, shape=21) +
            scale_fill_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ Best.Response.clp) +
            geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed',
                       color = '#C0C0C0') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('PD','DC'))) +
            ggtitle(my.cohort) +
            ylim(y.min, y.max +50)
          p1.1
          
          p.multi.4[[my.cohort]] = p1.1
          
          ## --------------------------------------------
          
          p2 = NULL 
          p2 = ggplot(my.df, aes(x = .weighted..recist., 
                                 y = Best.Response.clp,
                                 fill = Best.Response)) +
            geom_density_ridges(
              jittered_points = TRUE,
              position = position_points_jitter(width = 0.05, height = 0),
              point_shape = '|',
              point_size = 3, point_alpha = 1, alpha = 0.7
            ) +
            scale_fill_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ tissue) +
            ggtitle(my.cohort) +
            xlim(y.min-150, y.max +50)
          
          p2
          
          p.multi.2[[my.cohort]] = p2
          
          ## --------------------------------------------
          
          p2.1 = NULL 
          p2.1 = ggplot(my.df, aes(x = .weighted..recist., 
                                   y = tissue,
                                   fill = Best.Response.clp)) +
            geom_density_ridges(
              jittered_points = TRUE,
              position = position_points_jitter(width = 0.05, height = 0),
              point_shape = '|',
              point_size = 3, point_alpha = 1, alpha = 0.7
            ) +
            scale_fill_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ Best.Response.clp) +
            ggtitle(my.cohort) +
            xlim(y.min-150, y.max +50)
          
          p2.1
          
          p.multi.3[[my.cohort]] = p2.1
          
        }
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.01.pdf'),
            height = 10, width = 12)
        print(p.multi[['ici']] + p.multi[['pd1']] + 
                p.multi[['ipinivo']] + p.multi[['braf']] +
                plot_layout(ncol = 2))
        dev.off()
        
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.01.2.pdf'),
            height = 10, width = 12)
        print(p.multi.4[['ici']] + p.multi.4[['pd1']] + 
                p.multi.4[['ipinivo']] + p.multi.4[['braf']] +
                plot_layout(ncol = 2))
        dev.off()
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.ridge.01.pdf'),
            height = 10, width = 8)
        print(p.multi.2[['ici']] + p.multi.2[['pd1']] + 
                p.multi.2[['ipinivo']] + p.multi.2[['braf']] +
                plot_layout(ncol = 1))
        dev.off()
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.ridge.01.2.pdf'),
            height = 8, width = 10)
        print(p.multi.3[['ici']] + p.multi.3[['pd1']] + 
                p.multi.3[['ipinivo']] + p.multi.3[['braf']] +
                plot_layout(ncol = 2))
        dev.off()
        
        ## --------------------------------------------
        
        ## inter organ heterogeneity 
        if(TRUE) {
          
          dim(data.plot.3)
          
          ## add mixed response/uniform PD/uniform DC group from analysis above
          data.plot.4 = NULL 
          data.plot.4 = read.csv(
            paste0(clinical.file,'.ptByResponseOrgan.heatmap.csv'),row.names = 1
          )
          
          dim(data.plot.3) ## 1746   11
          sum(duplicated(data.plot.4$Deidentified.name)) ## 0
          data.plot.4 = merge(data.plot.3,
                              data.plot.4[,c('Deidentified.name',
                                             'Response.organ.mixed')],
                              by = 'Deidentified.name')
          dim(data.plot.4) ## 1746   12
          
          write.csv(data.plot.4,
                    file = paste0(clinical.file,'.ptByResOrganRECIST.line.01.csv'))
          
          ## --------------------------------------------
          
          ## lesions within patient
          if(TRUE) {
            
            p.multi = list()
            p.multi.2 = list()
            p.multi.3 = list()
            p.multi.4 = list()
            y.min = min(data.plot.4$.weighted..recist., na.rm=T)
            y.max = max(data.plot.4$.weighted..recist., na.rm=T)
            y.min
            y.max
            for(my.cohort in c('ici','pd1','ipinivo','braf')) {
              
              print(my.cohort)
              my.df = NULL 
              
              if(my.cohort == 'ici') {
                my.df = data.plot.4[data.plot.4$Cohort %in% c('pd1','ipinivo'),]
              } else {
                my.df = data.plot.4[data.plot.4$Cohort == my.cohort,]
              }
              
              my.df$Cohort.Deidentified.name = paste0(my.df$Cohort,'!',my.df$Deidentified.name)
              
              print(dim(my.df))
              print(length(unique(my.df$Deidentified.name))) ## 242
              print(length(unique( my.df$Cohort.Deidentified.name))) ## 242
              
              
              write.csv(my.df,
                        file = paste0(clinical.file,
                                      '.ptByResOrganRECIST.',
                                      my.cohort,'.line.01.csv')
              )
              
              ## --------------------------------------------
              
              ## sort organs by median tumor size % in PD group
              plot.order = NULL 
              plot.order = my.df[my.df$Best.Response.clp=='PD',]
              plot.order = aggregate(.weighted..recist. ~ tissue,
                                     data = plot.order,
                                     FUN = median,
                                     na.rm=T)
              plot.order$tissue = as.character(plot.order$tissue)
              plot.order = plot.order[order(plot.order$.weighted..recist.),]
              dim(plot.order)
              
              write.csv(plot.order,
                        file = paste0(clinical.file,
                                      '.ptByResOrganRECIST.',
                                      my.cohort,'.line.01.median.csv')
              )
              
              ## --------------------------------------------
              
              ## only plot patients with 2+ sites ....
              x = NULL 
              x = my.df
              x$Key = paste0(x$Cohort,'!',x$Deidentified.name)
              x = x[!is.na(x$.weighted..recist.),,drop=F]
              x = data.frame(table(x$Key))
              dim(x) ## 221   2
              x = x[x$Freq>=2,,drop=F]
              x$Var1 = as.character(x$Var1)
              dim(my.df) ## 1452   12
              my.df = my.df[paste0(my.df$Cohort,'!',my.df$Deidentified.name) %in%
                              x$Var1,,drop=F]
              dim(my.df) ## 666  12
              
              write.csv(my.df,
                        file = paste0(clinical.file,
                                      '.ptByResOrganRECIST.',
                                      my.cohort,'.line.01.2.csv')
              )
              
              ## --------------------------------------------
              
              ## sort organs by median tumor size % in PD group
              plot.order = NULL 
              plot.order = my.df[my.df$Best.Response.clp=='PD',]
              plot.order = aggregate(.weighted..recist. ~ tissue,
                                     data = plot.order,
                                     FUN = median,
                                     na.rm=T)
              plot.order$tissue = as.character(plot.order$tissue)
              plot.order = plot.order[order(plot.order$.weighted..recist.),]
              
              write.csv(plot.order,
                        file = paste0(clinical.file,
                                      '.ptByResOrganRECIST.',
                                      my.cohort,'.line.01.2.median.csv')
              )
              
              ## --------------------------------------------
              
              ## calculate sd within pt...
              if(TRUE) {
                
                my.comps = NULL 
                my.comps = list(c('mixedResponse','uniformProgression'),
                                c('mixedResponse',
                                  'uniformDiseaseControl'),
                                c('uniformProgression',
                                  'uniformDiseaseControl'))
                
                if(my.cohort == 'braf') {
                  my.comps = list(c('mixedResponse',
                                    'uniformDiseaseControl'))
                }
                
                ## --------------------------------------------
                
                
                x = NULL 
                x = my.df[!is.na(my.df$.weighted..recist.),]
                x = list(
                  higherthan100 = x[x$.weighted..recist. > 100,,drop=F],
                  within100 = x[x$.weighted..recist. <=100,,drop=F]
                )
                if(nrow(x$higherthan100) > 0) {
                  x$higherthan100$.weighted..recist.
                  x$higherthan100$.weighted..recist. = 100
                  x$higherthan100$.weighted..recist.
                }
                
                
                my.df.sd = NULL 
                my.df.sd = rbind(x$higherthan100,x$within100)
                my.df.sd = aggregate(.weighted..recist. ~ Deidentified.name,
                                     data = my.df.sd,
                                     FUN = sd,
                                     na.rm=T)
                my.df.sd$Deidentified.name = as.character(my.df.sd$Deidentified.name)
                my.df.sd$.weighted..recist..sd = my.df.sd$.weighted..recist.
                my.df.sd = merge(my.df.sd, 
                                 unique(my.df[,c('Deidentified.name',
                                                 'Response.organ.mixed')]),
                                 by = 'Deidentified.name')
                dim(my.df.sd) ## 111 4
                my.df.sd$Response.organ.mixed = factor(
                  my.df.sd$Response.organ.mixed,
                  levels = c('uniformProgression',
                             'mixedResponse',
                             'uniformDiseaseControl')
                )
                
                
                p0 = NULL 
                p0 = ggplot(my.df.sd, aes(Response.organ.mixed, 
                                          .weighted..recist..sd)) +
                  geom_boxplot(width = 0.6,outlier.shape = NA) +
                  geom_jitter(aes(color = Response.organ.mixed),
                              width = 0.1, height = 0) +
                  scale_color_manual(values = plot.colors) +
                  theme_pubr(x.text.angle = 90) +
                  theme(legend.position = 'none') +
                  stat_compare_means(method = 'wilcox.test',
                                     comparisons = my.comps) +
                  ggtitle(paste0(my.cohort,": tumor size pos% 100+ reset to 100"),
                          subtitle = paste0(nrow(my.df.sd),' pts'))
                
                
                print(p0) 
                
                p.multi.4[[my.cohort]] = p0
                
                write.csv(my.df.sd,
                          file = paste0(clinical.file, '.ptByResOrganRECISTsd.',
                                        my.cohort,'.box.01.2_4.csv'))
                
              }
              
              ## --------------------------------------------
              
              ## visualize individual pt...
              if(TRUE) {
                my.df$tissue = factor(my.df$tissue,
                                      levels = rev(plot.order$tissue))
                
                p1.1 = NULL 
                p1.1 = ggplot(my.df[my.df$Response.organ.mixed=='mixedResponse' &
                                      !is.na(my.df$.weighted..recist.),,drop=F],
                              aes(tissue, .weighted..recist.)) +
                  geom_boxplot(width = 0.6, outlier.shape = NA,
                               color = '#000000', lwd=0.5, fatten=2) +
                  geom_point(aes(fill = Organ.response.clp), 
                             size = 2, shape=21) +
                  geom_line(aes(group = Deidentified.name),
                            color = '#C0C0C0') +
                  scale_fill_manual(values = plot.colors) +
                  theme_pubr(x.text.angle = 90) +
                  facet_wrap( ~ Deidentified.name, ncol=6) +
                  # geom_hline(yintercept = c(0, -30,-70), 
                  #   linetype = 'dashed',
                  #            color = '#C0C0C0') +
                  stat_compare_means(method = 'wilcox.test',
                                     comparisons = list(c('PD','DC'))) +
                  ggtitle(my.cohort) +
                  ylim(y.min, y.max +50)
                p1.1
                
                p1.2 = NULL 
                p1.2 = ggplot(my.df[my.df$Response.organ.mixed=='uniformDiseaseControl' &
                                      !is.na(my.df$.weighted..recist.),,drop=F],
                              aes(tissue, .weighted..recist.)) +
                  geom_boxplot(width = 0.6, outlier.shape = NA,
                               color = '#000000', lwd=0.5, fatten=2) +
                  geom_point(aes(fill = Organ.response.clp), 
                             size = 2, shape=21) +
                  geom_line(aes(group = Deidentified.name),
                            color = '#C0C0C0') +
                  scale_fill_manual(values = plot.colors) +
                  theme_pubr(x.text.angle = 90) +
                  facet_wrap( ~ Deidentified.name, ncol=6) +
                  # geom_hline(yintercept = c(0, -30,-70), 
                  #   linetype = 'dashed',
                  #            color = '#C0C0C0') +
                  stat_compare_means(method = 'wilcox.test',
                                     comparisons = list(c('PD','DC'))) +
                  ggtitle(my.cohort) +
                  ylim(y.min, y.max +50)
                p1.2
                
                ## braf has only mixedResponse and uniformDiseaseControl, no uniformPD group
                if(my.cohort != 'braf') {
                  p1.3 = NULL 
                  p1.3 = ggplot(my.df[my.df$Response.organ.mixed=='uniformProgression' &
                                        !is.na(my.df$.weighted..recist.),,drop=F],
                                aes(tissue, .weighted..recist.)) +
                    geom_boxplot(width = 0.6, outlier.shape = NA,
                                 color = '#000000', lwd=0.5, fatten=2) +
                    geom_point(aes(fill = Organ.response.clp), 
                               size = 2, shape=21) +
                    geom_line(aes(group = Deidentified.name),
                              color = '#C0C0C0') +
                    scale_fill_manual(values = plot.colors) +
                    theme_pubr(x.text.angle = 90) +
                    facet_wrap( ~ Deidentified.name, ncol=6) +
                    # geom_hline(yintercept = c(0, -30,-70), 
                    #   linetype = 'dashed',
                    #            color = '#C0C0C0') +
                    stat_compare_means(method = 'wilcox.test',
                                       comparisons = list(c('PD','DC'))) +
                    ggtitle(my.cohort) +
                    ylim(y.min, y.max +50)
                  p1.3
                  
                }
                
                pdf(file = paste0(clinical.file, '.ptByResOrganRECIST.',
                                  my.cohort,'.line.01.2.pdf'),
                    width = 13, height = 26)
                if(my.cohort == 'braf') {
                  print(p1.1 + p1.2 + 
                          plot_layout(ncol = 1))
                } else {
                  print(p1.1 + p1.2 + p1.3 +
                          plot_layout(ncol = 1))
                }
                dev.off()
                
                
                
              }
              
            }
            
            pdf(file = paste0(clinical.file, '.ptByResOrganRECISTsd.box.01.2_1.pdf'),
                height = 12, width = 10)
            print(p.multi[['ici']] + p.multi[['pd1']] + 
                    p.multi[['ipinivo']] + p.multi[['braf']] +
                    plot_layout(ncol = 2))
            dev.off()
            
            
            pdf(file = paste0(clinical.file, '.ptByResOrganRECISTsd.box.01.2_2.pdf'),
                height = 12, width = 10)
            print(p.multi.2[['ici']] + p.multi.2[['pd1']] + 
                    p.multi.2[['ipinivo']] + p.multi.2[['braf']] +
                    plot_layout(ncol = 2))
            dev.off()
            
            pdf(file = paste0(clinical.file, '.ptByResOrganRECISTsd.box.01.2_3.pdf'),
                height = 12, width = 10)
            print(p.multi.3[['ici']] + p.multi.3[['pd1']] + 
                    p.multi.3[['ipinivo']] + p.multi.3[['braf']] +
                    plot_layout(ncol = 2))
            dev.off()
            
            pdf(file = paste0(clinical.file, '.ptByResOrganRECISTsd.box.01.2_4.pdf'),
                height = 12, width = 10)
            print(p.multi.4[['ici']] + p.multi.4[['pd1']] + 
                    p.multi.4[['ipinivo']] + p.multi.4[['braf']] +
                    plot_layout(ncol = 2))
            dev.off()
            
          }
          
        }
        
      }
      
      ## --------------------------------------------
      
      ## plotting - Fig. S4
      if(TRUE) {
        
        table(data.plot.2[,c('Cohort','Best.Response.clp')])
        
        data.plot.3 = NULL 
        data.plot.3 = data.plot.2[,c('Cohort','Deidentified.name',
                                     'Best.Response','Best.Response.clp',
                                     paste0(unlist(tissues),'..recist.'),
                                     paste0(unlist(tissues),'.2..recist.'))]
        data.plot.3 = reshape2::melt(data.plot.3)
        data.plot.3$tissue = gsub('..recist.|.2..recist.','',data.plot.3$variable)
        data.plot.3$recist = data.plot.3$value
        data.plot.3$lesion = 'NA'
        data.plot.3$lesion[which(data.plot.3$variable %in%
                                   paste0(data.plot.3$tissue,'..recist.'))] = 'lesion01'
        data.plot.3$lesion[which(data.plot.3$variable %in%
                                   paste0(data.plot.3$tissue,'.2..recist.'))] = 'lesion02'
        sum(is.na(data.plot.3$lesion)) # 0
        
        data.plot.3$Cohort.Deidentified.name.tissue = paste0(
          data.plot.3$Cohort,'!',
          data.plot.3$Deidentified.name,'!',data.plot.3$tissue
        )
        
        data.plot.3$Best.Response.clp = factor(
          data.plot.3$Best.Response.clp,levels = c('PD','DC','NA')
        )
        
        data.plot.3$Best.Response = factor(
          data.plot.3$Best.Response,levels = c('PD','SD','PR','CR','NA')
        )
        
        write.csv(data.plot.3,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.csv'))
        
        ## select paired lesions only
        if(TRUE) {
          
          data.plot.3.paired = NULL 
          data.plot.3.paired = data.frame(table(data.plot.3$Cohort.Deidentified.name.tissue[
            !is.na(data.plot.3$recist)
          ]))
          table(data.plot.3.paired$Freq)
          # 1   2 
          # 258 225 
          
          dim(data.plot.3.paired) ## 483   2
          data.plot.3.paired = data.plot.3.paired[data.plot.3.paired$Freq==2,]
          dim(data.plot.3.paired) ## 225   2
          
          data.plot.3.paired = data.plot.3[data.plot.3$Cohort.Deidentified.name.tissue %in%
                                             as.character(data.plot.3.paired$Var1),]
          dim(data.plot.3) # 3492   10
          dim(data.plot.3.paired) ## 450  10
          
          write.csv(data.plot.3.paired,
                    file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.paired.csv'))
          
          ## how many pts and how many tissues have paired data?
          x = NULL 
          x = unique(data.plot.3.paired[,c('Cohort','Cohort.Deidentified.name.tissue',
                                           'tissue')])
          dim(x) ## 225   3
          x = data.frame(table(x[,c('Cohort','tissue')]))
          x = reshape2::dcast(Cohort ~ tissue, 
                              value.var = 'Freq',
                              data = x,
                              fun.aggregate = sum)
          x
          
          write.csv(x,
                    file = paste0(clinical.file,
                                  '.ptByResOrganRECIST.box.02.paired.stats.csv'))
          
        }
        
        ## --------------------------------------------
        
        p.multi = list()
        y.min = min(data.plot.3$recist, na.rm=T)
        y.max = max(data.plot.3$recist, na.rm=T)
        y.min
        y.max
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          
          print(my.cohort)
          my.df = NULL 
          
          if(my.cohort == 'ici') {
            my.df = data.plot.3[data.plot.3$Cohort %in% c('pd1','ipinivo'),]
          } else {
            my.df = data.plot.3[data.plot.3$Cohort == my.cohort,]
          }
          
          print(dim(my.df))
          print(length(unique(my.df$Cohort.Deidentified.name.tissue)))
          
          ## --------------------------------------------
          
          p1 = NULL 
          p1 = ggplot(my.df, aes(lesion,recist)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#C0C0C0', lwd=0.5, fatten=2) +
            geom_point(aes(color = Best.Response),
                       size = 1) +
            geom_line(aes(group = Cohort.Deidentified.name.tissue),
                      size = 0.5) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('lesion01','lesion02')),
                               paired = T) +
            ggtitle(my.cohort) +
            ylim(y.min, y.max +50)
          p1
          
          p.multi[[my.cohort]] = p1
          
        }
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.pdf'),
            height = 10, width = 12)
        print(p.multi[['ici']] + p.multi[['pd1']] + 
                p.multi[['ipinivo']] + p.multi[['braf']] +
                plot_layout(ncol = 2))
        dev.off()
        
        ## --------------------------------------------
        
        p.multi = list()
        p.multi.2 = list()
        y.min = min(data.plot.3.paired$recist, na.rm=T)
        y.max = max(data.plot.3.paired$recist, na.rm=T)
        y.min
        y.max
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          
          print(my.cohort)
          my.df = NULL 
          
          if(my.cohort == 'ici') {
            my.df = data.plot.3.paired[data.plot.3.paired$Cohort %in% c('pd1','ipinivo'),]
          } else {
            my.df = data.plot.3.paired[data.plot.3.paired$Cohort == my.cohort,]
          }
          
          print(dim(my.df))
          print(length(unique(my.df$Cohort.Deidentified.name.tissue)))
          
          ## --------------------------------------------
          
          p1 = NULL 
          p1 = ggplot(my.df, aes(lesion,recist)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#C0C0C0', lwd=0.5, fatten=2) +
            geom_point(aes(color = Best.Response),
                       size = 1) +
            geom_line(aes(group = Cohort.Deidentified.name.tissue),
                      size = 0.5) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('lesion01','lesion02')),
                               paired = T) +
            ggtitle(my.cohort) +
            ylim(y.min, y.max +50)
          p1
          
          p.multi[[my.cohort]] = p1
          
          ## --------------------------------------------
          
          p1.1 = NULL 
          p1.1 = ggplot(my.df, aes(lesion,recist)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#C0C0C0', lwd=0.5, fatten=2) +
            geom_point(aes(color = Best.Response),
                       size = 1) +
            geom_line(aes(group = Cohort.Deidentified.name.tissue),
                      size = 0.5) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(Best.Response.clp ~  tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('lesion01','lesion02')),
                               paired = T) +
            ggtitle(my.cohort) +
            ylim(y.min, y.max +50)
          p1.1
          
          p.multi.2[[my.cohort]] = p1.1
          
        }
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.paired.pdf'),
            height = 10, width = 12)
        print(p.multi[['ici']] + p.multi[['pd1']] + 
                p.multi[['ipinivo']] + p.multi[['braf']] +
                plot_layout(ncol = 2))
        dev.off()
        
        pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.paired.2.pdf'),
            height = 10, width = 12)
        print(p.multi.2[['ici']] + p.multi.2[['pd1']] + 
                p.multi.2[['ipinivo']] + p.multi.2[['braf']] +
                plot_layout(ncol = 2))
        dev.off()
        
        ## --------------------------------------------
        
        ## I should use the absolute diff between lesion 01 and lesion 02 for statistics 
        if(TRUE) {
          
          p.multi = list()
          p.multi.2 = list()
          y.min = min(data.plot.3.paired$recist, na.rm=T)
          y.max = max(data.plot.3.paired$recist, na.rm=T)
          y.min
          y.max
          for(my.cohort in c('ici','pd1','ipinivo','braf')) {
            
            print(my.cohort)
            my.df = NULL 
            
            if(my.cohort == 'ici') {
              my.df = data.plot.3.paired[data.plot.3.paired$Cohort %in% c('pd1','ipinivo'),]
            } else {
              my.df = data.plot.3.paired[data.plot.3.paired$Cohort == my.cohort,]
            }
            
            print(dim(my.df))
            print(length(unique(my.df$Cohort.Deidentified.name.tissue)))
            
            ## --------------------------------------------
            
            my.df$Cohort.Deidentified.name = paste0(
              my.df$Cohort,'!',my.df$Deidentified.name
            )
            
            ## --------------------------------------------
            
            ## calculate the absolute diff between lesion 01 and lesion 02
            my.df.2 = NULL 
            my.df.2 = reshape2::dcast(Cohort.Deidentified.name + tissue ~ lesion ,
                                      value.var = 'recist',
                                      data = my.df,
                                      fun.aggregate = sum,
                                      na.rm=T)
            sum(is.na(my.df.2$lesion01) | my.df.2$lesion01 == 0) ## 0
            sum(is.na(my.df.2$lesion02) | my.df.2$lesion02 == 0) ## 0
            ## all sites have 2 lesions, as expected (since I only selected paired ones above)
            
            my.df.2$lesion0102.abs.diff = abs(my.df.2$lesion01 - my.df.2$lesion02)
            # my.df.2$lesion0102.abs.diff = my.df.2$lesion0102.abs.diff / 
            #   (    (abs(my.df.2$lesion01) + abs(my.df.2$lesion02)) /  2  )
            
            head(my.df.2)
            
            my.df.2$Cohort.Deidentified.name.tissue = paste0(my.df.2$Cohort,'!',
                                                             my.df.2$Deidentified.name,'!',
                                                             my.df.2$tissue)
            
            dim(my.df.2) ##  57  6
            my.df.2 = merge(my.df.2, unique(my.df[,c('Cohort','Deidentified.name',
                                                     'Best.Response',
                                                     'Best.Response.clp',
                                                     'Cohort.Deidentified.name')]),
                            by = 'Cohort.Deidentified.name')
            dim(my.df.2) ## 57 10
            
            write.csv(my.df.2,
                      file = paste0(clinical.file,'.ptByResOrganRECIST',
                                    '.',my.cohort,'.box.02.paired.absdiff.csv')
            )
            
            ## --------------------------------------------
            
            shapiro.test(log10( my.df.2$lesion0102.abs.diff + 1))
            
            p0 = NULL 
            p0 = ggplot(my.df.2, aes(log10( lesion0102.abs.diff + 1))) +
              geom_density()
            
            p0
            
            ## --------------------------------------------
            
            p1 = NULL 
            p1 = ggplot(my.df.2, aes(Best.Response.clp,
                                     log10(lesion0102.abs.diff + 1))) +
              geom_boxplot(width = 0.6, outlier.shape = NA,
                           color = '#000000', lwd=0.5, fatten=2) +
              geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                          size = 2) +
              scale_color_manual(values = plot.colors) +
              theme_pubr(x.text.angle = 90) +
              facet_grid(~ tissue) +
              # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
              stat_compare_means(method = 'wilcox.test',
                                 comparisons = list(c('PD','DC')),
                                 paired = F) +
              ggtitle(my.cohort)
            p1
            
            p.multi[[my.cohort]] = p1
            
            ## --------------------------------------------
            
            p2 = NULL 
            p2 = ggplot(my.df.2, aes(x = log10(lesion0102.abs.diff + 1), 
                                     y = tissue,
                                     fill = Best.Response.clp)) +
              geom_density_ridges(
                jittered_points = TRUE,
                position = position_points_jitter(width = 0.05, height = 0),
                point_shape = '|',
                point_size = 3, point_alpha = 1, alpha = 0.7
              ) +
              scale_fill_manual(values = plot.colors) +
              theme_pubr(x.text.angle = 90) +
              ggtitle(my.cohort) 
            
            p2
            
            p1 + p2
            
            p.multi.2[[my.cohort]] = p2
            
          }
          
          
          pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.paired.absdiff.pdf'),
              height = 10, width = 12)
          print(p.multi[['ici']] + p.multi[['pd1']] + 
                  p.multi[['ipinivo']] + p.multi[['braf']] +
                  plot_layout(ncol = 2))
          dev.off()
          
          pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.ridge.02.paired.absdiff.pdf'),
              height = 10, width = 12)
          print(p.multi.2[['ici']]  + p.multi.2[['pd1']] + 
                  p.multi.2[['ipinivo']] + p.multi.2[['braf']] +
                  plot_layout(ncol = 2))
          dev.off()
          
        }
        
      }
      
      ## --------------------------------------------
      
      ## plotting - Fig. S5
      if(TRUE) {
        
        ## ICI cohort only: baseline NLR, eosinophil count and line of immunotherapy.
        
        table(data.plot.2[,c('Cohort','Best.Response.clp')])
        
        data.plot.3 = NULL 
        data.plot.3 = data.plot.2[,c('Deidentified.name',
                                     'Best.Response','Best.Response.clp',
                                     data.plot.0.cols,
                                     paste0(unlist(tissues),'.weighted..recist.'))]
        dim(data.plot.3) ## 291  24
        
        ## select ici cohorts only
        data.plot.3 = data.plot.3[data.plot.3$Cohort %in% c('pd1','ipinivo'),]
        dim(data.plot.3) ## 242  24
        
        print(data.plot.0.cols)
        
        ## add clinical groups .......
        if(TRUE) {
          
          ## add eosino groups 
          eosino.pre.ici.median = NULL 
          eosino.pre.ici.median = median(data.plot.3$X.eosino.pre, na.rm=T)
          print(eosino.pre.ici.median) ## 1.9
          
          data.plot.3$X.eosino.pre.high = NA
          data.plot.3$X.eosino.pre.high[which(data.plot.3$X.eosino.pre > 
                                                eosino.pre.ici.median)] = 'Y'
          data.plot.3$X.eosino.pre.high[which(data.plot.3$X.eosino.pre <= 
                                                eosino.pre.ici.median)] = 'N'
          table(data.plot.3$X.eosino.pre.high)
          # N   Y 
          # 120 112 
          print(sum(is.na(data.plot.3$X.eosino.pre.high)))
          ## 10 
          
          ## add albumin groups 
          albumin.pre.ici.median = NULL 
          albumin.pre.ici.median = median(data.plot.3$Albumin.pre, na.rm=T)
          print(albumin.pre.ici.median) ## 3.9
          
          data.plot.3$Albumin.pre.high = NA
          data.plot.3$Albumin.pre.high[which(data.plot.3$Albumin.pre > 
                                               albumin.pre.ici.median)] = 'Y'
          data.plot.3$Albumin.pre.high[which(data.plot.3$Albumin.pre <= 
                                               albumin.pre.ici.median)] = 'N'
          table(data.plot.3$Albumin.pre.high)
          # N   Y 
          # 131 106
          print(sum(is.na(data.plot.3$Albumin.pre.high)))
          ## 5 
          
          ## add line of imtx groups 
          class(data.plot.3$Line.of.immunotherapy)
          data.plot.3$Line.of.immunotherapy.clp = NA
          data.plot.3$Line.of.immunotherapy.clp[
            which(data.plot.3$Line.of.immunotherapy %in% c(1))
          ] = 'FirstLine'
          data.plot.3$Line.of.immunotherapy.clp[
            which(data.plot.3$Line.of.immunotherapy %in% c(2,3,4))
          ] = 'TwoOrMoreLines'
          
          table(data.plot.3$Line.of.immunotherapy.clp)
          # FirstLine TwoOrMoreLines 
          # 199             43 
          print(sum(is.na(data.plot.3$Line.of.immunotherapy.clp)))
          
          ## collapse stages 
          data.plot.3$Mets.per.AJCC.01.clp = NA
          data.plot.3$Mets.per.AJCC.01.clp[which(data.plot.3$Mets.per.AJCC.01 %in% 
                                                   c('M1a','M1b') )] = 'M1aM1b'
          data.plot.3$Mets.per.AJCC.01.clp[which(data.plot.3$Mets.per.AJCC.01 %in% 
                                                   c('M1c','M1d') )] = 'M1cM1d'
          
          print(table(data.plot.3$Mets.per.AJCC.01.clp))
          
        }
        
        ## remove clinical cols that have numeric values ...
        data.plot.3 = data.plot.3[,!colnames(data.plot.3) %in%
                                    c('X.eosino.pre','Albumin.pre',
                                      'Line.of.immunotherapy')]
        
        data.plot.3 = reshape2::melt(data.plot.3)
        data.plot.3$tissue = gsub('.weighted..recist.','',data.plot.3$variable)
        data.plot.3$.weighted..recist. = data.plot.3$value
        
        data.plot.3$Best.Response.clp = factor(
          data.plot.3$Best.Response.clp,levels = c('PD','DC','NA')
        )
        
        data.plot.3$Best.Response = factor(
          data.plot.3$Best.Response,levels = c('PD','SD','PR','CR','NA')
        )
        
        write.csv(data.plot.3,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.03.csv'))
        
        ## --------------------------------------------
        
        ## select clinical variables to plot
        data.plot.0.cols.2 = NULL 
        data.plot.0.cols.2 = data.plot.0.cols
        data.plot.0.cols.2[which(data.plot.0.cols.2 == 'X.eosino.pre')] = 
          'X.eosino.pre.high'
        data.plot.0.cols.2[which(data.plot.0.cols.2 == 'Albumin.pre')] = 
          'Albumin.pre.high'
        data.plot.0.cols.2[which(data.plot.0.cols.2 == 'Line.of.immunotherapy')] = 
          'Line.of.immunotherapy.clp'
        data.plot.0.cols.2[which(data.plot.0.cols.2 == 'Mets.per.AJCC.01')] = 
          'Mets.per.AJCC.01.clp'
        
        print(data.plot.0.cols.2)
        
        ## --------------------------------------------
        
        y.min = min(data.plot.3$.weighted..recist., na.rm=T)
        y.max = max(data.plot.3$.weighted..recist., na.rm=T)
        y.min
        y.max
        
        my.out.dir = 'organResRECIST_by_clinicalVar_boxplots'
        if(! dir.exists(my.out.dir)) { dir.create(my.out.dir) }
        
        for(my.clinical.var in data.plot.0.cols.2) {
          
          print(my.clinical.var)
          my.df = NULL 
          
          if(my.clinical.var %in% c('Patient.identifier','Cohort','Cohort2')) { 
            next }
          
          print(table(data.plot.3[,my.clinical.var]))
          
          ## --------------------------------------------
          
          my.df = NULL 
          my.df = data.plot.3[,c(my.clinical.var,
                                 'Patient.identifier','Cohort',
                                 'Best.Response','Best.Response.clp',
                                 'tissue','.weighted..recist.')]
          colnames(my.df)[1] = 'Clinical.var'
          my.df = my.df[!is.na(my.df$Clinical.var),]
          dim(my.df)
          
          write.csv(my.df,
                    file = file.path(my.out.dir,
                                     paste0(clinical.file,'.',my.clinical.var,'.box.csv')))
          
          ## --------------------------------------------
          
          p1 = NULL ## split or no_split
          p1 = ggplot(my.df, aes(Clinical.var, .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            # facet_grid(Best.Response.clp ~ tissue) +
            facet_grid( ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c(
                                 names(table(my.df$Clinical.var))[1],
                                 names(table(my.df$Clinical.var))[2]
                               ))) +
            ggtitle(paste0('ici: ',my.clinical.var)) +
            ylim(y.min, y.max +50)
          p1
          
          pdf(file = file.path(my.out.dir,
                               paste0(clinical.file,'.',my.clinical.var,'.box.pdf')),
              width = 8, height = 6)
          print(p1)
          dev.off()
          
        }
        
        ## --------------------------------------------
        
        my.clinical.var.plot = data.plot.0.cols.2[! data.plot.0.cols.2 %in% c('Patient.identifier','Cohort','Cohort2')]
        
        print(my.clinical.var.plot)
        
        ## making plots to be included in supp figures
        ## set the levels right ...
        if(TRUE) {
          
          p1 = NULL 
          p1 = ggplot(data.plot.3[!is.na(data.plot.3$Age.group),],
                      aes(Age.group,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(Best.Response.clp ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('younger','older'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p1
          
          p2 = NULL 
          p2 = ggplot(data.plot.3[!is.na(data.plot.3$Melanoma.subtype),],
                      aes(Melanoma.subtype,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(Best.Response.clp ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('Cutaneous','Non_cutaneous'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p2
          
          p3 = NULL 
          p3 = ggplot(data.plot.3[!is.na(data.plot.3$Sex),],
                      aes(Sex,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(Best.Response.clp ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('Female','Male'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p3
          
          p4 = NULL 
          p4 = ggplot(data.plot.3[!is.na(data.plot.3$X.eosino.pre.high),],
                      aes(X.eosino.pre.high,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid(Best.Response.clp ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('N','Y'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p4
          
          pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.03.pdf'),
              height = 10, width = 12)
          print(p1 + p2 + p3 + p4 +
                  plot_layout(ncol = 2))
          dev.off()
          
        }
        
        if(TRUE) {
          
          p1 = NULL 
          p1 = ggplot(data.plot.3[!is.na(data.plot.3$Age.group),],
                      aes(Age.group,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid( ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('younger','older'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p1
          
          p2 = NULL 
          p2 = ggplot(data.plot.3[!is.na(data.plot.3$Melanoma.subtype),],
                      aes(Melanoma.subtype,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid( ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('Cutaneous','Non_cutaneous'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p2
          
          p3 = NULL 
          p3 = ggplot(data.plot.3[!is.na(data.plot.3$Sex),],
                      aes(Sex,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid( ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('Female','Male'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p3
          
          p4 = NULL 
          p4 = ggplot(data.plot.3[!is.na(data.plot.3$X.eosino.pre.high),],
                      aes(X.eosino.pre.high,
                          .weighted..recist.)) +
            geom_boxplot(width = 0.6, outlier.shape = NA,
                         color = '#000000', lwd=0.5, fatten=2) +
            geom_jitter(aes(color = Best.Response), width = 0.1, height = 0,
                        size = 1) +
            scale_color_manual(values = plot.colors) +
            theme_pubr(x.text.angle = 90) +
            facet_grid( ~ tissue) +
            # geom_hline(yintercept = c(0, -30,-70), linetype = 'dashed') +
            stat_compare_means(method = 'wilcox.test',
                               comparisons = list(c('N','Y'))) +
            ggtitle('ici') +
            ylim(y.min, y.max +50)
          p4
          
          pdf(file = paste0(clinical.file,'.ptByResOrganRECIST.box.03.2.pdf'),
              height = 8, width = 12)
          print(p1 + p2 + p3 + p4 +
                  plot_layout(ncol = 2))
          dev.off()
          
        }
        
      }
      
    }
    
    ## --------------------------------------------
    
    
  }
  
  ## --------------------------------------------
  
  ## compute stats
  if(TRUE) {
    
    ## fig S3
    if(TRUE) {
      
      ## PD vs DC 
      if(TRUE) {
        
        data.plot.3 = NULL 
        data.plot.3 = read.csv(paste0(clinical.file,'.ptByResOrganRECIST.box.01.csv'),
                               row.names = 1)
        
        ## --------------------------------------------
        
        ## compute stats in all samples 
        data.stats = NULL 
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          my.df.0 = NULL 
          if(my.cohort == 'ici') {
            my.df.0 = data.plot.3[data.plot.3$Cohort %in% c('pd1','ipinivo'),]
          } else {
            my.df.0 = data.plot.3[data.plot.3$Cohort == my.cohort,]
          }
          
          dim(my.df.0)
          
          data.stats.2 = NULL 
          for(my.tissue in names(tissues)) {
            print(my.tissue)
            
            my.df = NULL 
            my.df = my.df.0[my.df.0$tissue == tissues[[my.tissue]] &
                              !is.na(my.df.0$.weighted..recist.),]
            my.df = my.df[!is.na(my.df$Best.Response.clp),]
            dim(my.df)
            print(table(my.df[,'Best.Response.clp']))
            
            ## make sure there are 2 groups!
            if(length(table(my.df[,'Best.Response.clp'])) < 2) {
              print('not two groups left for stats!')
              next
            }
            
            my.stats = NULL 
            my.stats = tidy(wilcox.test(
              my.df$.weighted..recist.[my.df$Best.Response.clp =='PD'],
              my.df$.weighted..recist.[my.df$Best.Response.clp =='DC'])
            )
            
            data.stats.2 =
              rbind(data.stats.2,
                    data.frame(Cohort = my.cohort,
                               Clinical.var = 'Best.Response.clp',
                               Response = 'allsm_nosplit',
                               Tissue = my.tissue,
                               GroupA = 'PD',
                               GroupB = 'DC',
                               Sample.total = nrow(my.df),
                               GroupA.total = table(my.df$Best.Response.clp)[1],
                               GroupB.total = table(my.df$Best.Response.clp)[2],
                               my.stats,
                               stringsAsFactors = F)
              )
          }
          
          
          
          data.stats.2$p.adj = p.adjust(data.stats.2$p.value,
                                        method = 'fdr')
          data.stats = rbind(data.stats,
                             data.stats.2)
          
        }
        
        write.csv(data.stats,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.01.wilcox.csv'))
        
      }
      
      ## --------------------------------------------
      
      ## organ heterogeneity 
      if(TRUE) {
        
        data.stats = NULL 
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          my.df = NULL 
          my.df = read.csv(paste0(clinical.file,'.ptByResOrganRECISTsd.',
                                  my.cohort,'.box.01.2_4.csv'), row.names = 1)
          
          table(my.df$Response.organ.mixed)
          
          ## --------------------------------------------
          
          my.comps = NULL 
          my.comps = list(c('mixedResponse','uniformProgression'),
                          c('mixedResponse',
                            'uniformDiseaseControl'),
                          c('uniformProgression',
                            'uniformDiseaseControl'))
          
          if(my.cohort == 'braf') {
            my.comps = list(c('mixedResponse',
                              'uniformDiseaseControl'))
          }
          
          ## --------------------------------------------
          
          my.stats = NULL 
          for(my.comp in my.comps) {
            print(my.comp)
            
            my.stats = rbind(
              my.stats,
              data.frame(Cohort = my.cohort,
                         Comp = paste(my.comp, collapse = '_versus_'),
                         GroupA.total = sum(my.df$Response.organ.mixed %in% my.comp[1]),
                         GroupB.total = sum(my.df$Response.organ.mixed %in% my.comp[2]),
                         tidy(wilcox.test(my.df$.weighted..recist.[
                           my.df$Response.organ.mixed %in% my.comp[1]],
                           my.df$.weighted..recist.[my.df$Response.organ.mixed %in% my.comp[2]])
                         ),
                         stringsAsFactors = F
              )
            )
            
          }
          
          my.stats$p.adj = p.adjust(my.stats$p.value,method = 'fdr')
          
          data.stats = rbind(data.stats,
                             my.stats)
          
          
        }
        
        data.stats
        
        write.csv(data.stats,
                  file = paste0(clinical.file,'..ptByResOrganRECISTsd.box.01.2_4.wilcox.csv'))
      }
      
    }
    
    ## --------------------------------------------
    
    ## fig S4
    if(TRUE) {
      
      ## lesion 01 vs 02 from the same organ 
      if(TRUE) {
        
        data.plot.3 = NULL 
        data.plot.3 = read.csv(paste0(clinical.file,'.ptByResOrganRECIST.box.02.paired.csv'),
                               row.names = 1)
        
        ## --------------------------------------------
        
        data.stats = NULL 
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          my.df.0 = NULL 
          if(my.cohort == 'ici') {
            my.df.0 = data.plot.3[data.plot.3$Cohort %in% c('pd1','ipinivo'),]
          } else {
            my.df.0 = data.plot.3[data.plot.3$Cohort == my.cohort,]
          }
          
          dim(my.df.0)
          
          data.stats.2 = NULL 
          for(my.tissue in names(tissues)) {
            print(my.tissue)
            
            my.df = NULL 
            my.df = my.df.0[my.df.0$tissue == tissues[[my.tissue]],]
            dim(my.df)
            print(table(my.df$lesion))
            
            ## --------------------------------------------
            
            ## if only one pair; or pairs do not match
            if((table(my.df$lesion)[1] <= 1 | 
                table(my.df$lesion)[2] <= 1) | (
                  table(my.df$lesion)[1] != table(my.df$lesion)[2] 
                )) {
              print('if only one pair or pairs do not match, skip')
              next 
              
            }
            
            ## --------------------------------------------
            
            my.df = reshape2::dcast(Cohort.Deidentified.name.tissue ~ lesion,
                                    value.var = 'recist',
                                    data = my.df, 
                                    fun.aggregate = sum)
            my.stats = NULL 
            my.stats = tidy(wilcox.test(my.df$lesion01, my.df$lesion02, 
                                        paired = T))
            
            data.stats.2 =
              rbind(data.stats.2,
                    data.frame(Cohort = my.cohort,
                               Response = 'allsm_nosplit',
                               Tissue = my.tissue,
                               Sample.total = nrow(my.df) * 2,
                               Lesion0102.pairs = nrow(my.df),
                               my.stats,
                               stringsAsFactors = F)
              )
          }
          
          
          
          data.stats.2$p.adj = p.adjust(data.stats.2$p.value,
                                        method = 'fdr')
          data.stats = rbind(data.stats,
                             data.stats.2)
          
        }
        
        write.csv(data.stats,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.02.paired.wilcox.csv'))
        
      }
      
      ## --------------------------------------------
      
      ## absolute diff between lesion 01 and lesion 02 
      if(TRUE) {
        
        data.stats = NULL 
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          my.file.0 = NULL 
          my.file.0 = paste0(clinical.file, '.ptByResOrganRECIST.',
                             my.cohort,'.box.02.paired.absdiff.csv')
          print(my.file.0)
          
          my.df.0 = NULL 
          my.df.0 = read.csv(my.file.0, row.names = 1)
          
          dim(my.df.0)
          
          data.stats.2 = NULL 
          for(my.tissue in names(tissues)) {
            print(my.tissue)
            
            my.df = NULL 
            my.df = my.df.0[my.df.0$tissue == tissues[[my.tissue]],,drop=F]
            dim(my.df)
            print(table(my.df$Best.Response.clp))
            
            ## --------------------------------------------
            
            ## make sure there are 2 groups!
            if(length(table(my.df[,'Best.Response.clp'])) < 2) {
              print('not two groups left for stats!')
              next
            }
            
            ## --------------------------------------------
            
            ## log transform, same as the figure 
            my.df$lesion0102.abs.diff.log = log10(my.df$lesion0102.abs.diff + 1)
            
            my.stats = NULL 
            my.stats = tidy(wilcox.test(
              my.df$lesion0102.abs.diff.log[my.df$Best.Response.clp =='PD'],
              my.df$lesion0102.abs.diff.log[my.df$Best.Response.clp =='DC'])
            )
            
            data.stats.2 =
              rbind(data.stats.2,
                    data.frame(Cohort = my.cohort,
                               Clinical.var = 'Best.Response.clp',
                               Response = 'allsm_nosplit',
                               Tissue = my.tissue,
                               GroupA = 'PD',
                               GroupB = 'DC',
                               Sample.total = nrow(my.df),
                               GroupA.total = table(my.df$Best.Response.clp)[1],
                               GroupB.total = table(my.df$Best.Response.clp)[2],
                               my.stats,
                               stringsAsFactors = F)
              )
            
            
          }
          
          
          
          data.stats.2$p.adj = p.adjust(data.stats.2$p.value,
                                        method = 'fdr')
          data.stats = rbind(data.stats,
                             data.stats.2)
          
        }
        
        write.csv(data.stats,
                  file = paste0(clinical.file,
                                '.ptByResOrganRECIST.box.02.paired.absdiff.wilcox.csv'))
      }
      
    }
    
    ## --------------------------------------------
    
    ## fig S5
    if(TRUE) {
      
      ## PD vs DC boxplot stats: select clinical variables 
      if(TRUE) {
        data.plot.3 = NULL 
        data.plot.3 = read.csv(file = paste0(clinical.file,'.ptByResOrganRECIST.box.03.csv'),
                               row.names = 1)
        
        ## --------------------------------------------
        
        my.clinical.var.plot = c(
          'Age.group','Melanoma.subtype','Sex','X.eosino.pre.high'
        )
        
        ## --------------------------------------------
        
        ## compute stats in overall response PD and DC samples subsets
        data.stats = NULL 
        for(my.var in my.clinical.var.plot) {
          print(my.var)
          
          data.stats.2 = NULL 
          for(my.response in c('PD','DC')) {
            print(my.response)
            
            for(my.tissue in names(tissues)) {
              print(my.tissue)
              
              my.df = NULL 
              my.df = data.plot.3[data.plot.3$Best.Response.clp == my.response &
                                    data.plot.3$tissue == tissues[[my.tissue]] &
                                    !is.na(data.plot.3$.weighted..recist.),]
              dim(my.df)
              print(table(my.df[,my.var]))
              
              ## make sure there are 2 groups!
              if(length(table(my.df[,my.var])) < 2) { 
                print('not two groups left for stats!')
                next 
              }
              
              my.stats = NULL 
              my.stats = tidy(wilcox.test(
                my.df$.weighted..recist.[my.df[,my.var] ==
                                           names(table(my.df[,my.var]))[1]],
                my.df$.weighted..recist.[my.df[,my.var] ==
                                           names(table(my.df[,my.var]))[2]])
              )
              
              data.stats.2 =
                rbind(data.stats.2,
                      data.frame(Cohort = 'ici',
                                 Clinical.var = my.var,
                                 Response = my.response,
                                 Tissue = my.tissue,
                                 GroupA = names(table(my.df[,my.var]))[1],
                                 GroupB = names(table(my.df[,my.var]))[2],
                                 Sample.total = nrow(my.df),
                                 GroupA.total = table(my.df[,my.var])[1],
                                 GroupB.total = table(my.df[,my.var])[2],
                                 my.stats,
                                 stringsAsFactors = F)
                )
            }
            
          }
          
          data.stats.2$p.adj = p.adjust(data.stats.2$p.value,
                                        method = 'fdr')
          data.stats = rbind(data.stats,
                             data.stats.2)
          
        }
        
        write.csv(data.stats,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.03.wilcox.csv'))
        
        ## --------------------------------------------
        
        ## compute stats in all samples 
        data.stats = NULL 
        for(my.var in my.clinical.var.plot) {
          print(my.var)
          
          data.stats.2 = NULL 
          for(my.tissue in names(tissues)) {
            print(my.tissue)
            
            my.df = NULL 
            my.df = data.plot.3[data.plot.3$tissue == tissues[[my.tissue]] &
                                  !is.na(data.plot.3$.weighted..recist.),]
            dim(my.df)
            print(table(my.df[,my.var]))
            
            ## make sure there are 2 groups!
            if(length(table(my.df[,my.var])) < 2) { 
              print('not two groups left for stats!')
              next 
            }
            
            my.stats = NULL 
            my.stats = tidy(wilcox.test(
              my.df$.weighted..recist.[my.df[,my.var] ==
                                         names(table(my.df[,my.var]))[1]],
              my.df$.weighted..recist.[my.df[,my.var] ==
                                         names(table(my.df[,my.var]))[2]])
            )
            
            data.stats.2 =
              rbind(data.stats.2,
                    data.frame(Cohort = 'ici',
                               Clinical.var = my.var,
                               Response = 'allsm_nosplit',
                               Tissue = my.tissue,
                               GroupA = names(table(my.df[,my.var]))[1],
                               GroupB = names(table(my.df[,my.var]))[2],
                               Sample.total = nrow(my.df),
                               GroupA.total = table(my.df[,my.var])[1],
                               GroupB.total = table(my.df[,my.var])[2],
                               my.stats,
                               stringsAsFactors = F)
              )
          }
          
          
          
          data.stats.2$p.adj = p.adjust(data.stats.2$p.value,
                                        method = 'fdr')
          data.stats = rbind(data.stats,
                             data.stats.2)
          
        }
        
        write.csv(data.stats,
                  file = paste0(clinical.file,'.ptByResOrganRECIST.box.03.2.wilcox.csv'))
      }
      
      ## --------------------------------------------
      
      ## PD vs DC boxplot stats: all clinical variables 
      ## also make a heatmap instead of boxplots
      if(TRUE) {
        
        data.plot.3 = NULL 
        data.plot.3 = read.csv(paste0(clinical.file,'.ptByResOrganRECIST.box.03.csv'),
                               row.names = 1)
        
        my.out.dir = 'organResRECIST_by_clinicalVar_boxplots'
        
        ## --------------------------------------------
        
        ## compute stats in all samples 
        data.stats = NULL 
        for(my.clinical.var in data.plot.0.cols.2) {
          print(my.clinical.var)
          
          if(my.clinical.var %in% c('Patient.identifier','Cohort','Cohort2')) { 
            next }
          
          print(table(data.plot.3[,my.clinical.var]))
          
          ## --------------------------------------------
          
          my.file = NULL
          my.file = file.path(my.out.dir,
                              paste0(clinical.file,'.',my.clinical.var,'.box.csv'))
          print(my.file)
          
          my.df.0 = NULL 
          my.df.0 = read.csv(my.file, row.names = 1)
          dim(my.df.0)
          
          ## --------------------------------------------
          
          data.stats.2 = NULL 
          for(my.tissue in names(tissues)) {
            print(my.tissue)
            
            my.df = NULL 
            my.df = my.df.0[my.df.0$tissue == tissues[[my.tissue]] &
                              !is.na(my.df.0$.weighted..recist.),]
            dim(my.df)
            print(table(my.df[,'Clinical.var']))
            
            ## make sure there are 2 groups!
            if(length(table(my.df[,'Clinical.var'])) < 2) { 
              print('not two groups left for stats!')
              next 
            }
            
            my.stats = NULL 
            my.stats = tidy(wilcox.test(
              my.df$.weighted..recist.[my.df[,'Clinical.var'] ==
                                         names(table(my.df[,'Clinical.var']))[1]],
              my.df$.weighted..recist.[my.df[,'Clinical.var'] ==
                                         names(table(my.df[,'Clinical.var']))[2]])
            )
            
            
            ## for plotting log2fc heatmap downstream....
            my.mean1 = NULL 
            my.mean2 = NULL 
            my.log2fc = NULL 
            
            ## convert to rank...
            my.df$.weighted..recist.rank = rank(my.df$.weighted..recist.,
                                                na.last = NA)
            my.df$.weighted..recist.rank.log = log2(my.df$.weighted..recist.rank)
            
            my.mean1 = mean(my.df$.weighted..recist.rank.log[my.df[,'Clinical.var'] ==
                                                               names(table(my.df[,'Clinical.var']))[1]],
                            na.rm = T)
            my.mean2 = mean(my.df$.weighted..recist.rank.log[my.df[,'Clinical.var'] ==
                                                               names(table(my.df[,'Clinical.var']))[2]],
                            na.rm = T)
            my.log2fc = my.mean1 - my.mean2
            
            data.stats.2 =
              rbind(data.stats.2,
                    data.frame(Cohort = 'ici',
                               Clinical.var = my.clinical.var,
                               Response = 'allsm_nosplit',
                               Tissue = my.tissue,
                               GroupA = names(table(my.df[,'Clinical.var']))[1],
                               GroupB = names(table(my.df[,'Clinical.var']))[2],
                               Sample.total = nrow(my.df),
                               GroupA.total = table(my.df[,'Clinical.var'])[1],
                               GroupB.total = table(my.df[,'Clinical.var'])[2],
                               my.stats,
                               MeanA = my.mean1,
                               Mean2 = my.mean2,
                               Log2FC = my.log2fc,
                               stringsAsFactors = F)
              )
          }
          
          data.stats.2$p.adj = p.adjust(data.stats.2$p.value,
                                        method = 'fdr')
          data.stats = rbind(data.stats,
                             data.stats.2)
          
        }
        
        write.csv(data.stats,
                  file = paste0(clinical.file,
                                '.ptByResOrganRECIST.box.03.2.allclinvars.wilcox.csv'))
        
        print(data.stats[data.stats$p.adj < 0.10,,drop=F])
        
      }
      
    }
    
    
  }
  
}