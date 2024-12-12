

## Figure 1
if(TRUE) {
  
  ## prepare input
  if(TRUE) {
    
    data.plot = NULL 
    data.plot = clinical
    data.plot$Cohort2 = data.plot$Cohort
    data.plot$Cohort2[which(data.plot$Cohort %in% c('pd1','ipinivo'))] = 'ici'
    
    ## --------------------------------------------
    
    table(data.plot$Cohort)
    
    table(data.plot$Best.Response.by.RECIST)
    
    sum(is.na(data.plot$Best.Response.by.RECIST)) ## 0
    
    
    
  }
  
  ## --------------------------------------------
  
  ## draw figures
  if(TRUE) {
    
    ## disease control rate
    if(TRUE) {
      
      ## bargraph
      if(TRUE) {
        
        data.plot.2 = NULL 
        x = NULL
        x = data.frame(prop.table(table(data.plot[data.plot$Cohort2 =='ici',
                                                  c('Cohort2','Best.Response.by.RECIST')])))
        colnames(x)[1] = 'Cohort'
        data.plot.2 = rbind(
          x,
          data.frame(prop.table(table(data.plot[data.plot$Cohort=='pd1',
                                                c('Cohort','Best.Response.by.RECIST')]))),
          data.frame(prop.table(table(data.plot[data.plot$Cohort=='ipinivo',
                                                c('Cohort','Best.Response.by.RECIST')]))),
          data.frame(prop.table(table(data.plot[data.plot$Cohort=='braf',
                                                c('Cohort','Best.Response.by.RECIST')])))
        )
        
        data.plot.2$Pct = round(data.plot.2$Freq * 100, 
                                digits = 0)
        
        data.plot.2$Best.Response.by.RECIST = factor(data.plot.2$Best.Response.by.RECIST,
                                                     levels = c('PD','SD','PR','CR',
                                                                'NA'))
        
        write.csv(data.plot.2,
                  file = paste0(clinical.file,'.ptByResponse.bar.csv'))
        
        
        p1 = NULL 
        p1 = ggplot(data.plot.2, aes(Cohort,
                                     Pct)) +
          geom_bar(aes(fill = Best.Response.by.RECIST), stat = 'identity',
                   position = 'stack',
                   color = '#000000', width = 0.6) +
          scale_fill_manual(values = plot.colors) + 
          theme_pubr() 
        
        p1
        
        pdf(file = paste0(clinical.file,'.ptByResponse.bar.pdf'),
            height = 4, width = 6)
        print(p1)
        dev.off()
        
      }
      
      ## --------------------------------------------
      
      ## donut chart
      if(TRUE) {
        
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          ## --------------------------------------------
          
          data.plot.2 = NULL
          
          if(my.cohort == 'ici') {
            data.plot.2 = data.frame(prop.table(table(data.plot[
              data.plot$Cohort2 ==my.cohort, c('Cohort2','Best.Response.by.RECIST')])))
            colnames(data.plot.2)[1] = 'Cohort'
          } else {
            data.plot.2 = data.frame(prop.table(table(data.plot[
              data.plot$Cohort ==my.cohort, c('Cohort2','Best.Response.by.RECIST')])))
          }
          
          data.plot.2$DCR = NA
          data.plot.2$DCR[which(data.plot.2$Best.Response.by.RECIST %in% c('PD'))] = 'PD'
          data.plot.2$DCR[which(data.plot.2$Best.Response.by.RECIST %in% c('SD','PR','CR'))] = 'DC'
          data.plot.2 = data.plot.2 %>% group_by(Best.Response.by.RECIST, DCR) %>% 
            summarise(n = sum(Freq))
          print(data.plot.2)
          
          write.csv(data.plot.2,
                    file = paste0(clinical.file,'.',my.cohort,'.ptByResponse.donut.csv'))
          
          pdf(file = paste0(clinical.file,'.',my.cohort,'.ptByResponse.donut.pdf'),
              height = 4, width = 6)
          PieDonut(data.plot.2, aes(DCR, Best.Response.by.RECIST, count=n), 
                   title = my.cohort,
                   ratioByGroup = FALSE, r0 = 0.45, r1 = 0.9)
          dev.off()
          
        }
        
        
        
      }
      
    }
    
    ## --------------------------------------------
    
    ## disease control varied by line of therapy 
    if(TRUE) {
      
      ## bargraph
      if(TRUE) {
        
        data.plot.2 = NULL 
        
        for(my.cohort in c('ici','pd1','ipinivo','braf')) {
          print(my.cohort)
          
          ## --------------------------------------------
          
          my.df = NULL
          
          if(my.cohort == 'ici') {
            my.df = data.plot[
              data.plot$Cohort2 ==my.cohort, c('Cohort2','Best.Response.by.RECIST',
                                               'Line.of.immunotherapy')]
            colnames(my.df)[1] = 'Cohort'
          } else {
            my.df =data.plot[
              data.plot$Cohort==my.cohort,c('Cohort','Best.Response.by.RECIST',
                                            'Line.of.immunotherapy')]
          }
          
          print(table(my.df$Line.of.immunotherapy))
          
          if(my.cohort %in% c('ici','pd1','ipinivo')) {
            
            my.df$Line.of.immunotherapy.clp = NA
            my.df$Line.of.immunotherapy.clp[
              which(my.df$Line.of.immunotherapy %in% c(1))
            ] = 'FirstLine'
            my.df$Line.of.immunotherapy.clp[
              which(my.df$Line.of.immunotherapy %in% c(2,3,4))
            ] = 'TwoOrMoreLines'
            
            print(sum(is.na(my.df$Line.of.immunotherapy.clp)))
            
            ## convert to %
            my.df = rbind(
              data.frame(prop.table(table(
                my.df[my.df$Line.of.immunotherapy.clp == 'FirstLine',
                      c('Cohort','Best.Response.by.RECIST','Line.of.immunotherapy.clp')]))),
              data.frame(prop.table(table(
                my.df[my.df$Line.of.immunotherapy.clp == 'TwoOrMoreLines',
                      c('Cohort','Best.Response.by.RECIST','Line.of.immunotherapy.clp')])))
              
            )
            
          } else if(my.cohort %in% c('braf')) {
            my.df$Line.of.immunotherapy.clp = NA
            my.df$Line.of.immunotherapy.clp[
              which(my.df$Line.of.immunotherapy %in% c(0))
            ] = 'NoPriorImtx'
            my.df$Line.of.immunotherapy.clp[
              which(my.df$Line.of.immunotherapy %in% c(1,2,3,4))
            ] = 'YesPriorImtx'
            
            print(sum(is.na(my.df$Line.of.immunotherapy.clp)))
            
            ## convert to %
            my.df = rbind(
              data.frame(prop.table(table(
                my.df[my.df$Line.of.immunotherapy.clp == 'NoPriorImtx',
                      c('Cohort','Best.Response.by.RECIST','Line.of.immunotherapy.clp')]))),
              data.frame(prop.table(table(
                my.df[my.df$Line.of.immunotherapy.clp == 'YesPriorImtx',
                      c('Cohort','Best.Response.by.RECIST','Line.of.immunotherapy.clp')])))
              
            )
          }
          
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
                  file = paste0(clinical.file,'.ptByResponseByLineImtx.bar.csv'))
        
        
        p1 = NULL 
        p1 = ggplot(data.plot.2, aes(Line.of.immunotherapy.clp,
                                     Pct)) +
          geom_bar(aes(fill = Best.Response.by.RECIST), stat = 'identity',
                   position = 'stack',
                   color = '#000000', width = 0.6) +
          scale_fill_manual(values = plot.colors) + 
          facet_grid(~ Cohort, scales = 'free_x') +
          theme_pubr(x.text.angle = 90) 
        
        p1
        
        pdf(file = paste0(clinical.file,'.ptByResponseByLineImtx.bar.pdf'),
            height = 4, width = 10)
        print(p1)
        dev.off()
        
      }
      
    }
    
  }
  
}
