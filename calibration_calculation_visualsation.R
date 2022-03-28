
# Description of the function

calculate_calibration <- function( df_data, IS, cal.ref.pnt, delta_linearity, alpha_IR) {
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  # Read the batch list 
  batch.list=unique(df_data$batch)
  compound.batch.list=list()
  
  for (b in seq_along(batch.list)){
    compound.batch.list[[b]]=unique(df_data$compound[df_data$batch==batch.list[b]])
  }
  
  
# 0. Check the internal standard in all sample ####
  #(OPTIONAL because: )
  # Ask if there is an internal Standard? for each batch? 
  
  
  # calculate deviation, find outliers?
  
  # Ask if want to remove the internal standard
  df_data=subset(df_data, compound %!in% IS)
  
  for (b in seq_along(batch.list)){
    compound.batch.list[[b]] <- compound.batch.list[[b]][compound.batch.list[[b]]!=IS]
  }
  
  
# 1. Ion ratio ####
  # Calculate Ion Ration and replace NA by 0
  df_data <- df_data %>%
    mutate(IR = area2/area1,
           IR=if_else(is.na(IR)|is.infinite(IR), 0,IR) )   
  
  # //!\\ CHECK Area1> area2 ??
  
# 2. Calibration calculation ####
  # Summary of calibration results
  df_cal=subset(df_data, an_type == "cal")
  
  
  for(b  in seq_along(batch.list)){
    for(c in seq_along(compound.batch.list[[b]])){
      
      # Create Data Frame of calibration for batch b, compound c: 
      temp.cal.bc= subset( df_data, batch==batch.list[b] &
                             compound==compound.batch.list[[b]][c] &
                             an_type=="cal")
      
      # * Calibration linearity ----
      temp.cal.bc= temp.cal.bc %>%
        mutate( area.ref.mean = mean( area1[temp.cal.bc$cal.level==cal.ref.pnt]),
                linearity=area1 *cal.ref.pnt /cal.level /area.ref.mean,
                linearity.ck=if_else(linearity> (1-delta_linearity) & linearity< (1+delta_linearity), "1","0"),
                
                # * Ion Ratio reference  ----
                IR_ref0=mean( IR[temp.cal.bc$cal.level==cal.ref.pnt]),  # First take IR of the ref cal point for ref
                IR_ck0=if_else(IR> 0.5*IR_ref0 & IR< 1.5*IR_ref0, "1","0"), # cut out all cal not in .5-1.5 IR_ref0
                
                IR_ref=mean(subset( temp.cal.bc, linearity.ck=="1" & IR_ck0=="1" )$IR), # Calculate Ion Ration reference with the calibration in the linearity range and .55-1.45 IR_ref0
                IR_ck=if_else(IR> (1-alpha_IR)*IR_ref & IR< (1+alpha_IR)*IR_ref, "1","0")
        ) 
      
      
      subset(temp.cal.bc, linearity.ck=="1" )
      
      
      # linear model 
      l.model=lm(cal.level~area1, subset(temp.cal.bc,  linearity.ck=="1" & IR_ck=="1"  )   )     
      
      temp.cal.bc=temp.cal.bc %>%
        mutate(intercept=l.model[[1]][1],
               slope=l.model[[1]][2],
               R2=summary(l.model)$r.squared,
               detec.conc=area1*slope+intercept )
      
      df_cal[df_cal$batch==batch.list[b] & 
               df_cal$compound==compound.batch.list[[b]][c],
             c("area.ref.mean","linearity", "linearity.ck","IR_ck","intercept","slope","R2","detec.conc")] =  temp.cal.bc[, c("area.ref.mean","linearity","linearity.ck","IR_ck","intercept","slope","R2","detec.conc")]
      
    }
  }
  
  return(df_cal)
} # end of the function calibration 

# Optional: plot the calibation, all, per batch or specific compond


calibration_plot=function( df_cal) {#improve with ether select batch and or compound and or all 
  
  # Read the batch and compound list 
  batch.list=unique(df_cal$batch)
  compound.batch.list=list()
      
  for (b in seq_along(batch.list)){
    compound.batch.list[[b]]=unique(df_cal$compound[df_cal$batch==batch.list[b]])
  }
      
      
      
  for(b  in seq_along(batch.list)){
    for(c in seq_along(compound.batch.list[[b]])){
      
      
      
      
      #R2 to be displayed on graph
      
      temp.df.cal.bc= subset(  df_cal, batch==batch.list[b] &
                                 compound==compound.batch.list[[b]][c] )    
      
      R2 <- as.character(as.expression( substitute(
        italic(R)^2~"="~r2,
        list(r2 = format( unique(temp.df.cal.bc$R2) , digits = 3) )
      )))
      
      slope =  unique(temp.df.cal.bc$slope)
      
      intercept = unique(temp.df.cal.bc$intercept)
      
      PLOT=ggplot(  temp.df.cal.bc, aes(x=area1, y=cal.level, color=linearity.ck))+
        geom_point()+
        geom_abline( slope = slope, intercept =  intercept, size=1  )+
        geom_text(x = 1000, y = 10, hjust = 0 , parse = TRUE, colour="black", show.legend = FALSE, size=5,
                  label =  R2  )+
        ggtitle(label = paste(compound.batch.list[[b]][c],batch.list[b], sep = "_"  ))+
        theme_minimal()+
        theme(
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.x  = element_text(size=9),
          axis.text.y =  element_text(size=9))
      print(PLOT)
      
    }
  }
} # end calibration_plot 


# IS=c("13C_Caffeine","13C_Caffeine","13C_Caffeine")
# cal.ref.pnt=c(1, 1, 1)
# delta_linearity=0.3
# alpha_IR=0.3
# 
# df_cal=calculate_calibration(df_data, IS, cal.ref.pnt, delta_linearity, alpha_IR)
# calibration_plot(df_cal)