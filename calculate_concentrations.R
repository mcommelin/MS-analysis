#Analyse GC-MS, LC-MS and specific compound results

#LC calculation

#Process compounds - user specific?
rm(list=ls()) # cleaning console
graphics.off() # cleaning plots

'%!in%' <- function(x,y)!('%in%'(x,y))

meta <-  TRUE # true is meta in sample txt 
source("load_data.R")


IS=c("13C_Caffeine","13C_Caffeine","13C_Caffeine")
cal.ref.pnt=c(1, 1, 1)

df_data$an_type[grep("Std", df_data$sample_text)]="Std"

delta_linearity=0.3
alpha_IR=0.3

calculate.concentration <- function( df_data, IS, cal.ref.pnt, delta_linearity, alpha_IR) {  }



batch.list=unique(df_data$batch)
compound.batch.list=list()

for (b in seq_along(batch.list)){
compound.batch.list[[b]]=unique(df_data$compound[df_data$batch==batch.list[b]])
}

df_data$corr.factor=4.4 
df_data$dilution.factor=1

# The calculation involves 5 different steps : 
#   1. Find linearity range with the Calibration /Batch/compound
#     -compare the calibration per batch?

#   2. Calculate the IR range based on Area1/ Area2
#     -Check that all calibration in the linearity range is in the IR range
#     -Check that all samples and QC in the linearity range is in the IR range

#   3. Calculate the Limit of detection : Lowest calibration point, 
#                                         in IR range,
#                                         in linearity range
#     Retun table of lOQ /Batch/ Compound
#     Return vector of highest LOQ/ Compound over all batches

#   4. Calculate the concentrations in samples, QCs and Spikes
#     -Concentration in vial    C.vial
#     -Concentration in matrix  C.matrix 
#     -Add column for units???

#   5. Calculate the Recovery based on the concentration in QC
#     -table of QCs and recovery 
#     -stability of recovery /Batch/ Compound
#     -stability of recovery / Compound over all batches
#     -calculate a column of corrected concentration C.matrix.cor


# 0. Check the internal standard in all sample ####
#(OPTIONAL because: )
  # Ask if there is an internal Standard? for each batch? 

    
  # calculate deviation, find autliers?
    
  # Ask if want to remove the internal standard
    df_data=subset(df_data, compound %!in% IS)
    
    for (b in seq_along(batch.list)){
      compound.batch.list[[b]]=compound.batch.list[[b]][compound.batch.list[[b]]!=IS[b]]
    }

# 1. Ion ratio ####
  # Calculate Ion Ration and replace NA by 0
  df_data <- df_data %>%
    mutate(IR = area2/area1,
           IR=if_else(is.na(IR)|is.infinite(IR), 0,IR) )   
  
  # //!\\ CHECK Area1> area2 ??
    
# 2. Calibration linearity range ####
  
    
  # * Find calibration levels ---- 
    df_data<-df_data%>%
    mutate(cal.level = if_else(an_type == "cal" | an_type == "Cal",
                          as.numeric(str_extract(sample_text, "(\\d+\\.\\d*)|(\\d+)")), 0),
         
           cal.level = if_else(is.na(cal.level), 0, cal.level) )
  
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
        mutate( area.ref.mean= mean( area1[temp.cal.bc$cal.level==cal.ref.pnt[b]]),
                linearity=area1 *cal.ref.pnt[b] /cal.level /area.ref.mean,
                linearity.ck=if_else(linearity> (1-delta_linearity) & linearity< (1+delta_linearity), "1","0"),
                
                # * Ion Ration reference  ----
                IR_ref=mean(subset( temp.cal.bc, linearity.ck==1)$IR), # Calculate Ion Ration reference with the calibration in the linearity range
                IR_ck=if_else(IR> (1-alpha_IR)*IR_ref & IR< (1+alpha_IR)*IR_ref, "1","0"),
        ) 
     
     
     subset(temp.cal.bc, linearity.ck=="1" )
      
      
      # linear model 
      l.model=lm(cal.level~area1, subset(temp.cal.bc,  linearity.ck==1  )   )     
      
      temp.cal.bc=temp.cal.bc %>%
        mutate(intercept=l.model[[1]][1],
               slope=l.model[[1]][2],
               R2=summary(l.model)$r.squared,
               detec.conc=area1*slope+intercept )
        
      df_cal[df_cal$batch==batch.list[b] & 
             df_cal$compound==compound.batch.list[[b]][c],
             c("area.ref.mean","linearity", "linearity.ck","intercept","slope","R2","detec.conc")] =  temp.cal.bc[, c("area.ref.mean","linearity","linearity.ck","intercept","slope","R2","detec.conc")]
        

       
    }
    }

} # end of the function calibration 

  # Optional: plot the calibation, all, per batch or specific compond
    
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
    

# 3. Limit of quantification ####
  
  
  
# 4. Concentrations in samples and QCs  ####
  
  
# 5. QCs recovery ####




