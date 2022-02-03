#Analyse GC-MS, LC-MS and specific compound results

#LC calculation

#Process compounds - user specific?
rm(list=ls()) # cleaning console
graphics.off() # cleaning plots

'%!in%' <- function(x,y)!('%in%'(x,y))

meta <-  TRUE # true is meta in sample txt 

source("load_data.R")
df_data=load_raw_data("data_LC", "\t", meta)

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

# 0. Check the internal standard in all sample (OPTIONAL because: 
  # Ask if there is an internal Standard? for each batch? 
    IS=c("13C_Caffeine","13C_caffeine","13C_caffeine")
    
  # calculate deviation, find autliers?
    
  # Ask if want to remove the internal standard
    df_data=subset(df_data, compound %!in% IS)
    
    for (b in seq_along(df_data)){
      compound.batch.list[[b]]=compound.batch.list[[b]][compound.batch.list[[b]]!=IS]
    }
    compound.batch.list==IS
    subset( compound.batch.list[[b]], compound!= IS)
    
    
# 1. Calibration linearity range ####
  
# Find calibration levels: 
  df_data=df_data%>%
  mutate(cal.level = if_else(an_type == "cal" | an_type == "Cal",
                        as.numeric(str_extract(sample_text, "(\\d+\\.\\d*)|(\\d+)")), 0),
       
         cal.level = if_else(is.na(cal.level), 0, cal.level) )


  for(b  in 1:(length(Batchnames)-1)){
    for(c in seq_along(compound.batch[[i]])){
      # Create Data Frame calibration : 
      sub.cal.bc= subset( df_data, batch==batch.list[b] &
                          compound==compound.batch.list[[b]][c] &
                          an_type=="cal")
      
    }
    }


# 2. Ion ratio ####
  # Calculate Ion Ration and replace NA by 0
  df_data <- df_data %>%
    mutate(IR = Area1/Area2,
           IR=if_else(is.na(IR)|is.infinite(IR), 0,IR) )
  
  df_data$IR

# 3. Limit of quantification ####
  
  
# 4. Concentrations in samples and QCs  ####
  
  
# 3. QCs recovery ####




