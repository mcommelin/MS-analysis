#   3. Calculate the Limit of detection : Lowest calibration point, 
#                                         in IR range,
#                                         in linearity range
#     Retun table of lOQ /Batch/ Compound
#     Return vector of highest LOQ/ Compound over all batches


loq_calculation=function( df_cal) {
  
  # Read the batch and compound list 
  batch.list=unique(df_cal$batch)
  compound.batch.list=list()
  
  for (b in seq_along(batch.list)){
    compound.batch.list[[b]]=unique(df_cal$compound[df_cal$batch==batch.list[b]])
  }
  
  compound.list.all=unique( unlist( compound.batch.list))
    
  # LOQ in vial ----
    # Create a Data Frame for lod in vial [ng/mL] for each batch and each compounds:
    loq_vial<-data.frame(matrix(0,length(batch.list),length( compound.list.all)+1))
    names(loq_vial)<-c("batch",compound.list.all)
    loq_vial$batch=batch.list 
    
    
  for(b  in seq_along(batch.list)){
    for(c in seq_along(compound.list.all)){
      
      
      loq_vial[loq_vial$batch==batch.list[b], compound.list.all[c]]= min( df_cal$cal.level[  # minimum calibration point
        df_cal$batch==batch.list[b] & df_cal$compound==compound.list.all[c] &    # for batch b and compound c 
        df_cal$area1>100 &                                             # additional criteria based on expert judgment
        df_cal$linearity.ck=="1"  & df_cal$IR_ck=="1"] )                  # inside the linearity and IR range
  
    }
  }
    
    loq_vial <- do.call(data.frame,                      # Replace Inf in data by NA
                       lapply(loq_vial,
                              function(x) replace(x, is.infinite(x), NA)))
    
   return(loq_vial)
    
}
