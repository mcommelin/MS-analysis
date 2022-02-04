#Analyse GC-MS, LC-MS and specific compound results

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

#Process compounds - user specific?

# 4. Concentrations in samples and QCs  ####
sample_concentration <- function(df_data, df_cal) {
  
lm_fit <- df_cal %>%
  select(batch, compound, intercept, slope) %>%
  distinct(.keep_all = TRUE)

df_data <- df_data %>%
  left_join(lm_fit, by = c("batch", "compound")) %>%
  mutate(detec_conc = area1 * slope + intercept,
         sample_conc = detec_conc * correction_factor * dilution_factor)
return(df_data)
}

# 5. QCs recovery ####




