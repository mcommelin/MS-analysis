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

#     -table of QCs and recovery 
#     -stability of recovery /Batch/ Compound
#     -stability of recovery / Compound over all batches
#     -calculate a column of corrected concentration C.matrix.cor

# needed: real concentration of QC samples
# reference to original sample if available to compensate
recovery <- function(df_data, orig_alpha = 0.3) {

# how to obtain real_conc and origin
  # from meta_data_add of with a new meta_add???

# select original samples of QC's
  origin_s <- df_data %>%
    select(an_name, samp_conc) %>%
    rename(origin = an_name, samp_conc_o = samp_conc) %>%
    mutate(samp_conc_o = if_else(samp_conc_o < 0, 0, samp_conc_o))
  QC_meta <- df_data %>%
    filter(an_type == "QC") %>%
    left_join(origin_s, by = "origin") %>%
    select(an_name, origin, samp_conc_o)
  
df_recov <- df_data %>%
  filter(an_type == "curve" | an_type == "QC" | an_type == "point") %>%
  left_join(QC_meta, by = "an_name") %>% # should include real_conc and samp_conc_o
  distinct(an_code, .keep_all = T) %>%
  mutate(recov = (samp_conc - samp_conc_o) / real_conc,
         orig_eff = samp_conc_o / real_conc,
         recov = if_else(orig_eff > orig_alpha, NaN, recov)) # when original conc is >30% QC is not valid
return(df_recov)
}