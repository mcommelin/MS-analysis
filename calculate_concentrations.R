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
calculate_recovery <- function(df_data, orig_alpha = 0.3) {

# select original samples of QC's
  origin_s <- df_data %>%
    select(sample_text, sample_conc, compound) %>%
    rename(origin = sample_text, samp_conc_o = sample_conc) %>%
    mutate(samp_conc_o = if_else(samp_conc_o < 0, 0, samp_conc_o))
# calculate recovery
  df_recov <- df_data %>%
    filter(an_type == "QC") %>%
    left_join(origin_s, by = c("origin", "compound")) %>%
    distinct(tqs_code, compound, .keep_all = T) %>%
    mutate(samp_conc_o = if_else(is.na(samp_conc_o), 0, samp_conc_o),
           recov = (sample_conc - samp_conc_o) / real_conc,
           orig_eff = samp_conc_o / real_conc,
           recov = if_else(orig_eff > orig_alpha, NaN, recov)) %>% # when original conc is >30% QC is not valid
    select(tqs_code, compound, recov, orig_eff)
  df_data <- df_data %>%
    left_join(df_recov, by = c("tqs_code", "compound"))
  
return(df_data)
}

# summary overview of the recovery
# # should matrix type be included?
summary_recovery <- function(df_data, digits = 2) {

  recov_summary <- df_data %>%
    filter(!is.na(recov)) %>%
    group_by(batch, compound) %>%
    summarize(mean = round(mean(recov), digits = digits),
           sd = round(sd(recov), digits = digits),
           n = n(),
           n_out = sum(recov < 0.8 | recov > 1.2))

return(recov_summary)
}

# correct sample concentrations for recovery
# recovery_correction <- function(df_data, recov_summary) {
#   d1 <- as.matrix(lc_multi[5:35])
#   v1 <- as.matrix(recov_stats[2])
#   
#   lc_matrix <- as_tibble(sweep(d1, 2, v1, FUN = "/")) %>%
#     bind_cols(lc_multi[1:4]) %>%
#     mutate(across(starts_with("conc_"), ~ if_else(. < -500, -999, .))) %>%
#   
#   df_data <- df_data %>%
#     mutate(cor_sample_conc = )
#   
# }