# Template file to run full MS analysis

# Initialization -----------------
library(tidyverse)
source("load_data.R")
source("calculate_concentrations.R")
source("loq_calculation.R")
source("calibration_calculation_visualsation.R")

# Load data ------------
# load the raw data files obtained from MassLynx
df_data <- load_raw_data("test2")

# add meta data
df_data <- meta_data_add("test1_meta.txt", df_data, delim = "/")
df_data <- meta_data_add("meta_test2.txt", df_data, delim = " ")

# Calibration ----------------------
# give name of internal standard
IS <- "13C_Caffeine"

# select calibration samples and fit linear models
df_cal <- calculate_calibration(df_data, IS, 1, 0.3, 0.3)

#plot linear fit for each compound
calibration_plot(df_cal)

#loq calculation
df_loq <- loq_calculation(df_cal)

# Sample concentrations ----------------------
df_conc <- sample_concentration(df_data, df_cal)
