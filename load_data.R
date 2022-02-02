# load raw data form MS and GC machines and store is tidy table
# Initialization -----------------

library(tidyverse)

# Read header of colum names
# compile coumpound list
# make loop to load data for each compound
# Pesticide data ----------------------------------------------------------

## GLY & AMPA ----------------
# load or make compound list
compounds_gly <- c("GLY_qn", "GLY_ql", "GLY_IS", "AMPA_qn", "AMPA_ql", "AMPA_IS")
# load different analysis rounds and save to one file
batch <- dir("input/", pattern = "^GLY.*txt$", full.names = T)
batch_list <- vector("list", length = length(batch))
nm_file <- dir("input/", pattern = "^names_GLY.*", full.names = T)
for (j in seq_along(batch)) {
df_list <- vector("list", length = length(compounds_gly))
names_gly <- read_delim(nm_file[j], delim = "\t")
names(names_gly) <-c("number", "an_name")
n_samp <- max(names_gly$number)
for (i in seq_along(compounds_gly)) {
  skip = 3 + 3*i + (n_samp+1) * (i-1)
  df_list[[i]] <- read_delim(batch[j], delim = "\t", skip = skip, n_max = n_samp) %>%
    select(RT, Area)
  nms <- c(str_c("RT_", compounds_gly[i]), str_c("Area_", compounds_gly[i]))
  names(df_list[[i]]) <- nms
}
gly_batch <- read_delim(batch[j], delim = "\t", skip = 6, n_max = n_samp) %>%
  select('...1', Name)
names(gly_batch) <- c("number", "an_code")
gly_batch <- left_join(gly_batch, names_gly, by = "number")
batch_list[[j]] <- bind_cols(gly_batch, df_list) %>%
  mutate(week = str_extract(nm_file[j], "\\d+"))
}
gly_all <- bind_rows(batch_list) %>%
  mutate(an_type = if_else(an_name == "Blank", "blank", "sample"),
         an_type = if_else(an_name == "Sta 0.1 µg/ml", "point", an_type),
         an_type = if_else(str_detect(an_name, "Cal"), "curve", an_type),
         an_type = if_else(str_detect(an_name, "QC|XY"), "QC", an_type))
#add weight initial sample
gly_wn <- read_csv("input/GLY_weight_names.csv", skip = 3)
gly_all <- left_join(gly_all, gly_wn, by = "an_name")
## add ID variables to link pesticide data to point and event samples

#add sample type to pest data - to link to timeseries
gly_all <- gly_all %>%
  mutate(samp_type = str_extract(an_name, "_._"),
         samp_type = str_extract(samp_type, "[^_]")) 
# add pest_ID values + 10000 to distinct form an_ID

# same pest_ID for W and S sample of same timestamp
id <- rep(10001:10054, each = 2)
id <- id[1:107]
gly_ws <- gly_all %>%
  filter(str_detect(samp_type, "[WS]")) %>%
  mutate(date = ymd(str_extract(an_name, "^\\d*")),
         s_number = str_extract(an_name, "_\\d."),
         s_number = str_extract(s_number, "\\d+")) %>%
  arrange(date, s_number, samp_type) %>%
  mutate(pest_ID = id)
gly_c <- gly_all %>%
  filter(str_detect(samp_type, "[C]"))
gly_c <- mutate(gly_c, pest_ID = 1:nrow(gly_c) + max(gly_ws$pest_ID))

#make table to join pest_ID with event_samples
gly_s <- gly_ws %>%
  filter(str_detect(samp_type, "[S]")) %>%
  mutate(s_n2 = str_extract(an_name, "\\d+_\\d+"),
         s_n2 = str_extract(s_n2, "\\d+$"),
         s_n2 = if_else(is.na(s_n2), s_number, s_n2))
# think of some code to assign pest_id to event_sample. problem is range between s_number and s_n2
# no success yet, so did manually.
gly_temp <- gly_all %>%
  filter(is.na(samp_type))
gly_all <- bind_rows(gly_c, gly_ws, gly_temp) %>%
  arrange(week, number) %>%
  select(an_name, an_type, pest_ID, samp_type, everything())

# needed to link pest_ID to soil_samp
ss_pest <- select(gly_all, pest_ID, name, sub_loc)

gly_all <- select(gly_all, -(tube_code:sub_loc), -date, -s_number)

write_csv(gly_all, "data/LC_GLY.csv")

## multi-residue -----------------------------------------------------------
  # load different analysis rounds and save to one file
batch <- dir("input/", pattern = "^week.*", full.names = T)
batch_list <- vector("list", length = length(batch))
samples <- read_delim("input/week03.txt", delim = "\t", col_names = T, skip = 6, col_types = "ddccccc")
for (j in seq_along(batch)) {
  compounds_multi <- read_delim(batch[j], delim = "\t", col_names = F) %>%
    filter(str_detect(X1, "Compound \\d")) %>%
    mutate(X1 = str_replace_all(X1, " ", "_"))
  compounds <- str_replace_all(compounds_multi$X1, "Compound_\\d+:__|,", "")
  df_list <- vector("list", length = length(compounds))
  samples <- read_delim(batch[j], delim = "\t", col_names = T, skip = 6, col_types = "ddccccc")
  n_samp <- max(samples[ ,1], na.rm = T)
  for (i in seq_along(compounds)) {
    skip = 3 + 3*i + (n_samp+1) * (i-1)
    df_list[[i]] <- read_delim(batch[j], delim = "\t", skip = skip, n_max = n_samp, col_types = "ddcccc") %>%
      select(-(1:4))
    names(df_list[[i]]) <- c(str_c("RT_", compounds[i]), str_c("Area_", compounds[i]),
                             str_c("Area1_", compounds[i]))
  }
  batch_dat <- read_delim(batch[j], delim = "\t", skip = 6, n_max = n_samp) %>%
    select(2:4)
  names(batch_dat) <- c("number", "an_code", "an_name")
  batch_list[[j]] <- bind_cols(batch_dat, df_list) %>%
    mutate(week = str_extract(batch[j], "\\d+"))
}
#keep only samples of week03_water - remove duplicates in week03.
batch_list[[1]] <- anti_join(batch_list[[1]], batch_list[[2]], by = "an_code")

#combine all data, and add sample type, and analysis type for further calculations.
LC_all <- bind_rows(batch_list) %>%
  mutate(an_type = if_else(an_name == "Blank", "blank", "sample"),
         an_type = if_else(str_detect(an_name, "matrix|test"), "point", an_type),
         an_type = if_else(str_detect(an_name, "(S|s)olvent"), "point", an_type),
         an_type = if_else(str_detect(an_name, "ng/mL|B\\d Sta"), "curve", an_type),
         an_type = if_else(str_detect(an_name, "QC|BLANK|XY"), "QC", an_type),
         samp_type = str_extract(an_name, "_._"),
         samp_type = str_extract(samp_type, "[^_]")) %>%
  select(an_name, an_type, samp_type, week, everything()) 

#' add dilution ratio for batch 11 and 18
dil_ratio <- read_csv("input/dilution_batch12.csv")
LC_all <- LC_all %>%
  left_join(dil_ratio, by = "an_name") %>%
  mutate(dil_ratio = if_else(is.na(dil_ratio), 1, dil_ratio),
         dil_ratio = if_else(week == "03" | week == "04", 1, dil_ratio),
         dil_ratio = if_else(week == "18" & an_type == "sample",
                             as.numeric(str_extract(an_name, "\\d*$")), dil_ratio),
         an_name = if_else(week == "18" & an_type == "sample",
                           str_remove(an_name, "_\\d*$"), an_name))

# add aimed_weight and weight initial sample
aimed_w <- read_csv("input/LC_weight.csv") %>%
  rename(aimed_w_sample = aimed_w)
LC_all <- left_join(LC_all, aimed_w, by = "an_name") %>%
  select(an_name, an_type, samp_type, week, aimed_w_sample, everything())

#' fill NA for all samples that don't need dilution in batch12 + batch18 - 
#' MassLynx automatically finds values, but only diluted samples have correct values.
result_batches <- c("12", "18")
dil_batches <- vector("list", length = 2)
dil_batches[[1]] <- read_csv("input/dilute_samples_v2.csv")%>%
  filter(sum > 0) %>%
  rename_with(~str_replace(., "conc_", "")) %>%
  arrange(an_name)
# add 3 sample from batch 12 to batch18 due to LOQ repeat.
loq_18 <- tibble(an_name = c("20200820_C_B-LC", "20200820_C_D-LC", "20200924_S_1_9-LC")) %>%
  left_join(dil_batches[[1]], by = "an_name")
dil_batches[[2]]  <- read_csv("input/dilute_samples_batch18.csv") %>%
  filter(sum > 0) %>%
  full_join(loq_18) %>%
  select(an_name, all_of(str_subset(names(.), "qual_"))) %>%
  rename_with(~str_replace(., "qual_", "")) %>%
  arrange(an_name)
batch_na <- vector("list", length = 2)

for (k in seq_along(dil_batches)) {
  dil_batch <- dil_batches[[k]]
  result_batch <- filter(LC_all, week == result_batches[k]) %>%
    filter(an_type == "sample") %>%
    arrange(an_name, dil_ratio, samp_type) %>%
    semi_join(dil_batch, by = "an_name")
  compounds <- str_subset(names(dil_batches[[k]]), "^[A-Z]")
  batch_cols <- select(result_batch, all_of(str_subset(names(result_batch), "caffe|^[a-z]")))
  samples <- unique(batch_cols$an_name)
  comp_list <- vector("list", length = length(compounds))
  samp_list <- vector("list", length = length(samples))
  for (i in seq_along(compounds)) {
    dil_comp <- select(dil_batch, an_name, compounds[i])
    result_comp <- select(result_batch, an_name, str_subset(names(result_batch), str_c(compounds[i], "($|_.$)")))
    for (j in seq_along(samples)) {
      dil_samp <- filter(dil_comp, an_name == samples[j])[1,2]
      result_samp<- result_comp %>%
        filter(an_name == samples[j]) %>%
        select(-an_name)
      if (dil_samp[1,1] == 0 | is.na(dil_samp[1,1])) {
        result_samp[1:length(result_samp)] <- NA}
        samp_list[[j]] <- result_samp  
    }
  comp_list[[i]] <- bind_rows(samp_list)
  #check_list[[i]] <- sum(!is.na(comp_list[[i]]))/length(batch_samp)
  }
  batch_na[[k]] <- bind_cols(comp_list)%>%
    bind_cols(batch_cols)
  LC_all <- LC_all %>%
    anti_join(result_batch) %>%
    bind_rows(batch_na[[k]])
}

## add ID variables to link pesticide data to point and event samples
gly_all <- read_csv("data/LC_GLY.csv")
pest_ID_link <- gly_all %>%
  select(an_name, pest_ID) %>%
  filter(!is.na(pest_ID)) %>%
  mutate(an_name = str_replace(an_name, "GLY", "LC"))
LC_all <- LC_all %>%
  left_join(pest_ID_link, by = "an_name") %>%
  select(an_name, pest_ID, everything()) %>%
  rename_with(~str_replace(., "Oxathiopiprolin", "Oxathiapiprolin"))

write_csv(LC_all, "data/LC_multi.csv")
