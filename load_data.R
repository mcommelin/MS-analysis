#load all the available data from 2018 to 2020 in a tidy csv dataset
# Scheme_Catsod_data.pdf    table definitions
# overview_attributes.csv   units and description of attributes
# README.md                 general explanation of dataset

# Initialization -----------------
#' UPDATED Data UNTILL: 20210102

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(tmap)

# set locale(), standard timezone is Europe/Berline or UTC+1 - all time must be stored as UTC
# if a sensor is set to GMT+1, use 'Etc/GMT-1' as timezone, this is without DST.
nl_tz <- locale(tz = "Etc/GMT-1")

# Field units --------------------------------------------
#load shape files of fields - prepared from BRP files in: 'clip_BRP_to_cathment.R'
layers <- st_layers("input/maps/fields.gpkg")
layers <- layers$name
layers <- str_subset(layers, "fields")
years <- as.numeric(str_replace_all(layers, "fields", "20"))
df_list <- vector("list", length(layers))
gewas <- vector("list", length(layers))
for (i in seq_along(layers)) {
  df_list[[i]] <- st_read("input/maps/fields.gpkg", layer = layers[[i]])%>%
    mutate(year = years[[i]])
  df_list[[i]] <- df_list[[i]] %>%
    mutate(field_nr = (years[[i]] - 2017) * 1000 + 1:nrow(df_list[[i]])) %>%
    select(field_nr, year, CAT_GEWASCATEGORIE, GWS_GEWASCODE, GWS_GEWAS, geom)
  gewas[[i]] <- tibble(df_list[[i]]$GWS_GEWAS)
}
border <- st_read("input/maps/fields.gpkg", layer = "catch_outline")
#translate crop names
crops <- unique(bind_rows(gewas)) %>%
  mutate(crop_type = c("Wheat, winter-", "Maize, fodder-", "Barley, winter-",
                       "Conifers, open soil", "Grassland, permanent",
                       "Grass, temporary", "Potatoes, consumption", 
                       "Appels. Planted before current season", "Potatoes, seed NAK",
                       "Potatoes, seed TBM", "White cabbage, production",
                       "Beets, sugar-", "Appels. Planted in current season"))
names(crops) <- c("GWS_GEWAS", "crop_type")
#save fields as layer per year, with english crop names
field_owners <- read_csv("input/field_owners.csv") %>%
  select(-year)
fields_large <- tibble(large = c(rep("C", 9), rep("A", 5), rep("B", 7)),
                       field_nr = c(1009, 1020, 1013, 2024, 2019, 2023, 3024, 3012, 3022, 
                                    1021, 2017, 2009, 3014, 3015,
                                    1007, 1011, 2003, 2013, 3019, 3016, 3023))
append <- c(F, rep(T, length(layers) - 1))
for (i in seq_along(layers)) {
  df_list[[i]] <- df_list[[i]] %>%
    left_join(crops, by = "GWS_GEWAS") %>%
    left_join(field_owners, by = "field_nr") %>%
    left_join(fields_large, by = "field_nr") %>%
    mutate(area = round(st_area(geom), digits = 0))
  st_write(df_list[[i]], "data_closed/fields.gpkg", layers[i], append = F)
}
st_write(border, "data_closed/fields.gpkg", "catch_outline", append = F)
fields_all <- bind_rows(df_list)
st_write(fields_all, "data_closed/fields.gpkg", "fields_all", append = F)

# remove all gps and landowner data and store as csv
fields_ano <- fields_all %>%
  st_drop_geometry() %>%
  select(field_nr, year, large, crop_type, area)
write_csv(fields_ano, "data/fields.csv")

# Point samples -----------------------------------------------------------
#load all available gpx point data files. Combine to 1 gpkg and add a unique loc_id for each point.
pnts <- st_read("data/points.gpkg", layer = "points")
files <- dir("input/points/", pattern = "^.*\\.gpx$", full.names = T)
layers <- 'waypoints'
df_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df_list[[i]] <- st_read(files[[i]], layer = layers) %>%
    select(name, geometry)
}
pnts_new <- bind_rows(df_list)
# gpkg gives geometry column name 'geom' set this for all datasets you work with.
pnts_new <- st_sf(st_set_geometry(pnts_new, NULL), geom = st_geometry(pnts_new))
loc_id_max <- max(pnts$loc_id, na.rm = T)
pnts_new <- bind_rows(pnts, pnts_new) %>%
  distinct(geom, .keep_all = T) %>%
  filter(is.na(loc_id))
pnts_new <- pnts_new %>%
  mutate(loc_id = (1:nrow(pnts_new))+loc_id_max)
pnts <- bind_rows(pnts, pnts_new)
st_write(pnts, "data_closed/points.gpkg", "points", append = F)

#load TDR data and save to 1 table
files <- dir("input/points/", pattern = "^TDR.*\\.csv$", full.names = T)
dates <- ymd(str_extract(files, "(\\d+)"))
df_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df_list[[i]] <- read_csv(files[[i]], na = c("", "NA", "NULL"))
  nms <- str_subset(names(df_list[[i]]), "^TDR")
  df_list[[i]]<-  df_list[[i]] %>%
    pivot_longer(nms, names_to = "rep", values_to = 'VWC') %>%
    mutate(date = dates[[i]],
           rep = str_extract(rep, "\\d")) %>%
    left_join(pnts, by = "name") %>%
    select(-geom)
}
tdr <- bind_rows(df_list) %>%
  select(loc_id, date, sub_loc, rep, VWC, everything())

fields_all <- st_read("data/fields.gpkg", "fields_all")
pnts_rd <- st_transform(pnts, 28992) %>%
  select(loc_id)
years <- c(2018, 2019, 2020)
df_list <- vector("list", length = length(years))
for (i in seq_along(years)) {
  fields <- fields_all %>%
    filter(year == years[i]) %>%
    select(field_nr)
  tdr_i <- tdr %>%
    filter(year(date) == years[i]) %>%
    left_join(pnts_rd, by = 'loc_id') %>%
    st_as_sf()
  df_list[[i]] <- st_join(tdr_i, fields)
}
tdr <- bind_rows(df_list) %>%
  select(loc_id, date, sub_loc, rep, VWC, everything()) %>%
  st_drop_geometry() %>%
  tibble() %>%
  filter(!is.na(VWC)) %>%
  distinct()

write_csv(tdr, "data/TDR_VWC.csv")

#load KSAT data and save to 1 table
ksat <- read_csv("input/points/KSAT_2019.csv") %>%
  select(name:KSAT) %>%
  mutate(date = if_else(str_detect(name, "Punt"), "2019-11-05", "2019-08-29")) %>%
  left_join(pnts, by = 'name') %>%
  select(loc_id, sub_loc, date, everything(), -geom) %>%
  left_join(pnts_rd, by = 'loc_id') %>%
  st_as_sf()
fields <- fields_all %>%
  filter(year == 2019) %>%
  select(field_nr)
ksat <- tibble(st_join(ksat, fields)) %>%
  select(-geom) %>%
  distinct()

write_csv(ksat, "data/KSAT.csv")

# CR1000 data -------------------------------------------------------------
files <- dir("input/", pattern = "^CR1000.*", full.names = T)
df_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df_list[[i]] <- read_csv(files[[i]], na = c("", "NA", "NAN"), skip = 4, col_names = F, locale = nl_tz)
  names(df_list[[i]]) <- names(read_csv(files[[i]], skip = 1, n_max = 5))
  if (str_detect(files[[i]], "loc1")) {
      df_list[[i]] <-  mutate(df_list[[i]], location = 1)
    } else {
      df_list[[i]] <-  mutate(df_list[[i]], location = 2)
    }
}
cr1000 <- bind_rows(df_list)
names(cr1000) <- names(cr1000) %>%
  str_to_lower() %>%
  str_replace_all("([(])(.)([)])", "_\\2")
cr1000 <- select(cr1000, timestamp, location, rain_tot, everything()) %>%
  distinct()

# cr1000 rainfall
cr1000_rain <- cr1000 %>%
  select(-(pa_us_cs616_1:vwc_cs616_6)) %>%
  mutate(date = date(timestamp),
         loc_id = if_else(location == 1, 27, 30)) %>%
  select(timestamp, loc_id, everything(), -location)

#cr1000 VWC
pa_cols <- str_subset(names(cr1000), "pa_us")
vwc_cols <- str_subset(names(cr1000), "vwc")
cr_vwc <- cr1000 %>%
  select(-(rain_tot:batt_volt)) %>%
  pivot_longer(all_of(vwc_cols), names_to = "rep", values_to = "vwc",
               names_prefix = "vwc_cs616_", values_drop_na = T)%>%
  select(-(pa_us_cs616_1:pa_us_cs616_6))
cr_pa <- cr1000 %>%
  select(-(rain_tot:batt_volt)) %>%
  pivot_longer(all_of(pa_cols), names_to = "rep", values_to = "pa_us",
               names_prefix = "pa_us_cs616_", values_drop_na = T)%>%
  select(-(vwc_cs616_1:vwc_cs616_6))
cr1000_vars <- read_csv("input/vars_CR1000.csv", col_types = "dcccd")
cr1000_vwc <- left_join(cr_pa, cr_vwc, by = c('timestamp', 'location', 'rep')) %>%
  left_join(cr1000_vars, by = c('location', 'rep')) %>%
  filter(pa_us > 0.0001) %>%
  select(timestamp, loc_id, everything(), -(location:rep))

fields_all <- st_read("data/fields.gpkg", "fields_all")
pnts <- st_read("data/points.gpkg", layer = "points")
pnts_rd <- st_transform(pnts, 28992) %>%
  select(loc_id)
years <- c(2018, 2019, 2020)
df_list <- vector("list", length = length(years))
for (i in seq_along(years)) {
  fields <- fields_all %>%
    filter(year == years[i]) %>%
    select(field_nr)
  vwc_i <- cr1000_vwc %>%
    filter(year(timestamp) == years[i]) %>%
    left_join(pnts_rd, by = 'loc_id') %>%
    st_as_sf()
  df_list[[i]] <- st_join(vwc_i, fields)
}
cr1000_vwc <- bind_rows(df_list) %>%
  st_drop_geometry() %>%
  tibble() %>%
  distinct()

#save files
write_csv(cr1000_rain, "data/CR1000_rain.csv")
write_csv(cr1000_vwc, "data/CR1000_vwc.csv")

# Waterboard data 2014 - 2020 -------------------------------------------------------------
# time is GMT + 1
files <- dir("input/WB_raw/", pattern = "^WB.*\\.csv$", full.names = T)
nms <- c("timestamp", "Q", "Wh")
df_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df_list[[i]] <- read_csv(files[[i]], skip = 5, col_names = F, locale = nl_tz)
  names(df_list[[i]]) <- nms
}
wb_data <- bind_rows(df_list)
wb_data$timestamp <- dmy_hm(wb_data$timestamp, tz = "Etc/GMT-1")
options(pillar.sigfig = 5)
# load weekly data files for 2018 - present
files <- dir("input/WB_raw/", pattern = "p\\.csv$", full.names = T)
df_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df_list[[i]] <- read_csv(files[[i]], skip = 5, col_names = F, locale = nl_tz)
  names(df_list[[i]]) <- nms
}
df <- distinct(bind_rows(df_list))

# replace -999 with NA values
wb_data <- distinct(bind_rows(wb_data,df))
wb_data <- wb_data %>%
  mutate(Wh = replace(Wh, Wh == -999, NA),
         Q = replace(Q, Q == -999, NA))

#save file
write_csv(wb_data, "data/WB_data.csv")

# CR6 data ----------------------------------------------------------------
files <- dir("input/", pattern = "^C(R6|atsop).*\\.(csv|dat)$", full.names = T)
df_list <- vector("list", length(files))
for (i in seq_along(files)) {
  df_list[[i]] <- read_csv(files[[i]], skip = 4, col_names = F, locale = nl_tz)
  names(df_list[[i]]) <- names(read_csv(files[[i]], skip = 1, n_max = 5))
 }
#inspect timestamp colum to give correct format
df_list[[1]]$TIMESTAMP <- dmy_hm(df_list[[1]]$TIMESTAMP, tz = "Etc/GMT-1")
df_list[[4]]$TIMESTAMP <- mdy_hm(df_list[[4]]$TIMESTAMP, tz = "Etc/GMT-1")
  
cr6 <- bind_rows(df_list)
names(cr6) <- names(cr6) %>%
  str_to_lower() %>%
  str_replace_all("([(])(.)([)])", "_\\2")

cr6 <- select(cr6, timestamp, rain_tot, wat_level, level_mm, everything())
#remove duplicates
cr6 <- distinct(cr6)
# add date for dB relations
cr6 <- mutate(cr6, date = date(timestamp))

# One raw data file is made by merging (row_bind) all available data files together. 
# Because there is overlap. Identical entries are removed with ‘distinct()’.
# However still double timstamps occurred in the dataset indicating different entries 
# for some variables on the same timestamp – this is strange!!!
# After inspection, this can be explained by some rounding were the 5th digit sometimes differed 1 number. 
# Not clear how this could get into the data. Decided to remove with ‘!duplicated(timestamp)’
#
# remove duplicates which are distinct due to rounding error on 5th digit
cr6 <- cr6[ !duplicated(cr6$timestamp), ]

# check the timezones of CR6 data. Due to errors in datalogger setup timestamp can differ between 
# UTC+1 and UTC+2 because of DST. The WB_data is always UTC+1. So discharge peaks are compared to correct
# CR6 timestamps if needed. The timezone used in further calculations will be UTC+1.

# find timestamp where there are gaps in the dataset.
cr6_month <- cr6 %>%
  filter(date(timestamp) >= "2019-07-01" & date(timestamp) < "2019-10-01") %>%
  mutate(lag = as.numeric(timestamp - lag(timestamp))) %>%
  filter(lag >= 65) %>%
  select(timestamp, lag, everything())

# from the start of the data series up to 2019-07-12 all timestamps are UTC+1.
# on 2019-08-13, the time was set randomly, unclear why. On 2019-08-19, the time is corrected from 12:29 to 17:39 
# setting the time to UTC+2, by adding 5 hours and 9 minutes. To correct the time of this few days 4 hours and 9 minutes
# are added to the timestamps. This also matches the hydrograph on 2019-08-18 perfectly.
cr6_cor <- cr6 %>%
  mutate(timestamp = if_else(timestamp >= "2019-08-13 04:51" & timestamp <= "2019-08-19 11:29",
                             timestamp + hours(4) + minutes(9), timestamp))

#' from 2019-08-19 16:38 until present the cr6 timestamp is UTC+2, which will be corrected to UTC+1. During
#' analysis smaller anomalies in time are discovered, where wb_data and cr6 are out of sync a few minutes. 
#' This is not corrected in the raw data files, so must be taken into account during analysis!
cr6_cor <- cr6_cor %>%
  mutate(timestamp = if_else(timestamp >= "2019-08-19 16:39",
                             timestamp - hours(1), timestamp))

# check timezones by plotting hydrograph from cr6 and wb_data
ev_date <- "2020-10-03"
event <- cr6_cor %>%
  filter(date(timestamp) == ev_date) %>%
  mutate(Qcr6 = 1.428 * (wat_level)^1.55)
event_wb <- wb_data %>%
  filter(date(timestamp) == ev_date)
ggplot() +
  geom_line(data = event, aes(x = timestamp, y = Qcr6)) +
  geom_line(data = event_wb, aes(x = timestamp, y = Q), color = "blue")
#remove data before 2018-05-01, this is wrong due to energy drops of system, these points are collected in 2021-01.
cr6 <- cr6_cor %>%
  filter(date > "2018-05-01")

write_csv(cr6, "data/CR6.csv")

#check amount of data points per day
cr6 %>%
  count(day = floor_date(timestamp, "day")) %>%
  ggplot() +
  geom_line(aes(day, n))

# Pesticide applications --------------------------------------------------
#visually relate application fields to field units with following code:
fields_all <- st_read("data/fields.gpkg", "fields_all")
fields_all <- mutate(fields_all, area = st_area(fields_all))

#tmap_mode("view")

fields_map <- fields_all %>%
  filter(year == 2018)
#tm_shape(fields_map) + tm_fill('area')

# add field_nr to data
pest_appl_AB <- read_csv("input/Pest_application_AB.csv") %>%
  pivot_longer(c("field_nr", "field_nr_1", "field_nr_2"), names_to = "rep", values_to = "field_nr",
               values_drop_na = T) %>%
  select(date, field_nr, pesticide, dose, dose_unit)

pest_appl_C <- read_csv("input/Pest_application_C.csv") %>%
  rename(appl = weight, appl_unit = 'l/kg')

pest_appl <- bind_rows(pest_appl_AB, pest_appl_C)

write_csv(pest_appl, "data/Pest_application.csv")  

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

# EventSamples ------------------------------------------------------------
# timestamps are assigned to each sample see: Event_sample_selection.Rmd

#load csv sample sheet per year, and drop columns with calculations
# also adjust timestamp for 2019 and 2018 input, 2020 was already corrected in the input.
# from the start of the data series up to 2019-07-12 all timestamps are UTC+1.
# on 2019-08-13, the time was set randomly, unclear why. On 2019-08-19, the time is corrected at from 12:29 to 17:39 
# setting the time to UTC+2, by adding 5 hours and 9 minutes. To correct the time of this few days 4 hours and 9 minutes
# are added to the timestamps. This also matches the hydrograph on 2019-08-18 perfectly.

event_samples_2018 <- read_csv("input/Event_samples_2018.csv") %>%
  select(-(sed_g:vol_sed)) %>%
  mutate(timestamp = timestamp - hours(1)) %>%
  mutate(pest_samp = str_extract(notes, "^part.*")) %>%
  mutate(pest_samp = str_extract(pest_samp, "\\d+.\\d+"))
event_samples_2018 <- event_samples_2018 %>%
  mutate(an_ID = (1:nrow(event_samples_2018)+1000))

event_samples_2019 <- read_csv("input/Event_samples_2019.csv") %>%
  select(-(sed_g:vol_sed)) %>%
  mutate(timestamp = if_else(timestamp >= "2019-08-19 16:39",
                             timestamp - hours(2), timestamp - hours(1))) %>%
  mutate(timestamp = if_else(timestamp >= "2019-08-13 04:51" & timestamp <= "2019-08-19 11:29",
                             timestamp + hours(4) + minutes(9), timestamp))

#update names to uniform format
nms18 <- names(event_samples_2018)
nms18[4] <- "sample_code"
nms18[9] <- "cup_code"
names(event_samples_2018) <- nms18

nms19 <- names(event_samples_2019)
nms19[5:6] <- c("sample_code", "ISCO_bot_nr")
nms19[8] <- "cup_code"
nms19[12] <- "sample_nr"
names(event_samples_2019) <- nms19

#add columns to 2020
event_samples_2020 <- read_csv("input/Event_samples_2020.csv") %>%
  select(-(sed_g:vol_sed), -combine, - code) %>%
  mutate(ISCO = 3,
         sample_nr = as.numeric(str_extract(sample_code, "\\d*$"))) %>%
  arrange(timestamp)
event_samples_2020 <- event_samples_2020 %>%
  mutate(timestamp = if_else(date(timestamp) == "2020-11-16", timestamp - days(1), timestamp))

# combine datasets 2019 and 2020 and calculate sed_conc (2018 has a different handling pattern, include later
# because they are not needed for current analysis no further work is done for 2018)
rho_w <- 0.998
rho_sed <- 2.65
#adjust the TSS concentration for 2020-08-11 23:32, this is an outlier and most logical is a typing error in the raw data.
# The input sheet contains some anomalies due to errors in the lab work, following code adjust for that
event_samples <- bind_rows(event_samples_2019, event_samples_2020) %>%
  select(timestamp, ISCO, an_ID, everything()) %>%
  mutate(w_d_cup = if_else(w_d_cup == 266.32, 256.32, w_d_cup),
         w_sed = (w_d_cup - w_e_cup),
         vol_wat = if_else(is.na(w_d_ISCO), (w_f_ISCO - w_d_cup) / rho_w, 
                           (w_f_ISCO - w_d_ISCO) / rho_w)) %>%
  mutate(vol_wat = if_else(is.na(vol_wat), (w_f_cup - w_d_cup) / rho_w, vol_wat))
# update the 2 missing data samples - in 2019
mis_vals <- c(272.3, 271.36)
event_samples_na <- event_samples %>%
  filter(year(timestamp) == 2019, is.na(vol_wat)) %>%
  mutate(vol_wat = (mis_vals - w_d_cup) / rho_w)

#calculate sed_conc
event_samples <- event_samples %>%
  filter(!(year(timestamp) == 2019 & is.na(vol_wat))) %>%
  bind_rows(event_samples_na) %>%
  mutate(vol_sed = w_sed / rho_sed,
         sed_conc = if_else(w_sed < 0, 0, w_sed / ((vol_wat + vol_sed)/1000))) %>%
  select(-(w_sed:vol_sed))

# manually added pest_ID to this table.
event_samples_id <- read_csv("data/Event_samples.csv", lazy = FALSE) %>%
  select(pest_ID, sample_code)
event_samples <- left_join(event_samples, event_samples_id, by = "sample_code")

write_csv(event_samples,"data/Event_samples.csv")

#store sample analysis data
sample_analysis <- read_csv("input/Event_samples_analysis.csv") %>%
  select(an_ID, pH, OM, everything(), -(m.sed.v:OM.gr))
nms <- names(sample_analysis) %>%
  str_replace_all("(\\.)", "_")
nms[7:11] <- c("T_sample", "cup_code", "w_cup_dry", "w_cup_ig", "w_cup_empty")
names(sample_analysis) <- nms

sample_analysis <- mutate(sample_analysis, OM = (w_cup_dry - w_cup_ig) / (w_cup_dry - w_cup_empty))

write_csv(sample_analysis, "data/Event_samples_analysis.csv")

# Soil samples ------------------------------------------------------------
# added field numbers for soil samples 2019 based on visual comparison with field units
### fix.csv("input/Soil_samples.csv")
soil_samp <- read_csv("input/Soil_samples.csv")
tdr <- read_csv("data/TDR_VWC.csv") %>%
  select(loc_id, name, field_nr)
ss_point <- soil_samp %>%
  filter(is.na(field_nr)) %>%
  select(-field_nr, -field_nr_1, -field_nr_2) %>%
  left_join(tdr, by = "name") %>%
  distinct()
tdr <- tdr %>%
  select(-field_nr)
ss_field <- soil_samp %>%
  pivot_longer(c("field_nr", "field_nr_1", "field_nr_2"), names_to = "rep", values_to = "field_nr") %>%
  filter(!is.na(field_nr)) %>%
  left_join(tdr, by = 'name') %>%
  select(-rep) %>%
  distinct()
soil_samp <- bind_rows(ss_field, ss_point)
soil_samp <- left_join(soil_samp, ss_pest, by = c("name", "sub_loc")) %>%
  mutate(pest_ID = ifelse(date == "2019-09-30", NA, pest_ID))

write_csv(soil_samp, "data/Soil_samples.csv") 

#' add catchment sample analysis
catch_analysis <- read_csv("input/catchment_sample_analysis.csv")

write_csv(catch_analysis, "data/catchment_sample_analysis.csv")

# Texture -----------------------------------------------------------------

#'texture size distribution of sediment and soil samples was analyzed. 
#'Save with correct ID's to link to other tables.

tex_id <- read_csv("input/tex_range_id.csv")
tex_raw <- read_csv("input/textuur_analyse_01.csv") %>%
  mutate(tex_code = str_c(samp_name, "D", as.character(day(date))),
         tex_code = str_remove_all(tex_code, "-"),
         tex_code = str_replace(tex_code, "E", "S")) %>%
  select(-samp_name) %>%
  select(date, tex_code, everything())

write_csv(tex_id, "data/tex_range_id.csv")
write_csv(tex_raw, "data/texture_data.csv")

#' write code, to calculate the cumulative proportion for for each upper range.
#' calculate d50 and d90
#' this might belong in the general code.

# KNMI-data --------------------------------------------------------------
# the KNMI radar data is downloaded and clipped for Catsop in: 'KNMI_precipitation.R'
# timestamp is GMT + 0 / UTC
knmi_rain <- read_csv("input/KNMI_rain.csv", col_types = "Tcd", locale = default_locale()) %>%
  mutate(name = str_c("KNMI_", name))
pnts_new <- st_read("input/radar/pixels/KNMI_points.gpkg", "pixels") %>%
  mutate(name = str_c("KNMI_", as.character(1:6))) %>%
  select(-(row:y)) %>%
  st_transform(4326)
pnts <- st_read("data/points.gpkg", "points")
loc_id_max <- max(pnts$loc_id, na.rm = T)
pnts_new <- bind_rows(pnts, pnts_new) %>%
  distinct(geom, .keep_all = T) %>%
  filter(is.na(loc_id))
pnts_new <- pnts_new %>%
  mutate(loc_id = (1:nrow(pnts_new))+loc_id_max)
pnts <- bind_rows(pnts, pnts_new)
st_write(pnts, "data/points.gpkg", "points", append = F)

rain <- left_join(knmi_rain, pnts, by = 'name') %>%
  select(timestamp, loc_id, P, -name, - geom)
write_csv(rain, "data/KNMI_rain.csv")

#' long term temperature and rainfall for Maastricht Airport

knmi_ma_p <- read_csv("input/KNMI_P_MA.txt", skip = 13)
knmi_ma_t <- read_csv("input/KNMI_temp_MA.txt", skip = 13)
names(knmi_ma_p) <- tolower(names(knmi_ma_p))
names(knmi_ma_t) <- tolower(names(knmi_ma_t))
write_csv(knmi_ma_p, "data/knmi_ma_p.csv")
write_csv(knmi_ma_t, "data/knmi_ma_t.csv")

# Nutrient management ---------------
man_names <- read_csv("input/field_names_nutrients.csv", lazy = FALSE)
man_dat <- read_csv("input/nutrient_management.csv", lazy = FALSE) %>%
  filter(!is.na(field_name)) %>%
  select(-field_nr) %>%
  left_join(man_names, by = c("field_name", "year"))

write_csv(man_dat, "data/nutrient_management.csv")

