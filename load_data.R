# load raw data form MS and GC machines and store is tidy table
# Initialization -----------------

library(tidyverse)

dir <-  "data_LC"
delim <- "\t"
meta <-  FALSE
# function to load data from LC, GS or individual compound analysis
# beside the raw data also meta data is needed for the correct analysis of the samples
# Fill the raw data header names corresponding to the convention names for columns and make sure each area column has the word 'Area' in it.
tqs_code = "Name" #
sample_text = "Sample.Text" #
RT = "RT"

# function to load raw data.
load_raw_data <- function(dir, delim, meta) {

# find files in folder
data_files <- dir(dir, pattern = ".*txt$", full.names = T)
df_data <- vector("list", length = length(data_files))

for( i in seq_along(data_files)) {
# compile compound list
compounds <- read_delim(data_files[i], delim = delim, col_names = F) %>%
  filter(str_detect(X1, "Compound \\d")) %>%
  mutate(X1 = str_replace_all(X1, "[ -]", "_"))
compounds <- str_replace_all(compounds$X1, "Compound_\\d+:__|,", "")
#compounds

# Define how much lines to skip
header <- readLines(data_files[i]) %>%
  str_detect("Name")
skip <- which(header == TRUE)[1] - 1

# Prepare and load the data into df_data
data <- read.delim(data_files[i], skip = skip, sep = delim) %>%
  mutate(X = as.numeric(X)) %>%
  filter(!is.na(X)) %>%
  rename(number = X)
nms1 <- str_replace(names(data), "^.*Area.*$", "Area")
n_area <- sum(str_detect(nms1, "Area"))
nms1[nms1=="Area"]= str_c("Area", c(1:n_area))
names(data) <- nms1
n_samp <- max(data$number)
nms_data <- c(tqs_code, sample_text, RT, "Area1", "Area2")
nms_new <- c("tqs_code", "sample_text", "rt", "area1", "area2")
data <- data %>%
  select(all_of(nms_data)) %>%
  rename_with(., ~nms_new) %>%
  mutate(rt = as.numeric(rt),
         area1 = as.numeric(area1),
         area2 = as.numeric(area2),
         compound = "na",
         batch = str_remove(data_files[i], "\\..*$"),
         batch = str_remove(batch, str_c("^", dir, "/")))
  
# write compound name as variable in data

for (j in seq_along(compounds)) {
data$compound[((j-1)*n_samp+1):(j*n_samp)]= compounds[j]
}

if (meta == TRUE) {
# add batch name
# add analysis type
# add matrix type
data <- data %>%
  mutate(batch = str_remove(data_files[i], "\\..*$"),
         batch = str_remove(batch, str_c("^", dir, "/")),
         an_type = if_else(str_detect(sample_text, "[Bb]lank|blk"), "blank", "sample"),
         an_type = if_else(str_detect(sample_text, "matrix|test|[Ss]olvent"), "standard", an_type),
         an_type = if_else(str_detect(sample_text, "[Cc]al|ng/mL|B\\d"), "cal", an_type),
         an_type = if_else(str_detect(sample_text, "[Ss]ta|[Ss]td"), "standard", an_type),
         an_type = if_else(str_detect(sample_text, "QC|XY"), "QC", an_type))
}
df_data[[i]] <- data
}
df_data <- bind_rows(df_data)
return(df_data)
}

#test
df_data <- load_raw_data("data_LC", "\t", meta)

# to add - function to load meta data from separate file

# steps in code MC during data load:
# - add initial/aimed weight for sediment/soil matrix samples - think about this
# - add an ID to link the data to a dataset - not needed to write code
# - remove all redundant data from dilution batches. - separate function for dilutions