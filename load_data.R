# load raw data form MS and GC machines and store is tidy table
# Initialization -----------------

library(tidyverse)

# make loop to load data for each compound
dir <-  "data_LC"
delim <- "\t"

load_raw_data <- function(dir, delim) {

# test datasets:
data_files <- dir(dir, pattern = ".*txt$", full.names = T)
df_data <- vector("list", length = length(data_files))

for( i in seq_along(data_files)) {
# compile coumpound list
compounds <- read_delim(data_files[i], delim = delim, col_names = F) %>%
  filter(str_detect(X1, "Compound \\d")) %>%
  mutate(X1 = str_replace_all(X1, "[ -]", "_"))
compounds <- str_replace_all(compounds$X1, "Compound_\\d+:__|,", "")
compounds

# Read header of column names
header <- readLines(data_files[i]) %>%
  str_detect("Name")
skip <- which(header == TRUE)[1] - 1

# Prepare and load the data into Area_DF
data <- read.delim(data_files[i], skip = skip, sep = delim) %>%
  mutate(X = as.numeric(X),
         compound = "na") %>%
  filter(!is.na(X)) %>%
  rename(number = X)
nms1 <- str_replace(names(data), "^.*Area.*$", "Area")
n_area <- sum(str_detect(nms1, "Area"))
nms1[nms1=="Area"]= str_c("Area", c(1:n_area))
names(data) <- nms1

# write compound name as variable in data
n_samp <- max(data$number)
for (j in seq_along(compounds)) {
data$compound[((j-1)*n_samp+1):(j*n_samp)]= compounds[j]
}
# add batch name
# add analysis type
# add matrix type
data <- data %>%
  mutate(batch = str_remove(data_files[i], "\\..*$"),
         batch = str_remove(batch, str_c("^", dir, "/")),
         an_type = if_else(str_detect(Sample.Text, "[Bb]lank|blk"), "blank", "sample"),
         an_type = if_else(str_detect(Sample.Text, "matrix|test|[Ss]olvent"), "point", an_type),
         an_type = if_else(str_detect(Sample.Text, "[Cc]al|ng/mL|B\\d"), "curve", an_type),
         an_type = if_else(str_detect(Sample.Text, "[Ss]ta|[Ss]td"), "standard", an_type),
         an_type = if_else(str_detect(Sample.Text, "QC|XY"), "QC", an_type))
df_data[[i]] <- data
}
df_data <- bind_rows(df_data)
return(df_data)
}

df_data <- load_raw_data("data_LC", "\t")


# add matrix type
matrix_types <- tibble(name = c("Water", "Sediment", "Soil"),
                       id = c("W", "S", "C"))
# steps in code MC during data load:
# - add initial/aimed weight for sediment/soil matrix samples - think about this
# - add an ID to link the data to a dataset - not needed to write code
# - remove all redundant data from dilution batches. - separate function for dilutions