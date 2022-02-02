#Analyse GC-MS, LC-MS and specific compound results

#LC calculation

#Process compounds - user specific?

# Ion ratio
data <- df_data[[i]] %>%
  mutate(IR = Area1/Area2)
