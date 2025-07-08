
#Subset 1000 rows from each label into a new dataframe
df_Subset <- df_all %>%
  group_by(label) %>%
  slice_sample(n = 1000) %>%
  ungroup()
