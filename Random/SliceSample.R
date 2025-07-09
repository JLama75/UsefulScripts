
#Subset 1000 rows from each label into a new dataframe
df_Subset <- df_all %>%
  group_by(label) %>%
  slice_sample(n = 1000) %>%
  ungroup()

#To group based on same values in column "GENE". Values in colum Trait will be collapsed using ",". Only the minimum P-Values in column 'P' will be reported
library(dplyr)

library(dplyr)

df_collapsed <- df %>%
  group_by(Gene) %>%
  summarise(
    Trait = paste(unique(Trait), collapse = ","),
    Trait2     = paste(unique(Trait2), collapse = ","),
    P     = min(P, na.rm = TRUE),
    .groups = "drop"
  )
