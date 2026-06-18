library(dplyr); library(tidyr); library(readr); library(ggplot2)


INPUT_DIR="/data/original/TM"


Prepare_df <- function(df_name) {
  df = read_table(df_name, col_names = F)
  df = df[,c(4,8,9)]
  colnames(df) <- c("peak","log10P", "log10Q")
  return(df)
}

plot_pq_hist <- function(df, name,
                         pcol = "log10P", qcol = "log10Q",
                         bins = 50) {
  # backticks/all_of() handle the hyphenated column names safely
  
  
  long <- df %>%
    select(all_of(c(pcol, qcol))) %>%
    pivot_longer(cols = everything(),
                 names_to = "metric",
                 values_to = "value") %>%
    filter(!is.na(value))
  
  head(long)
  stats <- long %>%
    group_by(metric) %>%
    summarise(med = median(value),
              lo  = min(value),
              hi  = max(value),
              .groups = "drop") %>%
    mutate(label = sprintf("median = %.3g\nrange = [%.3g, %.3g]", med, lo, hi))
  head(stats)
  
  ggplot(long, aes(x = value)) +
    geom_histogram(bins = bins, boundary = 0,
                   fill = "#4C72B0", colour = "white", linewidth = 0.2) +
    # dashed line at the median of each facet
    geom_vline(data = stats, aes(xintercept = med),
               colour = "#C44E52", linetype = "dashed", linewidth = 0.5) +
    # median + range text, anchored to the top-right corner of each panel
    geom_text(data = stats, aes(x = Inf, y = Inf, label = label),
              hjust = 1.05, vjust = 1.4, size = 3, lineheight = 0.95) +
    facet_grid(cols = vars(metric), scales = "free_y",axes = "all") +
    labs(title = name, x = NULL, y = "Count") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90"),
          plot.title = element_text(face = "bold"))
}


# ---- function: one faceted histogram (P-value | Q-value) for a single dataframe ----

df1 <- Prepare_df(CASE_REP1)
df2 <- Prepare_df(CASE_REP2)
df3 <- Prepare_df(CASE_REP3)
df4 <- Prepare_df(CASE_REP4)
df5 <- Prepare_df(CASE_REP5_1)
df6 <- Prepare_df(CASE_REP5_2)

df1 <- Prepare_df(CTRL_REP1)
df2 <- Prepare_df(CTRL_REP2)
df3 <- Prepare_df(CTRL_REP3)
df4 <- Prepare_df(CTRL_REP4)
df5 <- Prepare_df(CTRL_REP5_1)
df6 <- Prepare_df(CTRL_REP5_2)

# ---- outside the function: combine 4-5 dataframes into ONE multi-page pdf ----
# Put the dataframes in a named list; names become the page titles
df_list <- list(
  TM_210 = df1,
  TM_213 = df2,
  TM_137 = df3,
  TM_155 = df4,
  TM_92_1 = df5,
  TM_92_2 = df6 
)

pdf("pq_histograms_CONTROL.pdf", width = 8, height = 4, onefile = TRUE)
for (nm in names(df_list)) {
  print(plot_pq_hist(df_list[[nm]], name = nm))
}
dev.off()
