#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

df_row <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
df_col <- read.table(args[2], sep = ",", header = FALSE, row.names = 1)

df <- df_row + df_col

write.table(df,
  file = "",
  sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)