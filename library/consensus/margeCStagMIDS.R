#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  df <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
} else {
  df_cs <- read.table("tmp_cs", sep = ",", header = FALSE)
  df_mids <- read.table("tmp_mids", sep = ",", header = FALSE, row.names = 1)
}

df_mut <- df_mids[which(df_mids != "M"), ]
df_cs <- df_cs[which(df_mids != "M")]

list_table <- apply(df_cs, 2, table)

for (idx in seq(length(list_table))) {
  values <- names(list_table[[idx]])
  values <- as.integer(values)
  for (val in values[values > 0]) {
    df[[idx]][df[[idx]] == val] <-
      list_table[[idx]][values == val]
  }
}

write.table(df,
  file = "",
  sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)