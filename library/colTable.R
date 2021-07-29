#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

df <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)

list_table <- apply(df, 2, table)

for (idx in seq(length(list_table))) {
  values <- names(list_table[[idx]])
  values <- as.integer(values)
  for (val in seq_along(values[values > 0])) {
    df[[idx]][df[[idx]] == val] <-
      list_table[[idx]][values == val]
  }
}

write.table(df,
  file = "",
  sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)