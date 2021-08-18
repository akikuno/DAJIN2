#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  df <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
} else {
  df <- read.table(".DAJIN_temp/tmp_expansion_barcode31_albino.csv", sep = ",", header = FALSE, row.names = 1)
}

list_table <- apply(df, 2, table)

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