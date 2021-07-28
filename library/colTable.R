#!/usr/bin/env Rscript

df <- read.csv(file("stdin"), header = FALSE)

df_id <- df[, 1]
df <- df[, -1]
list_table <- apply(df, 2, table)

for (idx in seq(length(list_table))) {
  values <- names(list_table[[idx]])
  values <- as.integer(values)
  for (val in seq_along(values[values > 0])) {
    df[[idx]][df[[idx]] == val] <-
      list_table[[idx]][values == val]
  }
}

df <- cbind(df_id, df)

write.table(df, "",
  sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE
)