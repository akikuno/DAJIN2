#!/usr/bin/env Rscript

###############################################################################
# Mask the Sequence error by subtraction from query to control
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  query <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
  control <- read.table(args[2], sep = ",", header = FALSE, row.names = 1)
} else {
  query <- read.table("tmp_query.csv", sep = ",", header = FALSE, row.names = 1)
  control <- read.table("tmp_control.csv", sep = ",", header = FALSE, row.names = 1)
}
colnames(query) <- seq(ncol(query))
colnames(control) <- seq(ncol(control))

###############################################################################
# MIDS frequency
###############################################################################

calculate_mids_frequency <- function(x) prop.table(table(x))

q_freq <- lapply(query, calculate_mids_frequency)
c_freq <- lapply(control, calculate_mids_frequency)

concat_qc <- function(idx, q_freq, c_freq) {
  df_q <- data.frame(q_freq[[idx]])
  df_c <- data.frame(c_freq[[idx]])
  df_qc <- merge(x = df_q, y = df_c, by = "x", all.x = TRUE)
  df_qc[is.na(df_qc)] <- 0
  return(df_qc)
}

qc_freq <- lapply(
  seq(length(q_freq)),
  function(idx) concat_qc(idx, q_freq, c_freq)
  )


detect_mutation <- function(idx, qc_freq) {
  tmp_freq <- qc_freq[[idx]]
  tmp_freq <- tmp_freq[tmp_freq$x != "M", ]
  tmp_freq$abssub <- abs(tmp_freq$Freq.y - tmp_freq$Freq.x)
  tmp_freq$mut <- tmp_freq$abssub > 0.1
  tmp_freq$mut <- tmp_freq$mut * !(tmp_freq$Freq.y > 0.1 & tmp_freq$Freq.x > 0.1)
  tmp_freq$prob <- tmp_freq$mut * tmp_freq$abssub
  tmp_freq$prob <- tmp_freq$prob / sum(tmp_freq$prob)
  tmp_return <- data.frame(x = tmp_freq$x, mut = tmp_freq$mut, prob = tmp_freq$prob)
  return(tmp_return)
}

qc_detect <- lapply(seq(length(q_freq)),
  function(idx) detect_mutation(idx, qc_freq)
  )

# idx <- 829
# idx <- 2024
# qc_detect[[idx]]
# tmp <- data.frame(x = c("S"), mut = 1, prob = 0.5)

query_corrected <- query
for (idx in seq(length(qc_detect))) {
  tmp_freq <- qc_detect[[idx]]
  if (sum(tmp_freq$mut) > 0) {
    for (mut in tmp_freq$x[tmp_freq$mut == FALSE])
      tmp_embed <- sample(tmp_freq$x,
        size = sum(query[, idx] == mut),
        replace = TRUE, prob = tmp_freq$prob)
      query_corrected[, idx][query[, idx] == mut] <- as.character(tmp_embed)
  } else {
    query_corrected[, idx] <- "M"
  }
}

###############################################################################
# Output
###############################################################################

write.table(query_corrected,
  file = "", sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)
