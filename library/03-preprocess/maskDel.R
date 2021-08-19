#!/usr/bin/env Rscript

###############################################################################
# Mask the Deletion that also exists in Control
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  query <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
  control <- read.table(args[2], sep = ",", header = FALSE, row.names = 1)
} else {
  query <- read.table("tmp_s", sep = ",", header = FALSE, row.names = 1)
  control <- read.table("tmp_c", sep = ",", header = FALSE, row.names = 1)
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

take_deletion_rate <- function(qc_freq) {
  qc_freq_ <- qc_freq
  del <- qc_freq_$x == "D"
  del_freq <- qc_freq_[del, ]$Freq.y / qc_freq_[del, ]$Freq.x
  return(del_freq)
}

qc_delfreq_subtract <- lapply(
  q_freq,
  function(x) take_deletion_rate(x)
  )


take_mis_rate <- function(idx, qc_freq) {
  qc_freq_ <- qc_freq[[idx]]
  qc_freq_ <- qc_freq_[qc_freq_$x != "D", ]
  freqx <- qc_freq_$Freq.x
  freqx <- freqx / (sum(freqx))
  qc_freq_$Freq.x <- freqx
  return(qc_freq_)
}

qc_misfreq <- lapply(
  seq(length(q_freq)),
  function(idx) take_mis_rate(idx, qc_freq)
  )

query_corrected <- query

for (idx in seq(ncol(query))) {
  tmp_delfreq <- qc_delfreq_subtract[[idx]]
  tmp_misfreq <- qc_misfreq[[idx]]

  if (tmp_delfreq < 1) {
    del_size <- sum(query[, idx] == "D")
    target_del_size <- as.integer(del_size * tmp_delfreq)

    tmp_embed <- sample(
      tmp_misfreq$x,
      size = target_del_size,
      replace = TRUE,
      prob = tmp_misfreq$Freq.x
      )
    tmp_embed <- as.character(tmp_embed)
    tmp_embed_D <- rep("D", del_size - target_del_size)
    query_corrected[, idx][query_corrected[, idx] == "D"] <- append(tmp_embed, tmp_embed_D)
  }
}

###############################################################################
# Output
###############################################################################

write.table(query_corrected,
  file = "", sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)
