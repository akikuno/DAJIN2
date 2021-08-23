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


# idx = 987 #?------------------------------------------------------------------------
# table(query[[idx]])#?------------------------------------------------------------------------
# idx = 829 #?------------------------------------------------------------------------
# table(control[[idx]])#?------------------------------------------------------------------------
idx = 829
subtract_mids <- function(idx, qc_freq) {
  tmp_freq <- qc_freq[[idx]]
  tmp_freq <- tmp_freq[tmp_freq$x != "M", ]
  tmp_subtract <- tmp_freq$Freq.y - tmp_freq$Freq.x
  # tmp_subtract[tmp_subtract < 0] <- 0
  # tmp_match <- 1 - sum(tmp_subtract)
  # tmp_subtract[tmp_freq$x == "M"] <- tmp_match
  # data.frame(x = tmp_freq$x, subtract = tmp_subtract)
  sum(tmp_subtract)
}

qc_subtract <- sapply(seq(length(q_freq)),
  function(idx) subtract_mids(idx, qc_freq)
  )

query_corrected <- query
for (idx in seq(ncol(query))) {
  tmp_freq <- qc_subtract[[idx]]
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
    tmp_embed_seq <- append(tmp_embed, tmp_embed_D)
    query_corrected[, idx][query[, idx] == "D"] <- sample(tmp_embed_seq)
  }
}

###############################################################################
# Output
###############################################################################

write.table(query_corrected,
  file = "", sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)
