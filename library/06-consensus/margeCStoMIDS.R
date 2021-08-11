#!/usr/bin/env Rscript

###############################################################################
# Input files
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  df_cs <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
  df_mids <- read.table(args[2], sep = ",", header = FALSE, row.names = 1)
} else {
  df_cs <- read.table("tmp_cs", sep = ",", header = FALSE, row.names = 1)
  df_mids <- read.table("tmp_mids", sep = ",", header = FALSE, row.names = 1)
}
colnames(df_cs) <- seq(ncol(df_cs))
df_fa <- read.table(".DAJIN_temp/fasta/control.fa", skip = 1, header = FALSE)
c_fa <- strsplit(df_fa$V1, "")[[1]]

###############################################################################
# Extract mutation
###############################################################################

locIDS <- which(df_mids != "M")
df_mutIDS <- df_mids[locIDS, ]
df_csIDS <- df_cs[locIDS]

l_csmut <-
  lapply(seq_along(df_mutIDS), function(x) {
  mut <- df_mutIDS[[x]]
  cs <- df_csIDS[[x]]
  csmut <- cs[grepl(mut, cs)]
  return(csmut)
})
names(l_csmut) <- seq_along(df_mutIDS)
l_table <- lapply(l_csmut, table)
l_max <- (lapply(l_table, which.max))
l_names <- lapply(l_max, names)

v_mut <- unlist(l_names)
v_mut <- sub("[IS]", "", v_mut)

###############################################################################
# Embed mutation to FASTA
###############################################################################

c_fa[locIDS] <- v_mut

###############################################################################
# Output embedded FASTA
###############################################################################

write.table(
  cbind(df_mids, c_fa),
  file = "",
  sep = ",",
  col.names = FALSE,
  row.names = TRUE,
  quote = FALSE
)