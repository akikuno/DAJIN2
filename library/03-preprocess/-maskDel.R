#!/usr/bin/env Rscript

###############################################################################
# Controlにも存在するDeletionをマスクする
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

# idx <- 1260
# q_freq[[1260]]
# c_freq[[1260]]

concat_q_c <- function(idx, q_freq, c_freq) {
  df_q <- data.frame(q_freq[[idx]])
  df_c <- data.frame(c_freq[[idx]])
  df_qc <- merge(x = df_q, y = df_c, by = "x", all.x = TRUE)
  df_qc[is.na(df_qc)] <- 0
  return(df_qc)
}

qc_freq <- lapply(
  seq(length(q_freq)),
  function(idx) concat_q_c(idx, q_freq, c_freq)
  )

take_deletion_rate <- function(idx, qc_freq) {
  qc_freq_ <- qc_freq[[idx]]
  del <- qc_freq_$x == "D"
  del_freq <- qc_freq_[del, ]$Freq.y / qc_freq_[del, ]$Freq.x
  return(del_freq)
}

qc_delfreq_subtract <- lapply(
  seq(length(q_freq)),
  function(idx) take_deletion_rate(idx, qc_freq)
  )

# qc_delfreq_subtract[[1260]]

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
# idx <- 829
# query[, idx]
# prop.table(table(query[, idx]))
# prop.table(table(query_corrected[, idx]))

###############################################################################
# Output
###############################################################################

write.table(query_corrected,
  file = "", sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)

# after <- prop.table(table(query[, idx]))

# list_corrected <-
#   lapply(seq_along(query),
#     function(idx) fn_correction(idx, query, list_error))

# prop_del <- function(idx, list) {
#   tmp <- list[[idx]]
#   tmp[names(tmp) == "D"]
# }

# control_prop_del <- unlist(lapply(seq_along(list_control), function(idx) prop_del(idx, list_control)))
# query_prop_del <- unlist(lapply(seq_along(list_query), function(idx) prop_del(idx, list_query)))

# ###############################################################################
# # controlのほうがDeletionが多い場合、controlに合わせる。
# ###############################################################################

# gt <- control_prop_del > query_prop_del
# query_prop_del[gt] <- control_prop_del[gt]

# ###############################################################################
# # queryからcontrolを引く
# ###############################################################################

# summary(query_prop_del -  control_prop_del)
# which.max(query_prop_del -  control_prop_del)
# table(query[1258:1262])
# query_prop_del[1258:1262]
# control_prop_del[1258:1262]
# (query_prop_del -  control_prop_del)[1258:1262]
# # plot(control_prop_del[1000:1200])
# # plot(control_prop_del[2126:2165])
# # which.max(control_prop_del)
# # summary(control_prop_del)
# ###############################################################################
# # controlとqueryにD頻度をあてる
# ###############################################################################

# patch_del_score <- function(mids, prop) {
#   mids_D <- matrix(0, nrow(mids), ncol(mids))
#   rownames(mids_D) <- rownames(mids)
#   for (i in 1:ncol(mids)) {
#     mids_D[, i][mids[, i] == "D"] <- prop[i]
#   }
#   return(mids_D)
# }
# query_del_score <- patch_del_score(query, query_prop_del)
# control_del_score <- patch_del_score(control, control_prop_del)

# #?-----------------------------------------------------------------------------
# patch_del_score <- function(mids) {
#   mids_D <- matrix(0, nrow(mids), ncol(mids))
#   rownames(mids_D) <- rownames(mids)
#   for (i in 1:ncol(mids)) {
#     mids_D[, i][mids[, i] == "D"] <- 1
#   }
#   return(mids_D)
# }
# query_del_score <- patch_del_score(query)
# control_del_score <- patch_del_score(control)

# #?-----------------------------------------------------------------------------

# ###############################################################################
# # Cosine similarity
# ###############################################################################

# cossim <- function(v1, v2) {
#   sum(v1 * v2) / sqrt(sum(v1^2) * sum(v2^2))
# }

# tmp <- rep(0, nrow(control_del_score))
# i=4032
# for (i in 1:nrow(control_del_score)) {
#   for (j in 1:nrow(query_del_score)) {
#     tmp_zero <- query_del_score[j, ] + control_del_score[i, ] == 0
#     tmp_que <- query_del_score[j, ][!tmp_zero]
#     tmp_cont <- control_del_score[i, ][!tmp_zero]
#     tmp_cossim <- cossim(tmp_que, tmp_cont)
#   }
#   tmp[i] <- which.max(tmp_cossim)
# }

# i <- which.max(tmp)
# max(tmp)
# summary(tmp)

# par(mfrow=c(2,1))
# plot(query_del_score[1000,], ylim = c(0,1))
# plot(control_del_score[i,], ylim = c(0,1))
# dev.off()

# tmp_len <- length(control_prop_del)
# size <- 100
# window1 <- as.integer(seq(1, tmp_len, by = tmp_len / size))
# window2 <- window1 + as.integer(tmp_len / size) + 5
# # window2 <- as.integer(seq(tmp_len / size, tmp_len, by = tmp_len / size))
# if (tail(window2, 1) != tmp_len)
#   window2[size] <- tmp_len

# list_cossim <- lapply(1:size, function(idx) {
#   cossim(
#     control_prop_del[window1[idx] : window2[idx]],
#     query_prop_del[window1[idx] : window2[idx]]
#     )
# })
# v_cossim <- unlist(list_cossim)
# # sum(unlist(v_cossim) < 0.9)

# ###############################################################################
# # Deletion error correction
# ###############################################################################

# query_corrected <- query

# for (i in 1:100) {
#   w1 <- window1[i]
#   w2 <- window2[i]
#   tmp_samp <- query[w1:w2]
#   if (v_cossim[i] > 0.9) {
#     tmp_samp[tmp_samp == "D"] <- "M"
#     query_corrected[w1:w2] <- tmp_samp
#   } else {
#     #! コサイン類似度が低いときの対応はどうするか。引き算する？無視する？
#   }
# }


# #? idx = 3
# fn_correction <- function(idx, query, list_error){

#   df_x <- data.frame(x = query[[idx]], id = rownames(query))
#   df_merge <- merge(
#     x = df_x,
#     y = list_error[[idx]],
#     by = "x",
#     all.x = TRUE)
#   #? sum(df_merge$x == "D")
#   df_tmp <- unique(df_merge[,-2])
#   df_tmp <- df_tmp[df_tmp$sc_ratio != 1, ]

#   for (i in seq_along(df_tmp$x)) {
#     sc_ratio <- df_tmp[i, ]$sc_ratio
#     x <- df_tmp[i, ]$x
#     if ( sc_ratio == 0) {
#       df_merge$x[df_merge$x == x] <- "M"
#     } else {
#       df_merge$x[df_merge$x == x] <- query(
#         c("M", x),
#         size = sum(df_merge$x == x),
#         replace = TRUE,
#         prob = c(sc_ratio, 1 - sc_ratio)
#         )
#     }
#   }
#   return(df_merge$x)
# }


# list_corrected <-
#   lapply(seq_along(query),
#     function(idx) fn_correction(idx, query, list_error))

# df_corrected <- as.data.frame(do.call(cbind, list_corrected))
# rownames(df_corrected) <- rownames(query)

# gt <- control_prop_del_mut > query_prop_del_mut
# query_prop_del_mut[gt] <- control_prop_del_mut[gt]

# png("tmp.png")
# par(mfrow=c(3,2))
# plot(control_prop_del, ylim=c(0,1), main="cont")
# plot(query_prop_del, ylim=c(0,1), main="samp")
# plot(control_prop_del_mut, ylim=c(0,1), main="cont_mut")
# plot(query_prop_del_mut, ylim=c(0,1), main="samp_mut")
# # plot(control_prop_del[(2146-50):(2146+50)], ylim=c(0,1), main="cont")
# # plot(query_prop_del[(2146-50):(2146+50)], ylim=c(0,1), main="samp")
# # plot(control_prop_del[(829-50):(829-30)], ylim=c(0,1), main="cont")
# # plot(query_prop_del[(829-50):(829-30)], ylim=c(0,1), main="samp")
# dev.off()


# ccf_ <- function(v1, v2) {
#   tmp <- ccf(v1, v2)
#   tmp <- data.frame(acf = tmp$acf, lag = tmp$lag)
#   max_acf <- tmp[which.max(tmp$acf),]$acf
#   return(max_acf)
# }

# tmp_ccf <- lapply(1:100, function(idx) {
#   tmpw1 <- w1[idx]
#   tmpw2 <- w2[idx]
#   tmp <- ccf_(control_prop_del[tmpw1:tmpw2], query_prop_del[tmpw1:tmpw2])
#   return(tmp)
# })
# plot(unlist(tmp_ccf), ylim=c(0,1))
# sum(unlist(tmp_ccf) < 0.9)
# summary(unlist(tmp_ccf))

# fn_max <- function(idx) {
#   tmp <- list_control[[idx]]
#   max(tmp[names(tmp) != "M"])
# }

# for (i in seq_along(list_control)) {
#   tmp_max <- fn_max(i)
# }
# tmp_max <- lapply(seq_along(list_control), fn_max)
# tmp_max <- unlist(tmp_max)
# summary(tmp_max)
# sum(tmp_max > 0.1)
# png("tmp.png")
# plot(tmp_max)
# dev.off()

# fn_match <- function(idx) {
#   tmp <- list_control[[idx]]
#   tmp[names(tmp) == "M"]
# }

# fn_sub <- function(idx) {
#   tmp <- list_control[[idx]]
#   tmp[names(tmp) == "S"]
# }

# fn_ins <- function(idx) {
#   tmp <- list_control[[idx]]
#   max(tmp[grepl("[0-9]",names(tmp))],0)
# }
# tmp_match <- unlist(lapply(seq_along(list_control), fn_match))
# tmp_sub <- unlist(lapply(seq_along(list_control), fn_sub))
# del <- unlist(lapply(seq_along(list_control), fn_del))
# tmp_ins <- unlist(lapply(seq_along(list_control), fn_ins))
# summary(tmp_sub)
# summary(del)
# summary(tmp_ins)
# png("tmp.png")
# par(mfrow=c(2,2))
# plot(tmp_match, ylim=c(0,1))
# plot(tmp_sub, ylim=c(0,1))
# plot(del, ylim=c(0,1))
# plot(tmp_ins, ylim=c(0,1))
# dev.off()



# sum(tmp_match < 0.8)
# sum((tmp_max > 0.1) | (tmp_match < 0.8))
# summary(tmp_match[(tmp_max > 0.1) | (tmp_match < 0.8)])
# which((tmp_max > 0.1) | (tmp_match < 0.8))

# fn_error_rate <- function(idx, list_query, list_control) {
#   df_s <- as.data.frame(list_query[[idx]])
#   df_c <- as.data.frame(list_control[[idx]])
#   df_sc <- merge(x = df_s, y = df_c, by = "x", all = TRUE)
#   df_sc[is.na(df_sc)] <- (1/10)^10

#   tmp <- 1 / (df_sc$Freq.x)
#   v_sub <- 1 - (df_sc$Freq.y * tmp)
#   v_sub[v_sub > 0.9] <- 1
#   v_sub[v_sub < 0] <- 0
#   names(v_sub) <- df_sc$x
#   v_sub <- v_sub[!is.infinite(v_sub)]
#   df_sub <- data.frame(x = names(v_sub), sc_ratio = v_sub)

#   df_merge <- merge(x = df_s$x, y = df_sub, by = "x")
#   return(df_merge)
# }

# list_error <- lapply(seq_along(list_query),
#   function(idx) fn_error_rate(idx, list_query, list_control))
