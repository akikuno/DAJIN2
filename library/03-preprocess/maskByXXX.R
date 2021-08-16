#!/usr/bin/env Rscript

###############################################################################
# Controlにも存在するDeletionをマスクする
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  sample <- read.table(args[1], sep = ",", header = FALSE, row.names = 1)
  control <- read.table(args[2], sep = ",", header = FALSE, row.names = 1)
} else {
  sample <- read.table("tmp_s.csv", sep = ",", header = FALSE, row.names = 1)
  control <- read.table("tmp_c.csv", sep = ",", header = FALSE, row.names = 1)
}

fn_prop <- function(x) prop.table(table(x))
list_sample <- lapply(sample, fn_prop)
list_control <- lapply(control, fn_prop)
list_control
idx = 829 #???????????????????????????????
idx = 1377 #???????????????????????????????
sample[[idx]]

fn_max <- function(idx) {
  tmp <- list_control[[idx]]
  max(tmp[names(tmp) != "M"])
}

for (i in seq_along(list_control)) {
  tmp_max <- fn_max(i)
}
tmp_max <- lapply(seq_along(list_control), fn_max)
tmp_max <- unlist(tmp_max)
summary(tmp_max)
sum(tmp_max > 0.1)
png("tmp.png")
plot(tmp_max)
dev.off()

fn_match <- function(idx) {
  tmp <- list_control[[idx]]
  tmp[names(tmp) == "M"]
}

# fn_sub <- function(idx) {
#   tmp <- list_control[[idx]]
#   tmp[names(tmp) == "S"]
# }
fn_del <- function(idx, list) {
  tmp <- list[[idx]]
  tmp[names(tmp) == "D"]
}
# fn_ins <- function(idx) {
#   tmp <- list_control[[idx]]
#   max(tmp[grepl("[0-9]",names(tmp))],0)
# }
# tmp_match <- unlist(lapply(seq_along(list_control), fn_match))
# tmp_sub <- unlist(lapply(seq_along(list_control), fn_sub))
# tmp_del <- unlist(lapply(seq_along(list_control), fn_del))
# tmp_ins <- unlist(lapply(seq_along(list_control), fn_ins))
# summary(tmp_sub)
# summary(tmp_del)
# summary(tmp_ins)
# png("tmp.png")
# par(mfrow=c(2,2))
# plot(tmp_match, ylim=c(0,1))
# plot(tmp_sub, ylim=c(0,1))
# plot(tmp_del, ylim=c(0,1))
# plot(tmp_ins, ylim=c(0,1))
# dev.off()


tmp_del_control <- unlist(lapply(seq_along(list_control), function(idx) fn_del(idx, list_control)))
tmp_del_sample <- unlist(lapply(seq_along(list_sample), function(idx) fn_del(idx, list_sample)))
which.max(tmp_del)
png("tmp.png")
par(mfrow=c(2,2))
plot(tmp_del_control[(2146-50):(2146+50)], ylim=c(0,1), main="cont")
plot(tmp_del_sample[(2146-50):(2146+50)], ylim=c(0,1), main="samp")
plot(tmp_del_control[(829-50):(829-30)], ylim=c(0,1), main="cont")
plot(tmp_del_sample[(829-50):(829-30)], ylim=c(0,1), main="samp")
dev.off()

sum(tmp_match < 0.8)
sum((tmp_max > 0.1) | (tmp_match < 0.8))
summary(tmp_match[(tmp_max > 0.1) | (tmp_match < 0.8)])
which((tmp_max > 0.1) | (tmp_match < 0.8))

fn_error_rate <- function(idx, list_sample, list_control) {
  df_s <- as.data.frame(list_sample[[idx]])
  df_c <- as.data.frame(list_control[[idx]])
  df_sc <- merge(x = df_s, y = df_c, by = "x", all = TRUE)
  df_sc[is.na(df_sc)] <- (1/10)^10

  tmp <- 1 / (df_sc$Freq.x)
  v_sub <- 1 - (df_sc$Freq.y * tmp)
  v_sub[v_sub > 0.9] <- 1
  v_sub[v_sub < 0] <- 0
  names(v_sub) <- df_sc$x
  v_sub <- v_sub[!is.infinite(v_sub)]
  df_sub <- data.frame(x = names(v_sub), sc_ratio = v_sub)

  df_merge <- merge(x = df_s$x, y = df_sub, by = "x")
  return(df_merge)
}

list_error <- lapply(seq_along(list_sample),
  function(idx) fn_error_rate(idx, list_sample, list_control))

#? idx = 3
fn_correction <- function(idx, sample, list_error){

  df_x <- data.frame(x = sample[[idx]], id = rownames(sample))
  df_merge <- merge(
    x = df_x,
    y = list_error[[idx]],
    by = "x",
    all.x = TRUE)
  #? sum(df_merge$x == "D")
  df_tmp <- unique(df_merge[,-2])
  df_tmp <- df_tmp[df_tmp$sc_ratio != 1, ]

  for (i in seq_along(df_tmp$x)) {
    sc_ratio <- df_tmp[i, ]$sc_ratio
    x <- df_tmp[i, ]$x
    if ( sc_ratio == 0) {
      df_merge$x[df_merge$x == x] <- "M"
    } else {
      df_merge$x[df_merge$x == x] <- sample(
        c("M", x),
        size = sum(df_merge$x == x),
        replace = TRUE,
        prob = c(sc_ratio, 1 - sc_ratio)
        )
    }
  }
  return(df_merge$x)
}


list_corrected <-
  lapply(seq_along(sample),
    function(idx) fn_correction(idx, sample, list_error))

df_corrected <- as.data.frame(do.call(cbind, list_corrected))
rownames(df_corrected) <- rownames(sample)

write.table(df_corrected,
  file = "",
  sep = ",", col.names = FALSE, row.names = TRUE, quote = FALSE
)