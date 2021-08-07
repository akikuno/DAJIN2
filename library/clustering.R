################################################################################
# Install required packages
################################################################################

if (!requireNamespace("dbscan", quietly = TRUE)) {
  install.packages("dbscan", repos = "https://cloud.r-project.org/")
}

################################################################################
# Input
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  df_score <- read.csv(args[1], header = FALSE, stringsAsFactors = FALSE)
  pca_hotelling <- read.csv(args[2], header = FALSE)
  threads <- as.integer(args[3])
} else {
  df_score <- read.csv(
    ".DAJIN_temp/clustering/tmp_score.csv",
    header = FALSE, stringsAsFactors = FALSE
  )
  threads <- as.integer(parallel::detectCores() - 1)
}
df_score <- df_score[, colSums(df_score) != 0]
dev.new()
png("tmp.png")
plot(t(df_score[1, ]))
plot(t(df_score[2, ]))
plot(t(df_score[3, ]))
dev.off()
################################################################################
# PCA as a preprocessing step
################################################################################

pca <- prcomp(df_score, scale = FALSE)

pc_num <- 1:5
pc_score <- as.data.frame(pca$x[, pc_num])
prop_var <- pca$sdev[pc_num]^2 / sum(pca$sdev[pc_num]^2)
pc_score <- sweep(pc_score, 2, prop_var, FUN = "*")

###? DEBUG
# library("tidyverse")
# df_pca <- tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
# ggplot(df_pca, aes(x = PC1, y = PC2)) +
#   geom_point()

ggplot(pc_score, aes(x = seq_along(.data$PC1), y = PC1)) +
  geom_point()

# sum(pc_score$PC1 > 0)
# sum(pc_score$PC1 < 0)

# ggplot(as_tibble(t(df_score[1:2, ])), aes(x = seq_along(.data$"1"), y = .data$"1")) +
#   geom_point()
# ggplot(as_tibble(t(df_score[1:2, ])), aes(x = seq_along(.data$"1"), y = .data$"2")) +
#   geom_point()

# # tmp <- abs(df_score[1, ] - df_score[2, ])
# # which.max(tmp)
# df_score[2, 6051] %>% head()

####### DEBUG

################################################################################
# Clustering
################################################################################

# ===========================================================
# Extract cluster size with the largest cluster size
# and the most frequent cluster numbers
# ===========================================================

cl_seq <-
  as.integer(seq(nrow(pc_score) * 0.2, nrow(pc_score) * 0.4, length.out = 10)) + 2

hdb <- parallel::mclapply(cl_seq,
  function(.x) length(unique(dbscan::hdbscan(df_score, minPts = .x)$cluster)),
  mc.cores = threads
)
hdb <- unlist(setNames(hdb, cl_seq))
hdb <- hdb[hdb != 1]

if (length(hdb) > 0) {
  cl_mode <- names(which.max(table(hdb)))
  cl_size <- as.integer(max(names(hdb[hdb == cl_mode])))
} else {
  cl_size <- max(cl_seq)
}

cl_hdb <- dbscan::hdbscan(pc_score, minPts = cl_size)$cluster + 1

# ===========================================================
# Repeat clustering
# ===========================================================

for (.cl in unique(cl_hdb)) {
  .df_score <- df_score[cl_hdb == .cl, ]
  .pca <- prcomp(.df_score, scale = FALSE)
  if (length(pc_num) > ncol(.pca$x)) next
  .pc_score <- as.data.frame(.pca$x[, pc_num])
  .prop_var <- .pca$sdev[pc_num]^2 / sum(.pca$sdev[pc_num]^2)
  .pc_score <- sweep(.pc_score, 2, .prop_var, FUN = "*")
  .minPts <- as.integer(nrow(.pc_score) * 0.4) + 2L
  .cl_hdb <- dbscan::hdbscan(.pc_score, minPts = .minPts)$cluster + 1L
  cl_hdb[cl_hdb == .cl] <- .cl_hdb + (.cl * 100)
}
cl_hdb <- as.integer(as.factor(cl_hdb))

################################################################################
# Output
################################################################################

write(cl_hdb, stdout(), ncolumns = 1)