#!/usr/bin/env Rscript

df <- read.table(file("stdin"), sep = ",", header = FALSE, row.names = 1)
# df <- read.table("tmp_maskMS.csv", sep = ",", header = FALSE, row.names = 1)

replacen <- function(vect) {
  tmp_table <- table(vect)
  tmp_table_omitN <- tmp_table[!grepl("N|n", names(tmp_table))]
  tmp_max <- which.max(tmp_table_omitN)
  tmp_names <- names(tmp_max)[1]
  if (length(tmp_names) == 0) {
    tmp_names <- "M"
  }
  vect[grepl("N|n", vect)] <- tmp_names
  return(vect)
}

df_replacen <- apply(df, 2, replacen)

# ? （ほかの部位のほうが異常にDが集積しているという判定がされるため, ）
# ? アルビノ点変異部位には効果がありませんでした.
# ?  一旦コメントアウトします（2021-08-06）
# # replace D
# # Dの頻度が”異常”かつDが最頻度ではない→最頻度に置換
# # それ以外→そのまま

# freqD <- function(vect) {
#   tmp_table <- table(vect)
#   tmp_table_D <- tmp_table[names(tmp_table) == "D"]
#   if (length(tmp_table_D) == 0) {
#     tmp_table_D <- 0
#   }
#   return(as.integer(tmp_table_D))
# }

# freq <- unlist(apply(df_replacen, 2, freqD))
# freq <- log(freq)
# hotelling <- (freq - mean(freq))^2/var(freq)
# hotelling_cols <- which(hotelling > qchisq(0.95, 1))

# replaceD <- function (vect) {
#   tmp_table <- table(vect)
#   tmp_max <- which.max(tmp_table)
#   tmp_names <- names(tmp_max)[1]
#   if (tmp_names != "D") {
#     vect[vect == "D"] <- tmp_names
#   }
#   return(vect)
# }

# if (length(hotelling_cols) > 0) {
#   df_replacen[, hotelling_cols] <- apply(df_replacen[, hotelling_cols], 2, replaceD)
# }

write.table(df_replacen, "", sep = ",", quote = FALSE, row.names = TRUE, col.names = FALSE)