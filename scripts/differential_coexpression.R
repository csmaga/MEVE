# Analysis of differential co-expression between AP and WO 

library(csdR)

setwd("/scratch/crs12448/AP_WO_21/Diff_Coex")

# Load in normalized and transformed counts
norm_counts<-read.csv("AP_WO_normalized_counts.csv",row.names="X")

norm_counts_AP<-as.matrix(norm_counts[c(1,2,5,7,9,11,14,15,17,19,20,21),])
norm_counts_WO<-as.matrix(norm_counts[-c(1,2,5,7,9,11,14,15,17,19,20,21),])

csd_results <- run_csd(
  x_1 = norm_counts_AP, x_2 = norm_counts_WO,
  n_it = 10, nThreads = 4L, verbose = FALSE
)

write.csv(csd_results, "csd_results.csv")

pairs_to_pick <- 100L
c_filter <- partial_argsort(csd_results$cVal, pairs_to_pick)
c_frame <- csd_results[c_filter, ]
s_filter <- partial_argsort(csd_results$sVal, pairs_to_pick)
s_frame <- csd_results[s_filter, ]
d_filter <- partial_argsort(csd_results$dVal, pairs_to_pick)
d_frame <- csd_results[d_filter, ]

write.csv(c_frame, "c_frame.csv")
write.csv(s_frame, "s_frame.csv")
write.csv(d_frame, "d_frame.csv")






