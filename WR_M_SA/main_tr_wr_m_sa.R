
rm(list = ls())

source("data_gen_wr_m_sa.R")
source("triply_robust_wr_m_sa.R")

n_rep <- 1000

est_rec <- array(NA, dim = c(n_rep,4))

glm_warning_ind <- numeric(n_rep)

for (i in 1:n_rep) {
  
  set.seed(i)
  
  df <- data_gen_wr(2000)
  
  tryCatch(est <- triply_robust_wr_est(df), 
           warning = function(w) {
             print(w)
             glm_warning_ind[i] <<- 1
           })
  
  est_rec[i,] <- as.matrix(est)
  
  if (i %% 1 == 0) {
    cat("\r", paste("Repetition", i, sep=" "))
    flush.console()
  }
}

colMeans(est_rec)

est_rec_final <- est_rec[glm_warning_ind == 0 & apply(est_rec,1,function(x) sum(x<0)) == 0,]
colMeans(est_rec_final)


colnames(est_rec_final) <- c("tau10", "tau01", "tau00", "tau11")

write.csv(est_rec_final, "WR_results/tr_m_sa.csv", row.names = F)


