
rm(list = ls())

source("data_gen_mw.R")
source("moment_pi+mu_mw.R")

n_rep <- 1000

est_rec <- array(NA, dim = c(8,3,n_rep))

glm_warning_ind <- numeric(n_rep)

for (i in 1:n_rep) {
  
  set.seed(i)
  
  df <- data_gen_mw(2000)
  
  tryCatch(est <- moment_pi_mu_mw_est(df), 
           warning = function(w) {
             print(w)
             glm_warning_ind[i] <<- 1
           })
  
  est_rec[,,i] <- as.matrix(est[,4:6])
  
  if (i %% 1 == 0) {
    cat("\r", paste("Repetition", i, sep=" "))
    flush.console()
  }
}

apply(est_rec, c(1,2), mean)

est_rec_final <- est_rec[,,glm_warning_ind == 0 & apply(est_rec,3,function(x) sum(x<0)) == 0]
apply(est_rec_final, c(1,2), mean)

col_names <- c("tpT_psT_ocT", "tpF_psT_ocT", "tpT_psF_ocT", "tpT_psT_ocF",
               "tpF_psF_ocT", "tpF_psT_ocF", "tpT_psF_ocF", "tpF_psF_ocF")

tau10_est <- data.frame(t(est_rec_final[,1,]))
tau00_est <- data.frame(t(est_rec_final[,2,]))
tau11_est <- data.frame(t(est_rec_final[,3,]))
colnames(tau10_est) <- col_names
colnames(tau00_est) <- col_names
colnames(tau11_est) <- col_names

write.csv(tau10_est, "MW_results/tau10_est_pi_mu.csv", row.names = F)
write.csv(tau00_est, "MW_results/tau00_est_pi_mu.csv", row.names = F)
write.csv(tau11_est, "MW_results/tau11_est_pi_mu.csv", row.names = F)



