
# rm(list = ls())
library(SuperLearner)
library(tidyverse)

dml_mw_est <- function(df) {
  # calculate the estimate using the dml estimator
  ## nuisance functions are estimated using machine learning methods
  
  create.SL.knn <- function(k = c(20, 30)){ 
    for (mm in seq(length(k))){ 
      eval(parse(text = paste('SL.knn.',k[mm], '<-function(...,k=', k[mm], ')SL.knn(...,k=k)', sep='')),
           envir=.GlobalEnv) 
      } 
    invisible(TRUE) 
  } 
  
  create.SL.knn(c(20, 30, 40, 50, 60, 70))
  
  SL.library <- c("SL.glmnet", "SL.glm", "SL.randomForest", "SL.polymars", "SL.mean")
  
  SL.library.2 <- c("SL.glmnet", "SL.glm", "SL.randomForest", "SL.polymars", "SL.mean")
  
  nuisance <- function(df_train, df_test) {
    # glmnet for tp
    m_tp <- SuperLearner(Y = df_train[,'z'], X = df_train[,c('x1','x2','x3','x4')], verbose = F,
                         SL.library = SL.library, method = "method.NNloglik", family = binomial())
    m_tp_m <- SuperLearner(Y = df_train[,'z'], X = df_train[,c('x1t','x2t','x3t','x4t')], verbose = F,
                           SL.library = SL.library, method = "method.NNloglik", family = binomial())
    
    # glmnet for ps
    m_ps <- SuperLearner(Y = df_train[,'d'], X = df_train[,c('z','x1','x2','x3','x4')], verbose = F,
                         SL.library = SL.library, method = "method.NNloglik", family = binomial())
    m_ps_m <- SuperLearner(Y = df_train[,'d'], X = df_train[,c('z','x1t','x2t','x3t','x4t')], verbose = F,
                           SL.library = SL.library, method = "method.NNloglik", family = binomial())
    
    # glmnet for oc
    m_oc <- SuperLearner(Y = df_train[,'y'], X = df_train[,c('z','d','x1','x2','x3','x4')], verbose = F,
                         family = gaussian(), SL.library = SL.library.2, method = "method.NNLS")
    m_oc_m <- SuperLearner(Y = df_train[,'y'], X = df_train[,c('z','d','x1t','x2t','x3t','x4t')], verbose = F,
                           family = gaussian(), SL.library = SL.library.2, method = "method.NNLS")
    
    # output
    ## tp
    tp_pred <- predict(m_tp, df_test[,c('x1','x2','x3','x4')], onlySL = TRUE)
    tp_est <- tp_pred$pred
    
    tp_m_pred <- predict(m_tp_m, df_test[,c('x1t','x2t','x3t','x4t')], onlySL = TRUE)
    tp_m_est <- tp_m_pred$pred
    
    ## ps
    x_ps1 <- cbind(1, df_test[,c('x1','x2','x3','x4')])
    names(x_ps1) <- c('z','x1','x2','x3','x4')
    ps1_pred <- predict(m_ps, x_ps1, onlySL = TRUE)
    ps1_est <- ps1_pred$pred
    
    x_ps0 <- cbind(0, df_test[,c('x1','x2','x3','x4')])
    names(x_ps0) <- c('z','x1','x2','x3','x4')
    ps0_pred <- predict(m_ps, x_ps0, onlySL = TRUE)
    ps0_est <- ps0_pred$pred
    
    x_ps1_m <- cbind(1, df_test[,c('x1t','x2t','x3t','x4t')])
    names(x_ps1_m) <- c('z','x1t','x2t','x3t','x4t')
    ps1_m_pred <- predict(m_ps_m, x_ps1_m, onlySL = TRUE)
    ps1_m_est <- ps1_m_pred$pred
    
    x_ps0_m <- cbind(0, df_test[,c('x1t','x2t','x3t','x4t')])
    names(x_ps0_m) <- c('z','x1t','x2t','x3t','x4t')
    ps0_m_pred <- predict(m_ps_m, x_ps0_m, onlySL = TRUE)
    ps0_m_est <- ps0_m_pred$pred
    
    ## mu
    x_oc11 <- cbind(1, cbind(1, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc11) <- c('z','d','x1','x2','x3','x4')
    oc11_pred <- predict(m_oc, x_oc11, onlySL = TRUE)
    mu11_est <- oc11_pred$pred
    
    x_oc00 <- cbind(0, cbind(0, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc00) <- c('z','d','x1','x2','x3','x4')
    oc00_pred <- predict(m_oc, x_oc00, onlySL = TRUE)
    mu00_est <- oc00_pred$pred
    
    x_oc10 <- cbind(1, cbind(0, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc10) <- c('z','d','x1','x2','x3','x4')
    oc10_pred <- predict(m_oc, x_oc10, onlySL = TRUE)
    mu10_est <- oc10_pred$pred
    
    x_oc01 <- cbind(0, cbind(1, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc01) <- c('z','d','x1','x2','x3','x4')
    oc01_pred <- predict(m_oc, x_oc01, onlySL = TRUE)
    mu01_est <- oc01_pred$pred
    
    x_oc11_m <- cbind(1, cbind(1, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc11_m) <- c('z','d','x1t','x2t','x3t','x4t')
    oc11_m_pred <- predict(m_oc_m, x_oc11_m, onlySL = TRUE)
    mu11_m_est <- oc11_m_pred$pred
    
    x_oc00_m <- cbind(0, cbind(0, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc00_m) <- c('z','d','x1t','x2t','x3t','x4t')
    oc00_m_pred <- predict(m_oc_m, x_oc00_m, onlySL = TRUE)
    mu00_m_est <- oc00_m_pred$pred
    
    x_oc10_m <- cbind(1, cbind(0, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc10_m) <- c('z','d','x1t','x2t','x3t','x4t')
    oc10_m_pred <- predict(m_oc_m, x_oc10_m, onlySL = TRUE)
    mu10_m_est <- oc10_m_pred$pred
    
    x_oc01_m <- cbind(0, cbind(1, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc01_m) <- c('z','d','x1t','x2t','x3t','x4t')
    oc01_m_pred <- predict(m_oc_m, x_oc01_m, onlySL = TRUE)
    mu01_m_est <- oc01_m_pred$pred
    
    results <- list(tp_est, tp_m_est, ps1_est, ps0_est, ps1_m_est, ps0_m_est, 
                    mu11_est, mu00_est, mu10_est, mu01_est, mu11_m_est, mu00_m_est, 
                    mu10_m_est, mu01_m_est)
    names(results) <- c("tp_x", "tp_x_m", "ps1_x", "ps0_x", "ps1_x_m", "ps0_x_m", 
                        "mu11_x", "mu00_x", "mu10_x", "mu01_x", "mu11_x_m", "mu00_x_m", 
                        "mu10_x_m", "mu01_x_m")
    
    return(results)
  }
  
  df$tp_x <- NA
  df$tp_x_m <- NA
  df$ps1_x <- NA
  df$ps0_x <- NA
  df$ps1_x_m <- NA
  df$ps0_x_m <- NA
  df$mu11_x <- NA
  df$mu11_x_m <- NA
  df$mu00_x <- NA
  df$mu00_x_m <- NA
  df$mu10_x <- NA
  df$mu10_x_m <- NA
  df$mu01_x <- NA
  df$mu01_x_m <- NA
  
  # cross-fitting (K = 5)
  K <- 5
  I1 <- 1:(nrow(df)/K)
  I2 <- (nrow(df)/K + 1):(2*nrow(df)/K)
  I3 <- (2*nrow(df)/K + 1):(3*nrow(df)/K)
  I4 <- (3*nrow(df)/K + 1):(4*nrow(df)/K)
  I5 <- (4*nrow(df)/K + 1):nrow(df)
  
  I_list <- list(I1, I2, I3, I4, I5)
  
  for (k in 1:K) {
    nuisance_results <- nuisance(df[-I_list[[k]],], df[I_list[[k]],])
    df$tp_x[I_list[[k]]] <- nuisance_results$tp_x
    df$tp_x_m[I_list[[k]]] <- nuisance_results$tp_x_m
    df$ps1_x[I_list[[k]]] <- nuisance_results$ps1_x
    df$ps0_x[I_list[[k]]] <- nuisance_results$ps0_x
    df$ps1_x_m[I_list[[k]]] <- nuisance_results$ps1_x_m
    df$ps0_x_m[I_list[[k]]] <- nuisance_results$ps0_x_m
    df$mu11_x[I_list[[k]]] <- nuisance_results$mu11_x
    df$mu11_x_m[I_list[[k]]] <- nuisance_results$mu11_x_m
    df$mu00_x[I_list[[k]]] <- nuisance_results$mu00_x
    df$mu00_x_m[I_list[[k]]] <- nuisance_results$mu00_x_m
    df$mu10_x[I_list[[k]]] <- nuisance_results$mu10_x
    df$mu10_x_m[I_list[[k]]] <- nuisance_results$mu10_x_m
    df$mu01_x[I_list[[k]]] <- nuisance_results$mu01_x
    df$mu01_x_m[I_list[[k]]] <- nuisance_results$mu01_x_m
  }
  
  # calculate the denominators
  ## different misspecification scenarios
  ### tp-T, ps-T
  p1_tpT_psT <- mean(df$z*(df$d-df$ps1_x)/df$tp_x) + mean(df$ps1_x)
  p0_tpT_psT <- mean((1-df$z)*(df$d-df$ps0_x)/(1-df$tp_x)) + mean(df$ps0_x)
  
  den10_tpT_psT <- (p1_tpT_psT - p0_tpT_psT)**2
  den00_tpT_psT <- (1 - p1_tpT_psT)**2
  den11_tpT_psT <- p0_tpT_psT**2
  
  ### tp-T, ps-F
  p1_tpT_psF <- mean(df$z*(df$d-df$ps1_x_m)/df$tp_x) + mean(df$ps1_x_m)
  p0_tpT_psF <- mean((1-df$z)*(df$d-df$ps0_x_m)/(1-df$tp_x)) + mean(df$ps0_x_m)
  
  den10_tpT_psF <- (p1_tpT_psF - p0_tpT_psF)**2
  den00_tpT_psF <- (1 - p1_tpT_psF)**2
  den11_tpT_psF <- p0_tpT_psF**2
  
  ### tp-F, ps-T
  p1_tpF_psT <- mean(df$z*(df$d-df$ps1_x)/df$tp_x_m) + mean(df$ps1_x)
  p0_tpF_psT <- mean((1-df$z)*(df$d-df$ps0_x)/(1-df$tp_x_m)) + mean(df$ps0_x)
  
  den10_tpF_psT <- (p1_tpF_psT - p0_tpF_psT)**2
  den00_tpF_psT <- (1 - p1_tpF_psT)**2
  den11_tpF_psT <- p0_tpF_psT**2
  
  ### tp-F, ps-F
  p1_tpF_psF <- mean(df$z*(df$d-df$ps1_x_m)/df$tp_x_m) + mean(df$ps1_x_m)
  p0_tpF_psF <- mean((1-df$z)*(df$d-df$ps0_x_m)/(1-df$tp_x_m)) + mean(df$ps0_x_m)
  
  den10_tpF_psF <- (p1_tpF_psF - p0_tpF_psF)**2
  den00_tpF_psF <- (1 - p1_tpF_psF)**2
  den11_tpF_psF <- p0_tpF_psF**2
  
  # calculate the numerators
  m_oc_all <- SuperLearner(Y = df[,'y'], X = df[,c('z','d','x1','x2','x3','x4')], verbose = F,
                           family = gaussian(), SL.library = SL.library.2, method = "method.NNLS")
  m_oc_m_all <- SuperLearner(Y = df[,'y'], X = df[,c('z','d','x1t','x2t','x3t','x4t')], verbose = F,
                             family = gaussian(), SL.library = SL.library.2, method = "method.NNLS")
  
  sigma <- sd(df[,'y']-m_oc_all$SL.predict)
  sigma_m <- sd(df[,'y']-m_oc_m_all$SL.predict)
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ## different misspecification scenarios
  ### define functions
  tauN10 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i, mu11_x_i, mu00_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j, mu11_x_j, mu00_x_j, sig) {
    
    g10_oioj <- 0.5*(z_i*d_i*(ps1_x_i-ps0_x_i)/(tp_x_i*ps1_x_i) *
                       (1-z_j)*(1-d_j)*(ps1_x_j-ps0_x_j)/((1-tp_x_j)*(1-ps0_x_j)) * 
                       (1*(y_i>y_j)-pnorm((mu11_x_i-mu00_x_j)/(sqrt(2)*sig))) + 
                       z_j*d_j*(ps1_x_j-ps0_x_j)/(tp_x_j*ps1_x_j) * 
                       (1-z_i)*(1-d_i)*(ps1_x_i-ps0_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                       (1*(y_j>y_i)-pnorm((mu11_x_j-mu00_x_i)/(sqrt(2)*sig))) +
                       ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i)-((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i)) *
                       ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j)-((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j)) *
                       pnorm((mu11_x_i-mu00_x_j)/(sqrt(2)*sig)) + 
                       ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j)-((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j)) *
                       ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i)-((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i)) *
                       pnorm((mu11_x_j-mu00_x_i)/(sqrt(2)*sig)))
    
    return(mean(g10_oioj))
  }
  
  tauN00 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i, mu10_x_i, mu00_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j, mu10_x_j, mu00_x_j, sig) {
    
    g00_oioj <- 0.5*(z_i*(1-d_i)/tp_x_i * (1-z_j)*(1-d_j)*(1-ps1_x_j)/((1-tp_x_j)*(1-ps0_x_j)) * 
                       (1*(y_i>y_j)-pnorm((mu10_x_i-mu00_x_j)/(sqrt(2)*sig))) + 
                       z_j*(1-d_j)/tp_x_j * (1-z_i)*(1-d_i)*(1-ps1_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                       (1*(y_j>y_i)-pnorm((mu10_x_j-mu00_x_i)/(sqrt(2)*sig))) +
                       (z_i*(ps1_x_i-d_i)/tp_x_i+1-ps1_x_i) * (z_j*(ps1_x_j-d_j)/tp_x_j+1-ps1_x_j) *
                       pnorm((mu10_x_i-mu00_x_j)/(sqrt(2)*sig)) + 
                       (z_j*(ps1_x_j-d_j)/tp_x_j+1-ps1_x_j) * (z_i*(ps1_x_i-d_i)/tp_x_i+1-ps1_x_i) *
                       pnorm((mu10_x_j-mu00_x_i)/(sqrt(2)*sig)))
    
    return(mean(g00_oioj))
  }
  
  tauN11 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i, mu11_x_i, mu01_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j, mu11_x_j, mu01_x_j, sig) {
    
    g11_oioj <- 0.5*(z_i*d_i*ps0_x_i/(tp_x_i*ps1_x_i) * (1-z_j)*d_j/(1-tp_x_j) * 
                       (1*(y_i>y_j)-pnorm((mu11_x_i-mu01_x_j)/(sqrt(2)*sig))) + 
                       z_j*d_j*ps0_x_j/(tp_x_j*ps1_x_j) * (1-z_i)*d_i/(1-tp_x_i) *
                       (1*(y_j>y_i)-pnorm((mu11_x_j-mu01_x_i)/(sqrt(2)*sig))) +
                       ((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i) *
                       ((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j) *
                       pnorm((mu11_x_i-mu01_x_j)/(sqrt(2)*sig)) + 
                       ((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j) *
                       ((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i) *
                       pnorm((mu11_x_j-mu01_x_i)/(sqrt(2)*sig)))
    
    return(mean(g11_oioj))
  }
  
  
  ### tp-T, ps-T, oc-T
  tauN10_tpT_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpT_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpT_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-F, ps-T, oc-T
  tauN10_tpF_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpF_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpF_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-T, ps-F, oc-T
  tauN10_tpT_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpT_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpT_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-T, ps-T, oc-F
  tauN10_tpT_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpT_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpT_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  ### tp-F, ps-F, oc-T
  tauN10_tpF_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpF_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpF_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-F, ps-T, oc-F
  tauN10_tpF_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpF_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpF_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  ### tp-T, ps-F, oc-F
  tauN10_tpT_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpT_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpT_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  ### tp-F, ps-F, oc-F
  tauN10_tpF_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpF_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpF_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  # final estimates
  ## tp-T, ps-T, oc-T
  tau10_tpT_psT_ocT <- tauN10_tpT_psT_ocT/den10_tpT_psT
  tau00_tpT_psT_ocT <- tauN00_tpT_psT_ocT/den00_tpT_psT
  tau11_tpT_psT_ocT <- tauN11_tpT_psT_ocT/den11_tpT_psT
  
  ## tp-F, ps-T, oc-T
  tau10_tpF_psT_ocT <- tauN10_tpF_psT_ocT/den10_tpF_psT
  tau00_tpF_psT_ocT <- tauN00_tpF_psT_ocT/den00_tpF_psT
  tau11_tpF_psT_ocT <- tauN11_tpF_psT_ocT/den11_tpF_psT
  
  ## tp-T, ps-F, oc-T
  tau10_tpT_psF_ocT <- tauN10_tpT_psF_ocT/den10_tpT_psF
  tau00_tpT_psF_ocT <- tauN00_tpT_psF_ocT/den00_tpT_psF
  tau11_tpT_psF_ocT <- tauN11_tpT_psF_ocT/den11_tpT_psF
  
  ## tp-T, ps-T, oc-F
  tau10_tpT_psT_ocF <- tauN10_tpT_psT_ocF/den10_tpT_psT
  tau00_tpT_psT_ocF <- tauN00_tpT_psT_ocF/den00_tpT_psT
  tau11_tpT_psT_ocF <- tauN11_tpT_psT_ocF/den11_tpT_psT
  
  ## tp-F, ps-F, oc-T
  tau10_tpF_psF_ocT <- tauN10_tpF_psF_ocT/den10_tpF_psF
  tau00_tpF_psF_ocT <- tauN00_tpF_psF_ocT/den00_tpF_psF
  tau11_tpF_psF_ocT <- tauN11_tpF_psF_ocT/den11_tpF_psF
  
  ## tp-F, ps-T, oc-F
  tau10_tpF_psT_ocF <- tauN10_tpF_psT_ocF/den10_tpF_psT
  tau00_tpF_psT_ocF <- tauN00_tpF_psT_ocF/den00_tpF_psT
  tau11_tpF_psT_ocF <- tauN11_tpF_psT_ocF/den11_tpF_psT
  
  ## tp-T, ps-F, oc-F
  tau10_tpT_psF_ocF <- tauN10_tpT_psF_ocF/den10_tpT_psF
  tau00_tpT_psF_ocF <- tauN00_tpT_psF_ocF/den00_tpT_psF
  tau11_tpT_psF_ocF <- tauN11_tpT_psF_ocF/den11_tpT_psF
  
  ## tp-F, ps-F, oc-F
  tau10_tpF_psF_ocF <- tauN10_tpF_psF_ocF/den10_tpF_psF
  tau00_tpF_psF_ocF <- tauN00_tpF_psF_ocF/den00_tpF_psF
  tau11_tpF_psF_ocF <- tauN11_tpF_psF_ocF/den11_tpF_psF
  
  tp_sp <- c("T", "F", "T", "T", "F", "F", "T", "F")
  ps_sp <- c("T", "T", "F", "T", "F", "T", "F", "F")
  oc_sp <- c("T", "T", "T", "F", "T", "F", "F", "F")
  
  results <- data.frame(tp = tp_sp, ps = ps_sp, oc = oc_sp,
                        tau10 = c(tau10_tpT_psT_ocT, tau10_tpF_psT_ocT, tau10_tpT_psF_ocT,
                                  tau10_tpT_psT_ocF, tau10_tpF_psF_ocT, tau10_tpF_psT_ocF,
                                  tau10_tpT_psF_ocF, tau10_tpF_psF_ocF),
                        tau00 = c(tau00_tpT_psT_ocT, tau00_tpF_psT_ocT, tau00_tpT_psF_ocT,
                                  tau00_tpT_psT_ocF, tau00_tpF_psF_ocT, tau00_tpF_psT_ocF,
                                  tau00_tpT_psF_ocF, tau00_tpF_psF_ocF),
                        tau11 = c(tau11_tpT_psT_ocT, tau11_tpF_psT_ocT, tau11_tpT_psF_ocT,
                                  tau11_tpT_psT_ocF, tau11_tpF_psF_ocT, tau11_tpF_psT_ocF,
                                  tau11_tpT_psF_ocF, tau11_tpF_psF_ocF))
  
  return(results)
}

