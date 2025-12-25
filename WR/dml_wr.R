
# rm(list = ls())
library(SuperLearner)
library(glmnet)
library(tidyverse)

dml_wr_est <- function(df) {
  # calculate the estimate using the dml estimator
  ## nuisance functions are estimated using machine learning methods
  
  SL.library <- c("SL.glmnet", "SL.glm", "SL.randomForest", "SL.polymars", "SL.mean")
  
  nuisance <- function(df_train, df_test) {
    # glmnet for tp
    m_tp <- SuperLearner(Y = df_train[,'z'], X = df_train[,c('x1','x2','x3','x4')], 
                         family = binomial(), SL.library = SL.library, method = "method.NNloglik")
    m_tp_m <- SuperLearner(Y = df_train[,'z'], X = df_train[,c('x1t','x2t','x3t','x4t')], 
                           family = binomial(), SL.library = SL.library, method = "method.NNloglik")
    
    # glmnet for ps
    m_ps <- SuperLearner(Y = df_train[,'d'], X = df_train[,c('z','x1','x2','x3','x4')], 
                         family = binomial(), SL.library = SL.library, method = "method.NNloglik")
    m_ps_m <- SuperLearner(Y = df_train[,'d'], X = df_train[,c('z','x1t','x2t','x3t','x4t')], 
                           family = binomial(), SL.library = SL.library, method = "method.NNloglik")
    
    # glmnet for oc
    m_oc <- glmnet(y = df_train[,'y'], x = as.matrix(df_train[,c('z','d','x1','x2','x3','x4')]),
                   family = "multinomial")
    m_oc_m <- glmnet(y = df_train[,'y'], x = as.matrix(df_train[,c('z','d','x1t','x2t','x3t','x4t')]),
                     family = "multinomial")
    
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
    pmat11_pred <- predict(m_oc, as.matrix(x_oc11), type = "response")[,,length(m_oc$lambda)]
    pY1_11_est <- pmat11_pred[,1]
    pY2_11_est <- pmat11_pred[,2]
    pY3_11_est <- pmat11_pred[,3]
    
    x_oc00 <- cbind(0, cbind(0, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc00) <- c('z','d','x1','x2','x3','x4')
    pmat00_pred <- predict(m_oc, as.matrix(x_oc00), type = "response")[,,length(m_oc$lambda)]
    pY1_00_est <- pmat00_pred[,1]
    pY2_00_est <- pmat00_pred[,2]
    pY3_00_est <- pmat00_pred[,3]
    
    x_oc10 <- cbind(1, cbind(0, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc10) <- c('z','d','x1','x2','x3','x4')
    pmat10_pred <- predict(m_oc, as.matrix(x_oc10), type = "response")[,,length(m_oc$lambda)]
    pY1_10_est <- pmat10_pred[,1]
    pY2_10_est <- pmat10_pred[,2]
    pY3_10_est <- pmat10_pred[,3]
    
    x_oc01 <- cbind(0, cbind(1, df_test[,c('x1','x2','x3','x4')]))
    names(x_oc01) <- c('z','d','x1','x2','x3','x4')
    pmat01_pred <- predict(m_oc, as.matrix(x_oc01), type = "response")[,,length(m_oc$lambda)]
    pY1_01_est <- pmat01_pred[,1]
    pY2_01_est <- pmat01_pred[,2]
    pY3_01_est <- pmat01_pred[,3]
    
    x_oc11_m <- cbind(1, cbind(1, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc11_m) <- c('z','d','x1t','x2t','x3t','x4t')
    pmat11_m_pred <- predict(m_oc_m, as.matrix(x_oc11_m), type = "response")[,,length(m_oc_m$lambda)]
    pY1_11_m_est <- pmat11_m_pred[,1]
    pY2_11_m_est <- pmat11_m_pred[,2]
    pY3_11_m_est <- pmat11_m_pred[,3]
    
    x_oc00_m <- cbind(0, cbind(0, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc00_m) <- c('z','d','x1t','x2t','x3t','x4t')
    pmat00_m_pred <- predict(m_oc_m, as.matrix(x_oc00_m), type = "response")[,,length(m_oc_m$lambda)]
    pY1_00_m_est <- pmat00_m_pred[,1]
    pY2_00_m_est <- pmat00_m_pred[,2]
    pY3_00_m_est <- pmat00_m_pred[,3]
    
    x_oc10_m <- cbind(1, cbind(0, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc10_m) <- c('z','d','x1t','x2t','x3t','x4t')
    pmat10_m_pred <- predict(m_oc_m, as.matrix(x_oc10_m), type = "response")[,,length(m_oc_m$lambda)]
    pY1_10_m_est <- pmat10_m_pred[,1]
    pY2_10_m_est <- pmat10_m_pred[,2]
    pY3_10_m_est <- pmat10_m_pred[,3]
    
    x_oc01_m <- cbind(0, cbind(1, df_test[,c('x1t','x2t','x3t','x4t')]))
    names(x_oc01_m) <- c('z','d','x1t','x2t','x3t','x4t')
    pmat01_m_pred <- predict(m_oc_m, as.matrix(x_oc01_m), type = "response")[,,length(m_oc_m$lambda)]
    pY1_01_m_est <- pmat01_m_pred[,1]
    pY2_01_m_est <- pmat01_m_pred[,2]
    pY3_01_m_est <- pmat01_m_pred[,3]
    
    results <- list(tp_est, tp_m_est, ps1_est, ps0_est, ps1_m_est, ps0_m_est, 
                    pY1_11_est, pY2_11_est, pY3_11_est, pY1_00_est, pY2_00_est, pY3_00_est, 
                    pY1_10_est, pY2_10_est, pY3_10_est, pY1_01_est, pY2_01_est, pY3_01_est, 
                    pY1_11_m_est, pY2_11_m_est, pY3_11_m_est, pY1_00_m_est, pY2_00_m_est, pY3_00_m_est, 
                    pY1_10_m_est, pY2_10_m_est, pY3_10_m_est, pY1_01_m_est, pY2_01_m_est, pY3_01_m_est)
    
    names(results) <- c("tp_x", "tp_x_m", "ps1_x", "ps0_x", "ps1_x_m", "ps0_x_m", 
                        "pY1_11_x", "pY2_11_x", "pY3_11_x", "pY1_00_x", "pY2_00_x", "pY3_00_x", 
                        "pY1_10_x", "pY2_10_x", "pY3_10_x", "pY1_01_x", "pY2_01_x", "pY3_01_x", 
                        "pY1_11_x_m", "pY2_11_x_m", "pY3_11_x_m", "pY1_00_x_m", "pY2_00_x_m", "pY3_00_x_m", 
                        "pY1_10_x_m", "pY2_10_x_m", "pY3_10_x_m", "pY1_01_x_m", "pY2_01_x_m", "pY3_01_x_m")
    
    return(results)
  }
  
  df$tp_x <- NA
  df$tp_x_m <- NA
  df$ps1_x <- NA
  df$ps0_x <- NA
  df$ps1_x_m <- NA
  df$ps0_x_m <- NA
  df$pY1_11_x <- NA
  df$pY2_11_x <- NA
  df$pY3_11_x <- NA
  df$pY1_11_x_m <- NA
  df$pY2_11_x_m <- NA
  df$pY3_11_x_m <- NA
  df$pY1_00_x <- NA
  df$pY2_00_x <- NA
  df$pY3_00_x <- NA
  df$pY1_00_x_m <- NA
  df$pY2_00_x_m <- NA
  df$pY3_00_x_m <- NA
  df$pY1_10_x <- NA
  df$pY2_10_x <- NA
  df$pY3_10_x <- NA
  df$pY1_10_x_m <- NA
  df$pY2_10_x_m <- NA
  df$pY3_10_x_m <- NA
  df$pY1_01_x <- NA
  df$pY2_01_x <- NA
  df$pY3_01_x <- NA
  df$pY1_01_x_m <- NA
  df$pY2_01_x_m <- NA
  df$pY3_01_x_m <- NA
  
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
    df$pY1_11_x[I_list[[k]]] <- nuisance_results$pY1_11_x
    df$pY2_11_x[I_list[[k]]] <- nuisance_results$pY2_11_x
    df$pY3_11_x[I_list[[k]]] <- nuisance_results$pY3_11_x
    df$pY1_11_x_m[I_list[[k]]] <- nuisance_results$pY1_11_x_m
    df$pY2_11_x_m[I_list[[k]]] <- nuisance_results$pY2_11_x_m
    df$pY3_11_x_m[I_list[[k]]] <- nuisance_results$pY3_11_x_m
    df$pY1_00_x[I_list[[k]]] <- nuisance_results$pY1_00_x
    df$pY2_00_x[I_list[[k]]] <- nuisance_results$pY2_00_x
    df$pY3_00_x[I_list[[k]]] <- nuisance_results$pY3_00_x
    df$pY1_00_x_m[I_list[[k]]] <- nuisance_results$pY1_00_x_m
    df$pY2_00_x_m[I_list[[k]]] <- nuisance_results$pY2_00_x_m
    df$pY3_00_x_m[I_list[[k]]] <- nuisance_results$pY3_00_x_m
    df$pY1_10_x[I_list[[k]]] <- nuisance_results$pY1_10_x
    df$pY2_10_x[I_list[[k]]] <- nuisance_results$pY2_10_x
    df$pY3_10_x[I_list[[k]]] <- nuisance_results$pY3_10_x
    df$pY1_10_x_m[I_list[[k]]] <- nuisance_results$pY1_10_x_m
    df$pY2_10_x_m[I_list[[k]]] <- nuisance_results$pY2_10_x_m
    df$pY3_10_x_m[I_list[[k]]] <- nuisance_results$pY3_10_x_m
    df$pY1_01_x[I_list[[k]]] <- nuisance_results$pY1_01_x
    df$pY2_01_x[I_list[[k]]] <- nuisance_results$pY2_01_x
    df$pY3_01_x[I_list[[k]]] <- nuisance_results$pY3_01_x
    df$pY1_01_x_m[I_list[[k]]] <- nuisance_results$pY1_01_x_m
    df$pY2_01_x_m[I_list[[k]]] <- nuisance_results$pY2_01_x_m
    df$pY3_01_x_m[I_list[[k]]] <- nuisance_results$pY3_01_x_m
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
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ## different misspecification scenarios
  ### define functions
  tauN10 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i, 
                     pY1_11_x_i, pY2_11_x_i, pY3_11_x_i, 
                     pY1_00_x_i, pY2_00_x_i, pY3_00_x_i, 
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j, 
                     pY1_11_x_j, pY2_11_x_j, pY3_11_x_j, 
                     pY1_00_x_j, pY2_00_x_j, pY3_00_x_j) {
    
    mu10_g_x_ij <- pY2_11_x_i*pY1_00_x_j + pY3_11_x_i*(pY1_00_x_j+pY2_00_x_j)
    mu10_g_x_ji <- pY2_11_x_j*pY1_00_x_i + pY3_11_x_j*(pY1_00_x_i+pY2_00_x_i)
    mu10_l_x_ij <- pY1_11_x_i*(pY2_00_x_j+pY3_00_x_j) + pY2_11_x_i*pY3_00_x_j
    mu10_l_x_ji <- pY1_11_x_j*(pY2_00_x_i+pY3_00_x_i) + pY2_11_x_j*pY3_00_x_i
    
    g10_oioj_g <- 0.5*(z_i*d_i*(ps1_x_i-ps0_x_i)/(tp_x_i*ps1_x_i) *
                         (1-z_j)*(1-d_j)*(ps1_x_j-ps0_x_j)/((1-tp_x_j)*(1-ps0_x_j)) * 
                         (1*(y_i>y_j)-mu10_g_x_ij) + 
                         z_j*d_j*(ps1_x_j-ps0_x_j)/(tp_x_j*ps1_x_j) * 
                         (1-z_i)*(1-d_i)*(ps1_x_i-ps0_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j>y_i)-mu10_g_x_ji) +
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i)-((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i)) *
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j)-((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j)) *
                         mu10_g_x_ij +
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j)-((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j)) *
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i)-((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i)) *
                         mu10_g_x_ji)
    
    g10_oioj_l <- 0.5*(z_i*d_i*(ps1_x_i-ps0_x_i)/(tp_x_i*ps1_x_i) *
                         (1-z_j)*(1-d_j)*(ps1_x_j-ps0_x_j)/((1-tp_x_j)*(1-ps0_x_j)) * 
                         (1*(y_i<y_j)-mu10_l_x_ij) + 
                         z_j*d_j*(ps1_x_j-ps0_x_j)/(tp_x_j*ps1_x_j) * 
                         (1-z_i)*(1-d_i)*(ps1_x_i-ps0_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j<y_i)-mu10_l_x_ji) +
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i)-((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i)) *
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j)-((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j)) *
                         mu10_l_x_ij +
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j)-((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j)) *
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i)-((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i)) *
                         mu10_l_x_ji)
    
    return(c(mean(g10_oioj_g), mean(g10_oioj_l)))
  }
  
  tauN00 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i, 
                     pY1_10_x_i, pY2_10_x_i, pY3_10_x_i, 
                     pY1_00_x_i, pY2_00_x_i, pY3_00_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j, 
                     pY1_10_x_j, pY2_10_x_j, pY3_10_x_j, 
                     pY1_00_x_j, pY2_00_x_j, pY3_00_x_j) {
    
    mu00_g_x_ij <- pY2_10_x_i*pY1_00_x_j + pY3_10_x_i*(pY1_00_x_j+pY2_00_x_j)
    mu00_g_x_ji <- pY2_10_x_j*pY1_00_x_i + pY3_10_x_j*(pY1_00_x_i+pY2_00_x_i)
    mu00_l_x_ij <- pY1_10_x_i*(pY2_00_x_j+pY3_00_x_j) + pY2_10_x_i*pY3_00_x_j
    mu00_l_x_ji <- pY1_10_x_j*(pY2_00_x_i+pY3_00_x_i) + pY2_10_x_j*pY3_00_x_i
    
    g00_oioj_g <- 0.5*(z_i*(1-d_i)/tp_x_i * (1-z_j)*(1-d_j)*(1-ps1_x_j)/((1-tp_x_j)*(1-ps0_x_j)) *
                         (1*(y_i>y_j)-mu00_g_x_ij) + 
                         z_j*(1-d_j)/tp_x_j * (1-z_i)*(1-d_i)*(1-ps1_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j>y_i)-mu00_g_x_ji) +
                         (z_i*(ps1_x_i-d_i)/tp_x_i+1-ps1_x_i) * (z_j*(ps1_x_j-d_j)/tp_x_j+1-ps1_x_j) *
                         mu00_g_x_ij + 
                         (z_j*(ps1_x_j-d_j)/tp_x_j+1-ps1_x_j) * (z_i*(ps1_x_i-d_i)/tp_x_i+1-ps1_x_i) *
                         mu00_g_x_ji)
    
    g00_oioj_l <- 0.5*(z_i*(1-d_i)/tp_x_i * (1-z_j)*(1-d_j)*(1-ps1_x_j)/((1-tp_x_j)*(1-ps0_x_j)) *
                         (1*(y_i<y_j)-mu00_l_x_ij) + 
                         z_j*(1-d_j)/tp_x_j * (1-z_i)*(1-d_i)*(1-ps1_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j<y_i)-mu00_l_x_ji) +
                         (z_i*(ps1_x_i-d_i)/tp_x_i+1-ps1_x_i) * (z_j*(ps1_x_j-d_j)/tp_x_j+1-ps1_x_j) *
                         mu00_l_x_ij + 
                         (z_j*(ps1_x_j-d_j)/tp_x_j+1-ps1_x_j) * (z_i*(ps1_x_i-d_i)/tp_x_i+1-ps1_x_i) *
                         mu00_l_x_ji)
    
    return(c(mean(g00_oioj_g), mean(g00_oioj_l)))
  }
  
  tauN11 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i, 
                     pY1_11_x_i, pY2_11_x_i, pY3_11_x_i, 
                     pY1_01_x_i, pY2_01_x_i, pY3_01_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j, 
                     pY1_11_x_j, pY2_11_x_j, pY3_11_x_j, 
                     pY1_01_x_j, pY2_01_x_j, pY3_01_x_j) {
    
    mu11_g_x_ij <- pY2_11_x_i*pY1_01_x_j + pY3_11_x_i*(pY1_01_x_j+pY2_01_x_j)
    mu11_g_x_ji <- pY2_11_x_j*pY1_01_x_i + pY3_11_x_j*(pY1_01_x_i+pY2_01_x_i)
    mu11_l_x_ij <- pY1_11_x_i*(pY2_01_x_j+pY3_01_x_j) + pY2_11_x_i*pY3_01_x_j
    mu11_l_x_ji <- pY1_11_x_j*(pY2_01_x_i+pY3_01_x_i) + pY2_11_x_j*pY3_01_x_i
    
    g11_oioj_g <- 0.5*(z_i*d_i*ps0_x_i/(tp_x_i*ps1_x_i) * (1-z_j)*d_j/(1-tp_x_j) *
                         (1*(y_i>y_j)-mu11_g_x_ij) +
                         z_j*d_j*ps0_x_j/(tp_x_j*ps1_x_j) * (1-z_i)*d_i/(1-tp_x_i) *
                         (1*(y_j>y_i)-mu11_g_x_ji) +
                         ((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i) *
                         ((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j) * mu11_g_x_ij +
                         ((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j) *
                         ((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i) * mu11_g_x_ji)
    
    g11_oioj_l <- 0.5*(z_i*d_i*ps0_x_i/(tp_x_i*ps1_x_i) * (1-z_j)*d_j/(1-tp_x_j) *
                         (1*(y_i<y_j)-mu11_l_x_ij) +
                         z_j*d_j*ps0_x_j/(tp_x_j*ps1_x_j) * (1-z_i)*d_i/(1-tp_x_i) *
                         (1*(y_j<y_i)-mu11_l_x_ji) +
                         ((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i) *
                         ((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j) * mu11_l_x_ij +
                         ((1-z_j)*(d_j-ps0_x_j)/(1-tp_x_j)+ps0_x_j) *
                         ((1-z_i)*(d_i-ps0_x_i)/(1-tp_x_i)+ps0_x_i) * mu11_l_x_ji)
    
    return(c(mean(g11_oioj_g), mean(g11_oioj_l)))
  }
  
  
  ### tp-T, ps-T, oc-T
  tauN10_tpT_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpT_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpT_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-F, ps-T, oc-T
  tauN10_tpF_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpF_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpF_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-T, ps-F, oc-T
  tauN10_tpT_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpT_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpT_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-T, ps-T, oc-F
  tauN10_tpT_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpT_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpT_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  ### tp-F, ps-F, oc-T
  tauN10_tpF_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpF_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpF_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-F, ps-T, oc-F
  tauN10_tpF_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpF_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpF_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  ### tp-T, ps-F, oc-F
  tauN10_tpT_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpT_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpT_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  ### tp-F, ps-F, oc-F
  tauN10_tpF_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpF_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpF_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  # final estimates
  ## tp-T, ps-T, oc-T
  tau10_tpT_psT_ocT <- tauN10_tpT_psT_ocT[1]/tauN10_tpT_psT_ocT[2]
  tau00_tpT_psT_ocT <- tauN00_tpT_psT_ocT[1]/tauN00_tpT_psT_ocT[2]
  tau11_tpT_psT_ocT <- tauN11_tpT_psT_ocT[1]/tauN11_tpT_psT_ocT[2]
  
  ## tp-F, ps-T, oc-T
  tau10_tpF_psT_ocT <- tauN10_tpF_psT_ocT[1]/tauN10_tpF_psT_ocT[2]
  tau00_tpF_psT_ocT <- tauN00_tpF_psT_ocT[1]/tauN00_tpF_psT_ocT[2]
  tau11_tpF_psT_ocT <- tauN11_tpF_psT_ocT[1]/tauN11_tpF_psT_ocT[2]
  
  ## tp-T, ps-F, oc-T
  tau10_tpT_psF_ocT <- tauN10_tpT_psF_ocT[1]/tauN10_tpT_psF_ocT[2]
  tau00_tpT_psF_ocT <- tauN00_tpT_psF_ocT[1]/tauN00_tpT_psF_ocT[2]
  tau11_tpT_psF_ocT <- tauN11_tpT_psF_ocT[1]/tauN11_tpT_psF_ocT[2]
  
  ## tp-T, ps-T, oc-F
  tau10_tpT_psT_ocF <- tauN10_tpT_psT_ocF[1]/tauN10_tpT_psT_ocF[2]
  tau00_tpT_psT_ocF <- tauN00_tpT_psT_ocF[1]/tauN00_tpT_psT_ocF[2]
  tau11_tpT_psT_ocF <- tauN11_tpT_psT_ocF[1]/tauN11_tpT_psT_ocF[2]
  
  ## tp-F, ps-F, oc-T
  tau10_tpF_psF_ocT <- tauN10_tpF_psF_ocT[1]/tauN10_tpF_psF_ocT[2]
  tau00_tpF_psF_ocT <- tauN00_tpF_psF_ocT[1]/tauN00_tpF_psF_ocT[2]
  tau11_tpF_psF_ocT <- tauN11_tpF_psF_ocT[1]/tauN11_tpF_psF_ocT[2]
  
  ## tp-F, ps-T, oc-F
  tau10_tpF_psT_ocF <- tauN10_tpF_psT_ocF[1]/tauN10_tpF_psT_ocF[2]
  tau00_tpF_psT_ocF <- tauN00_tpF_psT_ocF[1]/tauN00_tpF_psT_ocF[2]
  tau11_tpF_psT_ocF <- tauN11_tpF_psT_ocF[1]/tauN11_tpF_psT_ocF[2]
  
  ## tp-T, ps-F, oc-F
  tau10_tpT_psF_ocF <- tauN10_tpT_psF_ocF[1]/tauN10_tpT_psF_ocF[2]
  tau00_tpT_psF_ocF <- tauN00_tpT_psF_ocF[1]/tauN00_tpT_psF_ocF[2]
  tau11_tpT_psF_ocF <- tauN11_tpT_psF_ocF[1]/tauN11_tpT_psF_ocF[2]
  
  ## tp-F, ps-F, oc-F
  tau10_tpF_psF_ocF <- tauN10_tpF_psF_ocF[1]/tauN10_tpF_psF_ocF[2]
  tau00_tpF_psF_ocF <- tauN00_tpF_psF_ocF[1]/tauN00_tpF_psF_ocF[2]
  tau11_tpF_psF_ocF <- tauN11_tpF_psF_ocF[1]/tauN11_tpF_psF_ocF[2]
  
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

