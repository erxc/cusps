
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
    
    # glmnet for ps
    m_ps <- SuperLearner(Y = df_train[,'d'], X = df_train[,c('z','x1','x2','x3','x4')], 
                         family = binomial(), SL.library = SL.library, method = "method.NNloglik")
    
    # glmnet for oc
    m_oc <- glmnet(y = df_train[,'y'], x = as.matrix(df_train[,c('z','d','x1','x2','x3','x4')]),
                   family = "multinomial")
    
    # output
    ## tp
    tp_pred <- predict(m_tp, df_test[,c('x1','x2','x3','x4')], onlySL = TRUE)
    tp_est <- tp_pred$pred
    
    ## ps
    x_ps1 <- cbind(1, df_test[,c('x1','x2','x3','x4')])
    names(x_ps1) <- c('z','x1','x2','x3','x4')
    ps1_pred <- predict(m_ps, x_ps1, onlySL = TRUE)
    ps1_est <- ps1_pred$pred
    
    x_ps0 <- cbind(0, df_test[,c('x1','x2','x3','x4')])
    names(x_ps0) <- c('z','x1','x2','x3','x4')
    ps0_pred <- predict(m_ps, x_ps0, onlySL = TRUE)
    ps0_est <- ps0_pred$pred
    
    # sensitivity function
    eta_x <- 0.99*(1 - (ps1_est - ps0_est)/pmin(ps1_est, 1 - ps0_est))
    
    # principal scores combined with sensitivity function
    e10_x <- (ps1_est - ps0_est)/(1 - eta_x)
    e01_x <- eta_x*(ps1_est - ps0_est)/(1 - eta_x)
    e00_x <- 1 - ps0_est - (ps1_est - ps0_est)/(1 - eta_x)
    e11_x <- ps1_est - (ps1_est - ps0_est)/(1 - eta_x)
    
    
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
    
    
    results <- list(tp_est, ps1_est, ps0_est, eta_x, e10_x, e01_x, e00_x, e11_x,
                    pY1_11_est, pY2_11_est, pY3_11_est, pY1_00_est, pY2_00_est, pY3_00_est, 
                    pY1_10_est, pY2_10_est, pY3_10_est, pY1_01_est, pY2_01_est, pY3_01_est)
    
    names(results) <- c("tp_x", "ps1_x", "ps0_x", "eta_x", "e10_x", "e01_x", "e00_x", "e11_x",
                        "pY1_11_x", "pY2_11_x", "pY3_11_x", "pY1_00_x", "pY2_00_x", "pY3_00_x", 
                        "pY1_10_x", "pY2_10_x", "pY3_10_x", "pY1_01_x", "pY2_01_x", "pY3_01_x")
    
    return(results)
  }
  
  df$tp_x <- NA
  df$ps1_x <- NA
  df$ps0_x <- NA
  df$eta_x <- NA
  df$e10_x <- NA
  df$e01_x <- NA
  df$e00_x <- NA
  df$e11_x <- NA
  df$pY1_11_x <- NA
  df$pY2_11_x <- NA
  df$pY3_11_x <- NA
  df$pY1_00_x <- NA
  df$pY2_00_x <- NA
  df$pY3_00_x <- NA
  df$pY1_10_x <- NA
  df$pY2_10_x <- NA
  df$pY3_10_x <- NA
  df$pY1_01_x <- NA
  df$pY2_01_x <- NA
  df$pY3_01_x <- NA
  
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
    df$ps1_x[I_list[[k]]] <- nuisance_results$ps1_x
    df$ps0_x[I_list[[k]]] <- nuisance_results$ps0_x
    df$eta_x[I_list[[k]]] <- nuisance_results$eta_x
    df$e10_x[I_list[[k]]] <- nuisance_results$e10_x
    df$e01_x[I_list[[k]]] <- nuisance_results$e01_x
    df$e00_x[I_list[[k]]] <- nuisance_results$e00_x
    df$e11_x[I_list[[k]]] <- nuisance_results$e11_x
    df$pY1_11_x[I_list[[k]]] <- nuisance_results$pY1_11_x
    df$pY2_11_x[I_list[[k]]] <- nuisance_results$pY2_11_x
    df$pY3_11_x[I_list[[k]]] <- nuisance_results$pY3_11_x
    df$pY1_00_x[I_list[[k]]] <- nuisance_results$pY1_00_x
    df$pY2_00_x[I_list[[k]]] <- nuisance_results$pY2_00_x
    df$pY3_00_x[I_list[[k]]] <- nuisance_results$pY3_00_x
    df$pY1_10_x[I_list[[k]]] <- nuisance_results$pY1_10_x
    df$pY2_10_x[I_list[[k]]] <- nuisance_results$pY2_10_x
    df$pY3_10_x[I_list[[k]]] <- nuisance_results$pY3_10_x
    df$pY1_01_x[I_list[[k]]] <- nuisance_results$pY1_01_x
    df$pY2_01_x[I_list[[k]]] <- nuisance_results$pY2_01_x
    df$pY3_01_x[I_list[[k]]] <- nuisance_results$pY3_01_x
  }
  
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ### define functions
  tauN10 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, e10_x_i, tp_x_i, eta_x_i,
                     pY1_11_x_i, pY2_11_x_i, pY3_11_x_i, 
                     pY1_00_x_i, pY2_00_x_i, pY3_00_x_i, 
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, e10_x_j, tp_x_j, eta_x_j,
                     pY1_11_x_j, pY2_11_x_j, pY3_11_x_j, 
                     pY1_00_x_j, pY2_00_x_j, pY3_00_x_j) {
    
    mu10_g_x_ij <- pY2_11_x_i*pY1_00_x_j + pY3_11_x_i*(pY1_00_x_j+pY2_00_x_j)
    mu10_g_x_ji <- pY2_11_x_j*pY1_00_x_i + pY3_11_x_j*(pY1_00_x_i+pY2_00_x_i)
    mu10_l_x_ij <- pY1_11_x_i*(pY2_00_x_j+pY3_00_x_j) + pY2_11_x_i*pY3_00_x_j
    mu10_l_x_ji <- pY1_11_x_j*(pY2_00_x_i+pY3_00_x_i) + pY2_11_x_j*pY3_00_x_i
    
    g10_oioj_g <- 0.5*(z_i*d_i*e10_x_i/(tp_x_i*ps1_x_i) *
                         (1-z_j)*(1-d_j)*e10_x_j/((1-tp_x_j)*(1-ps0_x_j)) * 
                         (1*(y_i>y_j)-mu10_g_x_ij) + 
                         z_j*d_j*e10_x_j/(tp_x_j*ps1_x_j) * 
                         (1-z_i)*(1-d_i)*e10_x_i/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j>y_i)-mu10_g_x_ji) +
                         ((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         ((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         mu10_g_x_ij +
                         ((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         ((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         mu10_g_x_ji)
    
    g10_oioj_l <- 0.5*(z_i*d_i*e10_x_i/(tp_x_i*ps1_x_i) *
                         (1-z_j)*(1-d_j)*e10_x_j/((1-tp_x_j)*(1-ps0_x_j)) * 
                         (1*(y_i<y_j)-mu10_l_x_ij) + 
                         z_j*d_j*e10_x_j/(tp_x_j*ps1_x_j) * 
                         (1-z_i)*(1-d_i)*e10_x_i/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j<y_i)-mu10_l_x_ji) +
                         ((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         ((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         mu10_l_x_ij +
                         ((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         ((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         mu10_l_x_ji)
    
    return(c(mean(g10_oioj_g), mean(g10_oioj_l)))
  }
  
  tauN01 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, e01_x_i, tp_x_i, eta_x_i,
                     pY1_10_x_i, pY2_10_x_i, pY3_10_x_i, 
                     pY1_01_x_i, pY2_01_x_i, pY3_01_x_i, 
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, e01_x_j, tp_x_j, eta_x_j,
                     pY1_10_x_j, pY2_10_x_j, pY3_10_x_j, 
                     pY1_01_x_j, pY2_01_x_j, pY3_01_x_j) {
    
    mu01_g_x_ij <- pY2_10_x_i*pY1_01_x_j + pY3_10_x_i*(pY1_01_x_j+pY2_01_x_j)
    mu01_g_x_ji <- pY2_10_x_j*pY1_01_x_i + pY3_10_x_j*(pY1_01_x_i+pY2_01_x_i)
    mu01_l_x_ij <- pY1_10_x_i*(pY2_01_x_j+pY3_01_x_j) + pY2_10_x_i*pY3_01_x_j
    mu01_l_x_ji <- pY1_10_x_j*(pY2_01_x_i+pY3_01_x_i) + pY2_10_x_j*pY3_01_x_i
    
    g01_oioj_g <- 0.5*(z_i*(1-d_i)*e01_x_i/(tp_x_i*ps1_x_i) *
                         (1-z_j)*d_j*e01_x_j/((1-tp_x_j)*(1-ps0_x_j)) * 
                         (1*(y_i>y_j)-mu01_g_x_ij) + 
                         z_j*(1-d_j)*e01_x_j/(tp_x_j*ps1_x_j) * 
                         (1-z_i)*d_i*e01_x_i/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j>y_i)-mu01_g_x_ji) +
                         eta_x_i*((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                                    ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         eta_x_j*((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                                    ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         mu01_g_x_ij +
                         eta_x_j*((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                                    ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         eta_x_i*((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                                    ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         mu01_g_x_ji)
    
    g01_oioj_l <- 0.5*(z_i*(1-d_i)*e01_x_i/(tp_x_i*ps1_x_i) *
                         (1-z_j)*d_j*e01_x_j/((1-tp_x_j)*(1-ps0_x_j)) * 
                         (1*(y_i<y_j)-mu01_l_x_ij) + 
                         z_j*(1-d_j)*e01_x_j/(tp_x_j*ps1_x_j) * 
                         (1-z_i)*d_i*e01_x_i/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j<y_i)-mu01_l_x_ji) +
                         eta_x_i*((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                                    ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         eta_x_j*((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                                    ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         mu01_l_x_ij +
                         eta_x_j*((z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j))-
                                    ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         eta_x_i*((z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i))-
                                    ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         mu01_l_x_ji)
    
    return(c(mean(g01_oioj_g), mean(g01_oioj_l)))
  }
  
  tauN00 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, e00_x_i, tp_x_i, eta_x_i,
                     pY1_10_x_i, pY2_10_x_i, pY3_10_x_i, 
                     pY1_00_x_i, pY2_00_x_i, pY3_00_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, e00_x_j, tp_x_j, eta_x_j,
                     pY1_10_x_j, pY2_10_x_j, pY3_10_x_j, 
                     pY1_00_x_j, pY2_00_x_j, pY3_00_x_j) {
    
    mu00_g_x_ij <- pY2_10_x_i*pY1_00_x_j + pY3_10_x_i*(pY1_00_x_j+pY2_00_x_j)
    mu00_g_x_ji <- pY2_10_x_j*pY1_00_x_i + pY3_10_x_j*(pY1_00_x_i+pY2_00_x_i)
    mu00_l_x_ij <- pY1_10_x_i*(pY2_00_x_j+pY3_00_x_j) + pY2_10_x_i*pY3_00_x_j
    mu00_l_x_ji <- pY1_10_x_j*(pY2_00_x_i+pY3_00_x_i) + pY2_10_x_j*pY3_00_x_i
    
    g00_oioj_g <- 0.5*(z_i*(1-d_i)*e00_x_i/(tp_x_i*(1-ps0_x_i)) * 
                         (1-z_j)*(1-d_j)*e00_x_j/((1-tp_x_j)*(1-ps0_x_j)) *
                         (1*(y_i>y_j)-mu00_g_x_ij) + 
                         z_j*(1-d_j)*e00_x_j/(tp_x_j*(1-ps0_x_j)) * 
                         (1-z_i)*(1-d_i)*e00_x_i/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j>y_i)-mu00_g_x_ji) +
                         (((1-z_i)*(ps0_x_i-d_i)/(1-tp_x_i)+1-ps0_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) * 
                         (((1-z_j)*(ps0_x_j-d_j)/(1-tp_x_j)+1-ps0_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         mu00_g_x_ij + 
                         (((1-z_j)*(ps0_x_j-d_j)/(1-tp_x_j)+1-ps0_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) * 
                         (((1-z_i)*(ps0_x_i-d_i)/(1-tp_x_i)+1-ps0_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         mu00_g_x_ji)
    
    g00_oioj_l <- 0.5*(z_i*(1-d_i)*e00_x_i/(tp_x_i*(1-ps0_x_i)) * 
                         (1-z_j)*(1-d_j)*e00_x_j/((1-tp_x_j)*(1-ps0_x_j)) *
                         (1*(y_i<y_j)-mu00_l_x_ij) + 
                         z_j*(1-d_j)*e00_x_j/(tp_x_j*(1-ps0_x_j)) * 
                         (1-z_i)*(1-d_i)*e00_x_i/((1-tp_x_i)*(1-ps0_x_i)) *
                         (1*(y_j<y_i)-mu00_l_x_ji) +
                         (((1-z_i)*(ps0_x_i-d_i)/(1-tp_x_i)+1-ps0_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) * 
                         (((1-z_j)*(ps0_x_j-d_j)/(1-tp_x_j)+1-ps0_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         mu00_l_x_ij + 
                         (((1-z_j)*(ps0_x_j-d_j)/(1-tp_x_j)+1-ps0_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) * 
                         (((1-z_i)*(ps0_x_i-d_i)/(1-tp_x_i)+1-ps0_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         mu00_l_x_ji)
    
    return(c(mean(g00_oioj_g), mean(g00_oioj_l)))
  }
  
  tauN11 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, e11_x_i, tp_x_i, eta_x_i,
                     pY1_11_x_i, pY2_11_x_i, pY3_11_x_i, 
                     pY1_01_x_i, pY2_01_x_i, pY3_01_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, e11_x_j, tp_x_j, eta_x_j,
                     pY1_11_x_j, pY2_11_x_j, pY3_11_x_j, 
                     pY1_01_x_j, pY2_01_x_j, pY3_01_x_j) {
    
    mu11_g_x_ij <- pY2_11_x_i*pY1_01_x_j + pY3_11_x_i*(pY1_01_x_j+pY2_01_x_j)
    mu11_g_x_ji <- pY2_11_x_j*pY1_01_x_i + pY3_11_x_j*(pY1_01_x_i+pY2_01_x_i)
    mu11_l_x_ij <- pY1_11_x_i*(pY2_01_x_j+pY3_01_x_j) + pY2_11_x_i*pY3_01_x_j
    mu11_l_x_ji <- pY1_11_x_j*(pY2_01_x_i+pY3_01_x_i) + pY2_11_x_j*pY3_01_x_i
    
    g11_oioj_g <- 0.5*(z_i*d_i*e11_x_i/(tp_x_i*ps1_x_i) * (1-z_j)*d_j*e11_x_j/((1-tp_x_j)*ps1_x_j) *
                         (1*(y_i>y_j)-mu11_g_x_ij) +
                         z_j*d_j*e11_x_j/(tp_x_j*ps1_x_j) * (1-z_i)*d_i*e11_x_i/((1-tp_x_i)*ps1_x_i) *
                         (1*(y_j>y_i)-mu11_g_x_ji) +
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) * 
                         mu11_g_x_ij +
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) * 
                         mu11_g_x_ji)
    
    g11_oioj_l <- 0.5*(z_i*d_i*e11_x_i/(tp_x_i*ps1_x_i) * (1-z_j)*d_j*e11_x_j/((1-tp_x_j)*ps1_x_j) *
                         (1*(y_i<y_j)-mu11_l_x_ij) +
                         z_j*d_j*e11_x_j/(tp_x_j*ps1_x_j) * (1-z_i)*d_i*e11_x_i/((1-tp_x_i)*ps1_x_i) *
                         (1*(y_j<y_i)-mu11_l_x_ji) +
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) *
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) * 
                         mu11_l_x_ij +
                         ((z_j*(d_j-ps1_x_j)/tp_x_j+ps1_x_j) - 
                            (z_j*(d_j-ps1_x_j)/(tp_x_j*(1-eta_x_j))+ps1_x_j/(1-eta_x_j)) + 
                            ((1-z_j)*(d_j-ps0_x_j)/((1-tp_x_j)*(1-eta_x_j))+ps0_x_j/(1-eta_x_j))) *
                         ((z_i*(d_i-ps1_x_i)/tp_x_i+ps1_x_i) - 
                            (z_i*(d_i-ps1_x_i)/(tp_x_i*(1-eta_x_i))+ps1_x_i/(1-eta_x_i)) + 
                            ((1-z_i)*(d_i-ps0_x_i)/((1-tp_x_i)*(1-eta_x_i))+ps0_x_i/(1-eta_x_i))) * 
                         mu11_l_x_ji)
    
    return(c(mean(g11_oioj_g), mean(g11_oioj_l)))
  }
  
  
  tauN10 <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           e10_x_i = dfp$e10_x, tp_x_i = dfp$tp_x, eta_x_i = dfp$eta_x, 
           pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           e10_x_j = dfp$e10_x_2, tp_x_j = dfp$tp_x_2, eta_x_j = dfp$eta_x_2, 
           pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN01 <- 
    tauN01(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           e01_x_i = dfp$e01_x, tp_x_i = dfp$tp_x, eta_x_i = dfp$eta_x, 
           pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           e01_x_j = dfp$e01_x_2, tp_x_j = dfp$tp_x_2, eta_x_j = dfp$eta_x_2, 
           pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  tauN00 <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           e00_x_i = dfp$e00_x, tp_x_i = dfp$tp_x, eta_x_i = dfp$eta_x, 
           pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           e00_x_j = dfp$e00_x_2, tp_x_j = dfp$tp_x_2, eta_x_j = dfp$eta_x_2, 
           pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11 <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           e11_x_i = dfp$e11_x, tp_x_i = dfp$tp_x, eta_x_i = dfp$eta_x, 
           pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           e11_x_j = dfp$e11_x_2, tp_x_j = dfp$tp_x_2, eta_x_j = dfp$eta_x_2, 
           pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  # final estimates
  tau10 <- tauN10[1]/tauN10[2]
  tau01 <- tauN01[1]/tauN01[2]
  tau00 <- tauN00[1]/tauN00[2]
  tau11 <- tauN11[1]/tauN11[2]
  
  results <- data.frame(tau10 = tau10, tau01 = tau01, tau00 = tau00, tau11 = tau11)
  
  return(results)
}

