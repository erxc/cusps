
# rm(list = ls())
library(MASS)
library(tidyverse)

moment_p_mu_wr_est <- function(df) {
  # calculate the estimate using the triply robust estimator
  ## nuisance functions are estimated using parametric models

  m_tp <- glm(z ~ x1 + x2 + x3 + x4, data = df, family = binomial)
  m_tp_m <- glm(z ~ x1t + x2t + x3t + x4t, data = df, family = binomial)
  
  m_ps <- glm(d ~ z + x1 + x2 + x3 + x4, data = df, family = binomial)
  m_ps_m <- glm(d ~ z + x1t + x2t + x3t + x4t, data = df, family = binomial)
  
  m_oc <- polr(factor(y) ~ d + z + x1 + x2 + x3 + x4, data = df, Hess = T)
  m_oc_m <- polr(factor(y) ~ d + z + x1t + x2t + x3t + x4t, data = df, Hess = T)
  
  # calculate the denominators
  tp_x <- predict.glm(m_tp, newdata = data.frame(df[,c('x1','x2','x3','x4')]), type = "response")
  
  tp_x_m <- predict.glm(m_tp_m, newdata = data.frame(df[,c('x1t','x2t','x3t','x4t')]), type = "response")
  
  ps1_x <- 
    predict.glm(m_ps, newdata = data.frame(z = 1, x1 = df$x1, x2 = df$x2,
                                           x3 = df$x3, x4 = df$x4), type = "response")
  
  ps0_x <- 
    predict.glm(m_ps, newdata = data.frame(z = 0, x1 = df$x1, x2 = df$x2,
                                           x3 = df$x3, x4 = df$x4), type = "response")
  
  ps1_x_m <- 
    predict.glm(m_ps_m, newdata = data.frame(z = 1, x1t = df$x1t, x2t = df$x2t,
                                             x3t = df$x3t, x4t = df$x4t), type = "response")
  
  ps0_x_m <- 
    predict.glm(m_ps_m, newdata = data.frame(z = 0, x1t = df$x1t, x2t = df$x2t,
                                             x3t = df$x3t, x4t = df$x4t), type = "response")
  
  ## different misspecification scenarios
  ### tp-T, ps-T
  p1_tpT_psT <- mean(df$z*(df$d-ps1_x)/tp_x) + mean(ps1_x)
  p0_tpT_psT <- mean((1-df$z)*(df$d-ps0_x)/(1-tp_x)) + mean(ps0_x)
  
  den10_tpT_psT <- (p1_tpT_psT - p0_tpT_psT)**2
  den00_tpT_psT <- (1 - p1_tpT_psT)**2
  den11_tpT_psT <- p0_tpT_psT**2
  
  ### tp-T, ps-F
  p1_tpT_psF <- mean(df$z*(df$d-ps1_x_m)/tp_x) + mean(ps1_x_m)
  p0_tpT_psF <- mean((1-df$z)*(df$d-ps0_x_m)/(1-tp_x)) + mean(ps0_x_m)
  
  den10_tpT_psF <- (p1_tpT_psF - p0_tpT_psF)**2
  den00_tpT_psF <- (1 - p1_tpT_psF)**2
  den11_tpT_psF <- p0_tpT_psF**2
  
  ### tp-F, ps-T
  p1_tpF_psT <- mean(df$z*(df$d-ps1_x)/tp_x_m) + mean(ps1_x)
  p0_tpF_psT <- mean((1-df$z)*(df$d-ps0_x)/(1-tp_x_m)) + mean(ps0_x)
  
  den10_tpF_psT <- (p1_tpF_psT - p0_tpF_psT)**2
  den00_tpF_psT <- (1 - p1_tpF_psT)**2
  den11_tpF_psT <- p0_tpF_psT**2
  
  ### tp-F, ps-F
  p1_tpF_psF <- mean(df$z*(df$d-ps1_x_m)/tp_x_m) + mean(ps1_x_m)
  p0_tpF_psF <- mean((1-df$z)*(df$d-ps0_x_m)/(1-tp_x_m)) + mean(ps0_x_m)
  
  den10_tpF_psF <- (p1_tpF_psF - p0_tpF_psF)**2
  den00_tpF_psF <- (1 - p1_tpF_psF)**2
  den11_tpF_psF <- p0_tpF_psF**2
  
  # calculate the numerators
  pmat11_x <- predict(m_oc, newdata = data.frame(z = 1, d = 1, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  pmat11_x_m <- predict(m_oc_m, newdata = data.frame(z = 1, d = 1, x1t = df$x1t, x2t = df$x2t,
                                                     x3t = df$x3t, x4t = df$x4t), type = 'p')
  
  pmat00_x <- predict(m_oc, newdata = data.frame(z = 0, d = 0, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  pmat00_x_m <- predict(m_oc_m, newdata = data.frame(z = 0, d = 0, x1t = df$x1t, x2t = df$x2t,
                                                     x3t = df$x3t, x4t = df$x4t), type = 'p')
  
  pmat10_x <- predict(m_oc, newdata = data.frame(z = 1, d = 0, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  pmat10_x_m <- predict(m_oc_m, newdata = data.frame(z = 1, d = 0, x1t = df$x1t, x2t = df$x2t,
                                                     x3t = df$x3t, x4t = df$x4t), type = 'p')
  
  pmat01_x <- predict(m_oc, newdata = data.frame(z = 0, d = 1, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  pmat01_x_m <- predict(m_oc_m, newdata = data.frame(z = 0, d = 1, x1t = df$x1t, x2t = df$x2t,
                                                     x3t = df$x3t, x4t = df$x4t), type = 'p')
  
  df$tp_x <- tp_x
  df$tp_x_m <- tp_x_m
  df$ps1_x <- ps1_x
  df$ps0_x <- ps0_x
  df$ps1_x_m <- ps1_x_m
  df$ps0_x_m <- ps0_x_m
  df$pY1_11_x <- pmat11_x[,1]
  df$pY2_11_x <- pmat11_x[,2]
  df$pY3_11_x <- pmat11_x[,3]
  df$pY1_11_x_m <- pmat11_x_m[,1]
  df$pY2_11_x_m <- pmat11_x_m[,2]
  df$pY3_11_x_m <- pmat11_x_m[,3]
  df$pY1_00_x <- pmat00_x[,1]
  df$pY2_00_x <- pmat00_x[,2]
  df$pY3_00_x <- pmat00_x[,3]
  df$pY1_00_x_m <- pmat00_x_m[,1]
  df$pY2_00_x_m <- pmat00_x_m[,2]
  df$pY3_00_x_m <- pmat00_x_m[,3]
  df$pY1_10_x <- pmat10_x[,1]
  df$pY2_10_x <- pmat10_x[,2]
  df$pY3_10_x <- pmat10_x[,3]
  df$pY1_10_x_m <- pmat10_x_m[,1]
  df$pY2_10_x_m <- pmat10_x_m[,2]
  df$pY3_10_x_m <- pmat10_x_m[,3]
  df$pY1_01_x <- pmat01_x[,1]
  df$pY2_01_x <- pmat01_x[,2]
  df$pY3_01_x <- pmat01_x[,3]
  df$pY1_01_x_m <- pmat01_x_m[,1]
  df$pY2_01_x_m <- pmat01_x_m[,2]
  df$pY3_01_x_m <- pmat01_x_m[,3]
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ## different misspecification scenarios
  ### define functions
  tauN10 <- function(ps1_x_i, ps0_x_i, pY1_11_x_i, pY2_11_x_i, pY3_11_x_i, 
                     pY1_00_x_i, pY2_00_x_i, pY3_00_x_i, 
                     ps1_x_j, ps0_x_j, pY1_11_x_j, pY2_11_x_j, pY3_11_x_j, 
                     pY1_00_x_j, pY2_00_x_j, pY3_00_x_j) {
    
    mu10_g_x_ij <- pY2_11_x_i*pY1_00_x_j + pY3_11_x_i*(pY1_00_x_j+pY2_00_x_j)
    mu10_g_x_ji <- pY2_11_x_j*pY1_00_x_i + pY3_11_x_j*(pY1_00_x_i+pY2_00_x_i)
    mu10_l_x_ij <- pY1_11_x_i*(pY2_00_x_j+pY3_00_x_j) + pY2_11_x_i*pY3_00_x_j
    mu10_l_x_ji <- pY1_11_x_j*(pY2_00_x_i+pY3_00_x_i) + pY2_11_x_j*pY3_00_x_i
    
    g10_oioj_g <- 0.5*((ps1_x_i-ps0_x_i) * (ps1_x_j-ps0_x_j) * mu10_g_x_ij + 
                         (ps1_x_j-ps0_x_j) * (ps1_x_i-ps0_x_i) * mu10_g_x_ji)
    
    g10_oioj_l <- 0.5*((ps1_x_i-ps0_x_i) * (ps1_x_j-ps0_x_j) * mu10_l_x_ij + 
                         (ps1_x_j-ps0_x_j) * (ps1_x_i-ps0_x_i) * mu10_l_x_ji)
    
    return(c(mean(g10_oioj_g), mean(g10_oioj_l)))
  }
  
  tauN00 <- function(ps1_x_i, pY1_10_x_i, pY2_10_x_i, pY3_10_x_i, 
                     pY1_00_x_i, pY2_00_x_i, pY3_00_x_i, 
                     ps1_x_j, pY1_10_x_j, pY2_10_x_j, pY3_10_x_j, 
                     pY1_00_x_j, pY2_00_x_j, pY3_00_x_j) {
    
    mu00_g_x_ij <- pY2_10_x_i*pY1_00_x_j + pY3_10_x_i*(pY1_00_x_j+pY2_00_x_j)
    mu00_g_x_ji <- pY2_10_x_j*pY1_00_x_i + pY3_10_x_j*(pY1_00_x_i+pY2_00_x_i)
    mu00_l_x_ij <- pY1_10_x_i*(pY2_00_x_j+pY3_00_x_j) + pY2_10_x_i*pY3_00_x_j
    mu00_l_x_ji <- pY1_10_x_j*(pY2_00_x_i+pY3_00_x_i) + pY2_10_x_j*pY3_00_x_i
    
    g00_oioj_g <- 0.5*((1-ps1_x_i) * (1-ps1_x_j) * mu00_g_x_ij + 
                         (1-ps1_x_j) * (1-ps1_x_i) * mu00_g_x_ji)
    
    g00_oioj_l <- 0.5*((1-ps1_x_i) * (1-ps1_x_j) * mu00_l_x_ij + 
                         (1-ps1_x_j) * (1-ps1_x_i) * mu00_l_x_ji)
    
    return(c(mean(g00_oioj_g), mean(g00_oioj_l)))
  }
  
  tauN11 <- function(ps0_x_i, pY1_11_x_i, pY2_11_x_i, pY3_11_x_i, 
                     pY1_01_x_i, pY2_01_x_i, pY3_01_x_i, 
                     ps0_x_j, pY1_11_x_j, pY2_11_x_j, pY3_11_x_j, 
                     pY1_01_x_j, pY2_01_x_j, pY3_01_x_j) {
    
    mu11_g_x_ij <- pY2_11_x_i*pY1_01_x_j + pY3_11_x_i*(pY1_01_x_j+pY2_01_x_j)
    mu11_g_x_ji <- pY2_11_x_j*pY1_01_x_i + pY3_11_x_j*(pY1_01_x_i+pY2_01_x_i)
    mu11_l_x_ij <- pY1_11_x_i*(pY2_01_x_j+pY3_01_x_j) + pY2_11_x_i*pY3_01_x_j
    mu11_l_x_ji <- pY1_11_x_j*(pY2_01_x_i+pY3_01_x_i) + pY2_11_x_j*pY3_01_x_i
    
    g11_oioj_g <- 0.5*(ps0_x_i * ps0_x_j * mu11_g_x_ij + ps0_x_j * ps0_x_i * mu11_g_x_ji)
    
    g11_oioj_l <- 0.5*(ps0_x_i * ps0_x_j * mu11_l_x_ij + ps0_x_j * ps0_x_i * mu11_l_x_ji)
    
    return(c(mean(g11_oioj_g), mean(g11_oioj_l)))
  }
  
  
  ### tp-T, ps-T, oc-T
  tauN10_tpT_psT_ocT <- 
    tauN10(ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpT_psT_ocT <- 
    tauN00(ps1_x_i = dfp$ps1_x, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpT_psT_ocT <- 
    tauN11(ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-F, ps-T, oc-T
  tauN10_tpF_psT_ocT <- 
    tauN10(ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpF_psT_ocT <- 
    tauN00(ps1_x_i = dfp$ps1_x, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpF_psT_ocT <- 
    tauN11(ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-T, ps-F, oc-T
  tauN10_tpT_psF_ocT <- 
    tauN10(ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpT_psF_ocT <- 
    tauN00(ps1_x_i = dfp$ps1_x_m, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_m_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpT_psF_ocT <- 
    tauN11(ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-T, ps-T, oc-F
  tauN10_tpT_psT_ocF <- 
    tauN10(ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpT_psT_ocF <- 
    tauN00(ps1_x_i = dfp$ps1_x, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpT_psT_ocF <- 
    tauN11(ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  ### tp-F, ps-F, oc-T
  tauN10_tpF_psF_ocT <- 
    tauN10(ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN00_tpF_psF_ocT <- 
    tauN00(ps1_x_i = dfp$ps1_x_m, pY1_10_x_i = dfp$pY1_10_x, pY2_10_x_i = dfp$pY2_10_x, 
           pY3_10_x_i = dfp$pY3_10_x, pY1_00_x_i = dfp$pY1_00_x, pY2_00_x_i = dfp$pY2_00_x, 
           pY3_00_x_i = dfp$pY3_00_x,
           ps1_x_j = dfp$ps1_x_m_2, pY1_10_x_j = dfp$pY1_10_x_2, pY2_10_x_j = dfp$pY2_10_x_2, 
           pY3_10_x_j = dfp$pY3_10_x_2, pY1_00_x_j = dfp$pY1_00_x_2, pY2_00_x_j = dfp$pY2_00_x_2, 
           pY3_00_x_j = dfp$pY3_00_x_2)
  
  tauN11_tpF_psF_ocT <- 
    tauN11(ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x, pY2_11_x_i = dfp$pY2_11_x, 
           pY3_11_x_i = dfp$pY3_11_x, pY1_01_x_i = dfp$pY1_01_x, pY2_01_x_i = dfp$pY2_01_x, 
           pY3_01_x_i = dfp$pY3_01_x,
           ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_2, pY2_11_x_j = dfp$pY2_11_x_2, 
           pY3_11_x_j = dfp$pY3_11_x_2, pY1_01_x_j = dfp$pY1_01_x_2, pY2_01_x_j = dfp$pY2_01_x_2, 
           pY3_01_x_j = dfp$pY3_01_x_2)
  
  ### tp-F, ps-T, oc-F
  tauN10_tpF_psT_ocF <- 
    tauN10(ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpF_psT_ocF <- 
    tauN00(ps1_x_i = dfp$ps1_x, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpF_psT_ocF <- 
    tauN11(ps0_x_i = dfp$ps0_x, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           ps0_x_j = dfp$ps0_x_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  ### tp-T, ps-F, oc-F
  tauN10_tpT_psF_ocF <- 
    tauN10(ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpT_psF_ocF <- 
    tauN00(ps1_x_i = dfp$ps1_x_m, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_m_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpT_psF_ocF <- 
    tauN11(ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_01_x_j = dfp$pY1_01_x_m_2, pY2_01_x_j = dfp$pY2_01_x_m_2, 
           pY3_01_x_j = dfp$pY3_01_x_m_2)
  
  ### tp-F, ps-F, oc-F
  tauN10_tpF_psF_ocF <- 
    tauN10(ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
           pY3_11_x_j = dfp$pY3_11_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN00_tpF_psF_ocF <- 
    tauN00(ps1_x_i = dfp$ps1_x_m, pY1_10_x_i = dfp$pY1_10_x_m, pY2_10_x_i = dfp$pY2_10_x_m, 
           pY3_10_x_i = dfp$pY3_10_x_m, pY1_00_x_i = dfp$pY1_00_x_m, pY2_00_x_i = dfp$pY2_00_x_m, 
           pY3_00_x_i = dfp$pY3_00_x_m,
           ps1_x_j = dfp$ps1_x_m_2, pY1_10_x_j = dfp$pY1_10_x_m_2, pY2_10_x_j = dfp$pY2_10_x_m_2, 
           pY3_10_x_j = dfp$pY3_10_x_m_2, pY1_00_x_j = dfp$pY1_00_x_m_2, pY2_00_x_j = dfp$pY2_00_x_m_2, 
           pY3_00_x_j = dfp$pY3_00_x_m_2)
  
  tauN11_tpF_psF_ocF <- 
    tauN11(ps0_x_i = dfp$ps0_x_m, pY1_11_x_i = dfp$pY1_11_x_m, pY2_11_x_i = dfp$pY2_11_x_m, 
           pY3_11_x_i = dfp$pY3_11_x_m, pY1_01_x_i = dfp$pY1_01_x_m, pY2_01_x_i = dfp$pY2_01_x_m, 
           pY3_01_x_i = dfp$pY3_01_x_m,
           ps0_x_j = dfp$ps0_x_m_2, pY1_11_x_j = dfp$pY1_11_x_m_2, pY2_11_x_j = dfp$pY2_11_x_m_2, 
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

