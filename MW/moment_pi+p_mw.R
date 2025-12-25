
# rm(list = ls())
library(tidyverse)

moment_pi_p_mw_est <- function(df) {
  # calculate the estimate using the triply robust estimator
  ## nuisance functions are estimated using parametric models

  m_tp <- glm(z ~ x1 + x2 + x3 + x4, data = df, family = binomial)
  m_tp_m <- glm(z ~ x1t + x2t + x3t + x4t, data = df, family = binomial)
  
  m_ps <- glm(d ~ z + x1 + x2 + x3 + x4, data = df, family = binomial)
  m_ps_m <- glm(d ~ z + x1t + x2t + x3t + x4t, data = df, family = binomial)
  
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
  df$tp_x <- tp_x
  df$tp_x_m <- tp_x_m
  df$ps1_x <- ps1_x
  df$ps0_x <- ps0_x
  df$ps1_x_m <- ps1_x_m
  df$ps0_x_m <- ps0_x_m
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ## different misspecification scenarios
  ### define functions
  tauN10 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j) {
    
    g10_oioj <- 0.5*(z_i*d_i*(ps1_x_i-ps0_x_i)/(tp_x_i*ps1_x_i) *
                       (1-z_j)*(1-d_j)*(ps1_x_j-ps0_x_j)/((1-tp_x_j)*(1-ps0_x_j)) * (1*(y_i>y_j)) + 
                       z_j*d_j*(ps1_x_j-ps0_x_j)/(tp_x_j*ps1_x_j) * 
                       (1-z_i)*(1-d_i)*(ps1_x_i-ps0_x_i)/((1-tp_x_i)*(1-ps0_x_i)) * (1*(y_j>y_i)))
    
    return(mean(g10_oioj))
  }
  
  tauN00 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j) {
    
    g00_oioj <- 0.5*(z_i*(1-d_i)/tp_x_i * (1-z_j)*(1-d_j)*(1-ps1_x_j)/((1-tp_x_j)*(1-ps0_x_j)) * 
                       (1*(y_i>y_j)) + 
                       z_j*(1-d_j)/tp_x_j * (1-z_i)*(1-d_i)*(1-ps1_x_i)/((1-tp_x_i)*(1-ps0_x_i)) *
                       (1*(y_j>y_i)))
    
    return(mean(g00_oioj))
  }
  
  tauN11 <- function(z_i, d_i, y_i, ps1_x_i, ps0_x_i, tp_x_i,
                     z_j, d_j, y_j, ps1_x_j, ps0_x_j, tp_x_j) {
    
    g11_oioj <- 0.5*(z_i*d_i*ps0_x_i/(tp_x_i*ps1_x_i) * (1-z_j)*d_j/(1-tp_x_j) * (1*(y_i>y_j)) + 
                       z_j*d_j*ps0_x_j/(tp_x_j*ps1_x_j) * (1-z_i)*d_i/(1-tp_x_i) * (1*(y_j>y_i)))
    
    return(mean(g11_oioj))
  }
  
  
  ### tp-T, ps-T, oc-T
  tauN10_tpT_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN00_tpT_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN11_tpT_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2)
  
  ### tp-F, ps-T, oc-T
  tauN10_tpF_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN00_tpF_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN11_tpF_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, 
           ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  ### tp-T, ps-F, oc-T
  tauN10_tpT_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN00_tpT_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN11_tpT_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2)
  
  ### tp-T, ps-T, oc-F
  tauN10_tpT_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN00_tpT_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN11_tpT_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_2)
  
  ### tp-F, ps-F, oc-T
  tauN10_tpF_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN00_tpF_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN11_tpF_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  ### tp-F, ps-T, oc-F
  tauN10_tpF_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN00_tpF_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN11_tpF_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x, ps0_x_i = dfp$ps0_x, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_2, ps0_x_j = dfp$ps0_x_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  ### tp-T, ps-F, oc-F
  tauN10_tpT_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN00_tpT_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2)
  
  tauN11_tpT_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_2)
  
  ### tp-F, ps-F, oc-F
  tauN10_tpF_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN00_tpF_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2)
  
  tauN11_tpF_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, y_i = dfp$y, ps1_x_i = dfp$ps1_x_m, ps0_x_i = dfp$ps0_x_m, 
           tp_x_i = dfp$tp_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, y_j = dfp$y_2, ps1_x_j = dfp$ps1_x_m_2, ps0_x_j = dfp$ps0_x_m_2, 
           tp_x_j = dfp$tp_x_m_2)
  
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

