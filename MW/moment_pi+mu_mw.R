
# rm(list = ls())
library(tidyverse)

moment_pi_mu_mw_est <- function(df) {
  # calculate the estimate using the triply robust estimator
  ## nuisance functions are estimated using parametric models

  m_tp <- glm(z ~ x1 + x2 + x3 + x4, data = df, family = binomial)
  m_tp_m <- glm(z ~ x1t + x2t + x3t + x4t, data = df, family = binomial)
  
  m_ps <- glm(d ~ z + x1 + x2 + x3 + x4, data = df, family = binomial)
  m_ps_m <- glm(d ~ z + x1t + x2t + x3t + x4t, data = df, family = binomial)
  
  m_oc <- lm(y ~ d + z + x1 + x2 + x3 + x4, data = df)
  m_oc_m <- lm(y ~ d + z + x1t + x2t + x3t + x4t, data = df)
  
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
  mu11_x <- predict.lm(m_oc, newdata = data.frame(z = 1, d = 1, x1 = df$x1, x2 = df$x2,
                                                  x3 = df$x3, x4 = df$x4))
  mu11_x_m <- predict.lm(m_oc_m, newdata = data.frame(z = 1, d = 1, x1t = df$x1t, x2t = df$x2t,
                                                      x3t = df$x3t, x4t = df$x4t))
  
  mu00_x <- predict.lm(m_oc, newdata = data.frame(z = 0, d = 0, x1 = df$x1, x2 = df$x2,
                                                  x3 = df$x3, x4 = df$x4))
  mu00_x_m <- predict.lm(m_oc_m, newdata = data.frame(z = 0, d = 0, x1t = df$x1t, x2t = df$x2t,
                                                      x3t = df$x3t, x4t = df$x4t))
  
  mu10_x <- predict.lm(m_oc, newdata = data.frame(z = 1, d = 0, x1 = df$x1, x2 = df$x2,
                                                  x3 = df$x3, x4 = df$x4))
  mu10_x_m <- predict.lm(m_oc_m, newdata = data.frame(z = 1, d = 0, x1t = df$x1t, x2t = df$x2t,
                                                      x3t = df$x3t, x4t = df$x4t))
  
  mu01_x <- predict.lm(m_oc, newdata = data.frame(z = 0, d = 1, x1 = df$x1, x2 = df$x2,
                                                  x3 = df$x3, x4 = df$x4))
  mu01_x_m <- predict.lm(m_oc_m, newdata = data.frame(z = 0, d = 1, x1t = df$x1t, x2t = df$x2t,
                                                      x3t = df$x3t, x4t = df$x4t))
  
  sigma <- summary(m_oc)$sigma
  sigma_m <- summary(m_oc_m)$sigma
  
  df$tp_x <- tp_x
  df$tp_x_m <- tp_x_m
  df$ps1_x <- ps1_x
  df$ps0_x <- ps0_x
  df$ps1_x_m <- ps1_x_m
  df$ps0_x_m <- ps0_x_m
  df$mu11_x <- mu11_x
  df$mu11_x_m <- mu11_x_m
  df$mu00_x <- mu00_x
  df$mu00_x_m <- mu00_x_m
  df$mu10_x <- mu10_x
  df$mu10_x_m <- mu10_x_m
  df$mu01_x <- mu01_x
  df$mu01_x_m <- mu01_x_m
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ## different misspecification scenarios
  ### define functions
  tauN10 <- function(z_i, d_i, tp_x_i, mu11_x_i, mu00_x_i, 
                     z_j, d_j, tp_x_j, mu11_x_j, mu00_x_j, sig) {
    
    g10_oioj <- 0.5*((z_i*d_i/tp_x_i-(1-z_i)*d_i/(1-tp_x_i)) *
                       (z_j*d_j/tp_x_j-(1-z_j)*d_j/(1-tp_x_j)) * 
                       pnorm((mu11_x_i-mu00_x_j)/(sqrt(2)*sig)) + 
                       (z_j*d_j/tp_x_j-(1-z_j)*d_j/(1-tp_x_j)) * 
                       (z_i*d_i/tp_x_i-(1-z_i)*d_i/(1-tp_x_i)) *
                       pnorm((mu11_x_j-mu00_x_i)/(sqrt(2)*sig)))

    return(mean(g10_oioj))
  }
  
  tauN00 <- function(z_i, d_i, tp_x_i, mu10_x_i, mu00_x_i,
                     z_j, d_j, tp_x_j, mu10_x_j, mu00_x_j, sig) {
    
    g00_oioj <- 0.5*((1-d_i*z_i/tp_x_i)*(1-d_j*z_j/tp_x_j) * 
                       pnorm((mu10_x_i-mu00_x_j)/(sqrt(2)*sig)) + 
                       (1-d_j*z_j/tp_x_j)*(1-d_i*z_i/tp_x_i) *
                       pnorm((mu10_x_j-mu00_x_i)/(sqrt(2)*sig)))
    
    return(mean(g00_oioj))
  }
  
  tauN11 <- function(z_i, d_i, tp_x_i, mu11_x_i, mu01_x_i,
                     z_j, d_j, tp_x_j, mu11_x_j, mu01_x_j, sig) {
    
    g11_oioj <- 0.5*(d_i*(1-z_i)/(1-tp_x_i) * d_j*(1-z_j)/(1-tp_x_j) *
                       pnorm((mu11_x_i-mu01_x_j)/(sqrt(2)*sig)) + 
                       d_j*(1-z_j)/(1-tp_x_j) * d_i*(1-z_i)/(1-tp_x_i) *
                       pnorm((mu11_x_j-mu01_x_i)/(sqrt(2)*sig)))

    return(mean(g11_oioj))
  }
  
  
  ### tp-T, ps-T, oc-T
  tauN10_tpT_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpT_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpT_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-F, ps-T, oc-T
  tauN10_tpF_psT_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpF_psT_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpF_psT_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-T, ps-F, oc-T
  tauN10_tpT_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpT_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpT_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-T, ps-T, oc-F
  tauN10_tpT_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpT_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpT_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  ### tp-F, ps-F, oc-T
  tauN10_tpF_psF_ocT <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN00_tpF_psF_ocT <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x, mu00_x_i = dfp$mu00_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_2, mu00_x_j = dfp$mu00_x_2, sig = sigma)
  
  tauN11_tpF_psF_ocT <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x, mu01_x_i = dfp$mu01_x,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_2, mu01_x_j = dfp$mu01_x_2, sig = sigma)
  
  ### tp-F, ps-T, oc-F
  tauN10_tpF_psT_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpF_psT_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpF_psT_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  ### tp-T, ps-F, oc-F
  tauN10_tpT_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpT_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpT_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_2, mu11_x_j = dfp$mu11_x_m_2, mu01_x_j = dfp$mu01_x_m_2, sig = sigma_m)
  
  ### tp-F, ps-F, oc-F
  tauN10_tpF_psF_ocF <- 
    tauN10(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu11_x_j = dfp$mu11_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN00_tpF_psF_ocF <- 
    tauN00(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu10_x_i = dfp$mu10_x_m, mu00_x_i = dfp$mu00_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
           tp_x_j = dfp$tp_x_m_2, mu10_x_j = dfp$mu10_x_m_2, mu00_x_j = dfp$mu00_x_m_2, sig = sigma_m)
  
  tauN11_tpF_psF_ocF <- 
    tauN11(z_i = dfp$z, d_i = dfp$d, 
           tp_x_i = dfp$tp_x_m, mu11_x_i = dfp$mu11_x_m, mu01_x_i = dfp$mu01_x_m,
           z_j = dfp$z_2, d_j = dfp$d_2, 
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

