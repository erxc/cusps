
# rm(list = ls())
library(tidyverse)
library(MASS)

triply_robust_wr_est <- function(df) {
  # calculate the estimate using the triply robust estimator
  ## nuisance functions are estimated using parametric models

  m_tp <- glm(z ~ x1 + x2 + x3 + x4, data = df, family = binomial)
  
  m_ps <- glm(d ~ z + x1 + x2 + x3 + x4, data = df, family = binomial)
  
  m_oc <- polr(factor(y) ~ d + z + x1 + x2 + x3 + x4, data = df, Hess = T)
  
  # calculate the denominators
  tp_x <- predict.glm(m_tp, newdata = data.frame(df[,c('x1','x2','x3','x4')]), type = "response")
  
  ps1_x <- 
    predict.glm(m_ps, newdata = data.frame(z = 1, x1 = df$x1, x2 = df$x2,
                                           x3 = df$x3, x4 = df$x4), type = "response")
  
  ps0_x <- 
    predict.glm(m_ps, newdata = data.frame(z = 0, x1 = df$x1, x2 = df$x2,
                                           x3 = df$x3, x4 = df$x4), type = "response")
  
  # sensitivity function
  eta_x <- 0.99*(1 - (ps1_x - ps0_x)/pmin(ps1_x, 1 - ps0_x))
  
  # principal scores combined with sensitivity function
  e10_x <- (ps1_x - ps0_x)/(1 - eta_x)
  e01_x <- eta_x*(ps1_x - ps0_x)/(1 - eta_x)
  e00_x <- 1 - ps0_x - (ps1_x - ps0_x)/(1 - eta_x)
  e11_x <- ps1_x - (ps1_x - ps0_x)/(1 - eta_x)
  
  p1 <- mean(df$z*(df$d-ps1_x)/tp_x) + mean(ps1_x)
  p0 <- mean((1-df$z)*(df$d-ps0_x)/(1-tp_x)) + mean(ps0_x)
  
  e10 <- mean(e10_x)
  e01 <- mean(e01_x)
  e00 <- mean(e00_x)
  e11 <- mean(e11_x)
  
  den10 <- e10**2
  den01 <- e01**2
  den00 <- e00**2
  den11 <- e11**2
  
  # calculate the numerators
  pmat11_x <- predict(m_oc, newdata = data.frame(z = 1, d = 1, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  
  pmat00_x <- predict(m_oc, newdata = data.frame(z = 0, d = 0, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  
  pmat10_x <- predict(m_oc, newdata = data.frame(z = 1, d = 0, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  
  pmat01_x <- predict(m_oc, newdata = data.frame(z = 0, d = 1, x1 = df$x1, x2 = df$x2,
                                                 x3 = df$x3, x4 = df$x4), type = 'p')
  
  df$tp_x <- tp_x
  df$eta_x <- eta_x
  df$ps1_x <- ps1_x
  df$ps0_x <- ps0_x
  df$e10_x <- e10_x
  df$e01_x <- e01_x
  df$e00_x <- e00_x
  df$e11_x <- e11_x
  df$pY1_11_x <- pmat11_x[,1]
  df$pY2_11_x <- pmat11_x[,2]
  df$pY3_11_x <- pmat11_x[,3]
  df$pY1_00_x <- pmat00_x[,1]
  df$pY2_00_x <- pmat00_x[,2]
  df$pY3_00_x <- pmat00_x[,3]
  df$pY1_10_x <- pmat10_x[,1]
  df$pY2_10_x <- pmat10_x[,2]
  df$pY3_10_x <- pmat10_x[,3]
  df$pY1_01_x <- pmat01_x[,1]
  df$pY2_01_x <- pmat01_x[,2]
  df$pY3_01_x <- pmat01_x[,3]
  
  ## create a dataframe for paired observations
  ### this creates a dataset with n*(n-1) pairs
  dfp <- df %>% 
    setNames(paste0(names(.), '_2')) %>% 
    crossing(df) %>% 
    filter(id != id_2)
  
  ## different misspecification scenarios
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

