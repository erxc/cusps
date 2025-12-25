
# rm(list = ls())

data_gen_wr <- function(n) {
  # id
  id <- 1:n
  
  # 4 covariates: 3 normal, 1 bernoulli
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rbinom(n, 1, 0.5)
  
  # treatment assignment/probability
  tp <- plogis(-1*x1 + 0.5*x2 - 0.25*x3 - 0.1*x4)
  z <- rbinom(n, 1, tp)
  
  # principal scores/intermediate variables
  p1 <- plogis(-1 + 2*1 + 1*x1 - 0.8*x2 + 0.6*x3 - 1*x4)
  p0 <- plogis(-1 + 2*0 + 1*x1 - 0.8*x2 + 0.6*x3 - 1*x4)
  
  eta <- 0.99*(1 - (p1 - p0)/pmin(p1, 1 - p0))
  
  e10 <- (p1 - p0)/(1 - eta)
  e01 <- eta*(p1 - p0)/(1 - eta)
  e00 <- 1 - p0 - (p1 - p0)/(1 - eta)
  e11 <- p1 - (p1 - p0)/(1 - eta)
  
  s <- rep(NA, n)
  for (i in 1:n) {
    s[i] <- sample(1:4, 1, prob = c(e10[i], e01[i], e00[i], e11[i]))
  }  
  
  d <- rep(NA, n)
  for (i in 1:n) {
    if (z[i] == 1 & (s[i] == 1 | s[i] == 4)) {d[i] <- 1}
    if (z[i] == 0 & (s[i] == 2 | s[i] == 4)) {d[i] <- 1}
    if (z[i] == 1 & (s[i] == 2 | s[i] == 3)) {d[i] <- 0}
    if (z[i] == 0 & (s[i] == 1 | s[i] == 3)) {d[i] <- 0}
  }
  
  # outcome: ordinal normal
  sigma <- 1 # standard deviation
  p_cat1 <- plogis(-1 + 2*z - 1*d + 1*x1 - 1*x2 + 1.2*x3 - 0.8*x4)
  p_cat2 <- plogis(1 + 2*z - 1*d + 1*x1 - 1*x2 + 1.2*x3 - 0.8*x4) - 
    plogis(-1 + 2*z - 1*d + 1*x1 - 1*x2 + 1.2*x3 - 0.8*x4)
  p_cat3 <- 1 - (p_cat1 + p_cat2)
  
  y <- rep(NA, n)
  for (i in 1:n) {
    y[i] <- sample(1:3, 1, prob = c(p_cat1[i], p_cat2[i], p_cat3[i]))
  }
  
  
  # output data
  data_out <- data.frame(id)
  data_out$x1 <- x1
  data_out$x2 <- x2
  data_out$x3 <- x3
  data_out$x4 <- x4
  data_out$z <- z
  data_out$s <- s
  data_out$d <- d
  data_out$y <- y
  
  return(data_out)
}



