
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
  ps <- plogis(-1 + 2*z + 1*x1 - 0.8*x2 + 0.6*x3 - 1*x4)
  d <- rbinom(n, 1, ps)
  
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
  
  # create 4 mutated covariates to simulate misspecification
  x1t <- exp(0.5*x1)
  x2t <- x2/(1+x1)
  x3t <- (x2*x3/25+0.6)**3
  x4t <- x3*x4/sqrt(2)
  
  # output data
  data_out <- data.frame(id)
  data_out$x1 <- x1
  data_out$x2 <- x2
  data_out$x3 <- x3
  data_out$x4 <- x4
  data_out$x1t <- x1t
  data_out$x2t <- x2t
  data_out$x3t <- x3t
  data_out$x4t <- x4t
  data_out$z <- z
  data_out$d <- d
  data_out$y <- y
  
  return(data_out)
}

# data_gen_mw(100)
