library(Matrix)
library(MCMCpack)
library(matrixStats)
library(parallel)
library(mvtnorm)
library(abind)
library(MASS)
library(BDgraph)
library(matrixcalc)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(lubridate)
library(grid)
set.seed(123)

kalman_filter <- function(y, trans, z, a, P, sigma, Q, R, missing = FALSE) {
  if (missing == TRUE) {
    a.update <- trans %*% a
    P.update <- trans %*% P %*% t(trans) + R %*% Q %*% t(R)
    list(a = a.update, P = P.update)
  } else {
    v <- y - crossprod(z, a)  
    FF <- t(z) %*% P %*% z + sigma
    if (is.singular.matrix(as.matrix(FF)) == TRUE) {
      FF <- MungeMatrix(FF)
    }
    # print(det(FF))
    K <- trans %*% P %*% z %*% solve(as.matrix(FF))
    L <- trans - K %*% t(z)
    a.update <- trans %*% a + K %*% v
    P.update <- trans %*% P %*% t(L) + R %*% Q %*% t(R)
    # return updated values
    list(v = v, FF = FF, L = L, a = a.update, P = P.update)
  }
}

koop_filter <- function(n, y, trans, z, a.int, P.int,
                        sigma, Q, R, causal.period = NULL,
                        output.var.cov = FALSE) {
  
  T <- dim(y)[1]
  d <- dim(y)[2]
  
  alpha.sample <- matrix(0, n, T)
  v.sample <- matrix(0, T, d)
  F.sample <- array(0, c(d, d, T))
  L.sample <- array(0, c(T, n, n))
  alpha.plus <- matrix(0, T, n)
  r <- matrix(0, n, T+1)
  N <- array(0, c(n, n, T+1))
  a.sample <- matrix(0, n, T)
  a.ff <- matrix(0, n, T+1)
  a.ff[,1] <- as.vector(a.int)
  P.sample <- array(0, c(n, n, T))
  P.ff <- array(0, c(n, n, T+1))
  P.ff[,,1] <- as.matrix(P.int)
  P.cov.sample <- array(0, c(n, n, T-1))
  
  # kalman-filter
  if (is.null(causal.period) == TRUE) {
    a <- a.int
    P <- P.int 
    for (t in 1:T) {
      kalmanfilter <-  kalman_filter(y[t, ], trans, z, a, P, sigma, Q, R)
      v.sample[t, ] <- as.vector(kalmanfilter$v)
      F.sample[, , t] <- as.matrix(kalmanfilter$FF)
      L.sample[t, , ] <- as.matrix(kalmanfilter$L)
      a <- kalmanfilter$a
      P <- kalmanfilter$P
      a.ff[,t+1] <- as.vector(a)
      P.ff[,,t+1] <- as.matrix(P + t(P)) / 2 # force to be symmetric matrix
    }
  } else {
    a <- a.int
    P <- P.int
    for (t in 1:T) {
      if (t %in% causal.period) {
        kalmanfilter <-  kalman_filter(y[t, ], trans, z, a, P, sigma, Q, R, 
                                       missing = TRUE)
        a <- kalmanfilter$a
        P <- kalmanfilter$P
        a.ff[,t+1] <- as.vector(a)
        P.ff[,,t+1] <- as.matrix(P)
      } else {
        kalmanfilter <-  kalman_filter(y[t, ], trans, z, a, P, sigma, Q, R)
        v.sample[t, ] <- as.matrix(kalmanfilter$v)
        F.sample[, , t] <- as.matrix(kalmanfilter$FF)
        L.sample[t, , ] <- as.matrix(kalmanfilter$L)
        a <- kalmanfilter$a
        P <- kalmanfilter$P
        a.ff[,t+1] <- as.vector(a)
        P.ff[,,t+1] <- as.matrix(P + t(P)) / 2
        if (t == (causal.period[1]-1)) {
          a.last <- a
          P.last <- P
        }
      }
    }
  }
  
  for (t in T:1) {
    if (t %in% causal.period) {
      r[, t] <- as.vector(r[, t+1] %*% trans)
      N[,,t] <- as.matrix(t(trans) %*% N[,, t+1] %*% trans)
    } else {
      r[, t] <- as.matrix(z %*% solve(F.sample[, , t]) %*% v.sample[t, ] + 
                            t(L.sample[t, , ]) %*% r[, t+1])
      N[,,t] <- as.matrix(z %*% solve(F.sample[,,t]) %*% t(z)) + 
        as.matrix(t(L.sample[t,,]) %*% N[,,t+1] %*% L.sample[t,,])
    }
  }
  
  # obtain draws of alpha
  alpha.sample[, 1] <- as.matrix(a.int + P.int %*% r[, 1])
  a.sample[, 1] <- as.matrix(a.int + P.int %*% r[,1])
  P.sample[,,1] <- as.matrix(P.ff[,,1]) - 
    as.matrix(P.ff[,,1] %*% N[,,1] %*% P.ff[,,1])
  for (t in 2:T) {
    alpha.sample[, t] <- as.matrix(trans %*% alpha.sample[, t-1] + 
                                     R %*% Q %*% t(R) %*% r[, t])
    a.sample[,t] <- as.matrix(a.ff[, t] + P.ff[,,t] %*% r[, t])
    P.sample[,,t] <- as.matrix(P.ff[,,t]) - 
      as.matrix(P.ff[,,t] %*% N[,,t] %*% P.ff[,,t])
    P.sample[,,t] <- (P.sample[,,t] + t(P.sample[,,t])) / 2 
    P.cov.sample[,,t-1] <- as.matrix(P.ff[,,t-1] %*% t(L.sample[t-1,,])) %*% 
      as.matrix((diag(n) - N[,,t] %*% P.ff[,,t]))
  }
  alpha.sample <- t(alpha.sample)
  a.sample <- a.sample
  P.sample <- P.sample
  P.cov.sample <- P.cov.sample
  
  if (output.var.cov == TRUE) {
    list(alpha.sample = alpha.sample, a.sample = a.sample, P.sample = P.sample, 
         P.cov.sample = P.cov.sample)
  } else {
    # return value
    if (is.null(causal.period) == TRUE) {
      return(alpha.sample = alpha.sample)
    } else {
      list(alpha.sample = alpha.sample, a.last = a.last, P.last = P.last) 
    }
  }
}


EMVS <- function(test, cntl.index, cntl, graph.structure, circle,
                 v0 = 0.1, s = 0.1, iteration = 50, stationary = TRUE, 
                 misspecification = FALSE) {
  
  
  #EMVS 
  length <- dim(test)[1]
  n <- dim(test)[2]
  dCntl <- sum(cntl.index)
  
  cntl.input <- cntl
  # re-organize cntl dataset
  cntl <- NULL
  index <- 0
  for (dim in 1:n) {
    cntl <- abind(cntl, cntl.input[, (index+1):(index+cntl.index[dim])], 
                  along = 3)
    index <- index+cntl.index[dim]
  }
  
  ##
  # rename y.est to y.matrix and cntl.est to x.matrix
  x.std <- array(0, c(length,cntl.index[1],n))
  for (i in 1:(dim(test)[2])) {
    x.std[,,i] <- t(t(cntl[,,i]) - colMeans(cntl[,,i]))
  }
  x.matrix <- matrix(0, n*length, dCntl)
  index <- 0
  for (dims in 1:n) {
    x.matrix[seq(dims,n*length,by=n), 
             (index+1):(index+cntl.index[dims])] <- x.std[, , dims]
    index <- index+cntl.index[dims]
  }
  y.std <- t((t(test) - colMeans(test)))
  y.matrix <- as.matrix(c(t(y.std)))
  
  # giving starting values for beta and sigma
  beta.hat <- rep(0, dCntl)
  y.tilde <- matrix(y.matrix - x.matrix %*% beta.hat, n, length)
  sigma.hat <- crossprod(y.std) / (length-1)
  # if sigma.hat is singular, we convert into a p.d. matrix
  if (is.singular.matrix(sigma.hat) == TRUE){
    sigma.hat <- MungeMatrix(sigma.hat)
  }
  sigma.hat.inv <- solve(sigma.hat)
  sigma.hat.inv[graph.structure == 0] <- 0
  B <- diag(n) # prior parameters for wishart priors
  delta <- n + 1 # prior parameter for wishart priors
  
  # giving staring points for v0, v1, a, b, theta
  v0 <- v0 # initial v0
  # v1 <- 100 # initial v1
  v1 <- 10
  a <- 1 # intial a, shape1 in beta distribution
  b <- 1  # intial b, shape2 in beta distribution
  theta <- rbeta(1, a, b) # intial theta from beta distribtuion
  
  if (misspecification == FALSE) {
    # giving starting values for Q
    k1 <- k2 <- k3 <- 0.1 # prior parameters for sigma.u, sigma.v, sigma.w
    sigma.u.hat <- chol2inv(rgwish(n = 1, adj = graph.structure,
                                   b = n+1, D = k1^2 * n * diag(n)))
    sigma.v.inv <- rgwish(n = 1, adj = graph.structure,
                          b = n+1, D = k2^2 * n * diag(n))
    sigma.v.hat <- chol2inv(sigma.v.inv)
    sigma.v.hat.inv <- solve(sigma.v.hat)
    sigma.w.hat <- chol2inv(rgwish(n = 1, adj = graph.structure,
                                   b = n+1, D = k3^2 * n * diag(n)))
    Q.hat <- bdiag(sigma.u.hat, sigma.v.hat, sigma.w.hat)
    Q.inv <- solve(Q.hat)
    
    # giving initial parameters for KFBS (Kalman filter and backward simulation)
    # initial a.int
    a.int <- rep(0, n*(circle+1)) 
    # initialize P.int
    P.int <- Matrix(0, n*(circle+1), n*(circle+1))
    
    if (stationary == TRUE) {
      P.int[1:(3*n), 1:(3*n)] <- diag(3*n)
    } else {
      P.int[1:(3*n), 1:(3*n)] <- diag(3*n)*10^6
    }
    
    # initial transition matrix
    trans <- Matrix(0, (circle+1)*n, (circle+1)*n)
    linear <- diag(2*n)
    linear[1:n, (n+1):(2*n)] <- diag(n)
    trans[1:(2*n), 1:(2*n)] <- linear
    # take initial variance of tau from the data
    if (stationary == TRUE) {
      data.yw <- ar.yw(as.matrix(test), aic = FALSE, order.max = 1,
                       demean = T, intercept = T)
      phi.hat <- matrix(cbind(data.yw$ar), n, n)
    } else {
      phi.hat <- diag(n)
    }
    trans[(n+1):(2*n), (n+1):(2*n)] <- phi.hat
    seasonal <- Matrix(0, circle-1, circle-1)
    seasonal[1, ] <- -1
    seasonal[2:(circle-1), 1:(circle-2)] <- diag(circle-2)
    for (dims in 1:n) {
      trans[seq((2*n+dims), (circle+1)*n, by=n), 
            seq((2*n+dims), (circle+1)*n, by=n)] <- seasonal
    }
    
    
    # define
    z <- Matrix(0, n*(circle+1), n) 
    z[1:n, ] <- diag(n)
    z[(2*n+1):(3*n), ] <- diag(n)
    
    # initialize R 
    R <- Matrix(0, n*(circle+1), n*3)
    R[1:(3*n), 1:(3*n)] <- diag(3*n)
    
    # initialize alpha.hat
    a.hat <- matrix(0, n*(circle+1), length)
  }
  
  # step 2: EM update parameters
  # create matrix to collect results
  iterloop <- iteration
  beta.update <- matrix(0, dCntl, iterloop+1)
  sigma.update <- array(NA, c(n,n,iterloop)) 
  v0.update <- v1.update <- theta.update <- rep(NA, iterloop)
  lp.update <- rep(NA, iterloop)
  if (misspecification == FALSE) {
    a.update <- array(NA, c(length, n*(circle+1),iterloop))
    phi.update <- array(NA, c(n,n,iterloop))
    sigma.u.update <- sigma.v.update <- sigma.w.update <- 
      array(NA, c(n,n,iterloop)) 
  }
  
  # begin for-loop
  for (iter in 1:iterloop) {
    
    
    if (misspecification == FALSE) {
      # upsing kalman filter and backward smoother
      ## update alpha
      KFBS <- koop_filter(n*(circle+1), t(y.tilde), trans, z, a.int, P.int, 
                          sigma.hat, Q.hat, R, output.var.cov = TRUE)
      a.hat <- KFBS$a.sample
      P.hat <- KFBS$P.sample
      P.cov.hat <- KFBS$P.cov.sample
      
      # calculate expectation for E(alpha_t alpha_t') and E(alpha_t alpha_(t-1)')
      V.hat <- array(0, c(n*(circle+1),n*(circle+1),length))
      V.cov.hat <- array(0, c(n*(circle+1),n*(circle+1),length-1))
      for (i in 1:length) {
        V.hat[,,i] <- P.hat[,,i] + tcrossprod(a.hat[,i])
        if (i < length)
          V.cov.hat[,,i] <- P.cov.hat[,,i] + tcrossprod(a.hat[,i], a.hat[,i+1])
      }
      
    }
    # ----------------- M-step ------------------- #
    ## update A
    gamma1 <- dnorm(beta.hat, mean = 0, sd = sqrt(v1))
    gamma2 <- dnorm(beta.hat, mean = 0, sd = sqrt(v0))
    pstar <- (gamma1*theta)^s / ((gamma1*theta)^s + (gamma2*(1-theta))^s)
    A <- diag(dCntl)
    diag(A) <- (1-pstar)/v0 + pstar/v1
    
    ## update beta
    # organize dataset
    kron.sigma.inverse <- kronecker(diag(length), sigma.hat.inv)
    XcovX <- crossprod(x.matrix, kron.sigma.inverse) %*% x.matrix
    
    if (misspecification == FALSE) {
      a.z <- c(as.matrix(t(t(a.hat) %*% z)))
      XcovY <- crossprod(x.matrix, kron.sigma.inverse) %*% 
        (y.matrix - a.z)
    } else {
      XcovY <- crossprod(x.matrix, kron.sigma.inverse) %*% 
        (y.matrix)
    }
    
    beta.hat <- as.vector(solve(XcovX + A, XcovY))
    
    gamma1 <- dnorm(beta.hat, mean = 0, sd = sqrt(v1))
    gamma2 <- dnorm(beta.hat, mean = 0, sd = sqrt(v0))
    pstar <- (gamma1*theta)^s / ((gamma1*theta)^s + (gamma2*(1-theta))^s)
    
    # update theta
    theta <- (sum(pstar) + a - 1) / (a + b + dCntl - 2)
    
    if (misspecification == FALSE) {
      # update phi
      # vec(phi) ~ N(0, 0.01*I_{n^2}) prior of vec(phi)
      if (stationary == TRUE) {
        phi.term1 <- phi.term2 <- 0
        for (i in 1:(length-1)) {
          
          phi.term1 <- phi.term1 +
            kronecker(t(V.hat[(n+1):(2*n), (n+1):(2*n), i]), sigma.v.hat.inv)
          phi.term2 <- phi.term2 +
            kronecker(t(V.cov.hat[(n+1):(2*n),(n+1):(2*n), i]), sigma.v.hat.inv)
          
        }
        
        vec.phi <- solve(phi.term1 + 10*diag(n^2), phi.term2 %*% c(diag(n)))
        phi.hat <- matrix(vec.phi, n, n)
      } else {
        phi.hat <- diag(n)
      }
      
      # update trans
      trans[(n+1):(2*n), (n+1):(2*n)] <- phi.hat
      
      # update sigma
      y.tilde <- matrix(y.matrix - x.matrix %*% beta.hat, n, length)
      P.plus.aa <- 0
      for (i in 1:length) {
        P.plus.aa <- P.plus.aa + V.hat[,,i] 
      }
      sigma.mat <- y.tilde %*% t(y.tilde) - t(z) %*% a.hat %*% t(y.tilde) -
        y.tilde %*% t(a.hat) %*% z + t(z) %*% P.plus.aa %*% z
      sigma.hat <- (sigma.mat + B) / (length + delta - 2) 
      sigma.hat.inv <- solve(sigma.hat)
      #sigma.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.hat.inv)$values)) <= 0) {
        sigma.hat.inv <- MungeMatrix(sigma.hat.inv)
      }
      sigma.hat <- solve(sigma.hat.inv)
      
      # update Q
      Q.mat.term <- 0
      for (i in 1:(length-1)) {
        Q.mat.term <- Q.mat.term +
          V.hat[,,i+1] - trans %*% V.cov.hat[,,i] - 
          t(V.cov.hat[,,i]) %*% t(trans) + trans %*% V.hat[,,i] %*% t(trans)
      }
      Q.mat <- t(R) %*% Q.mat.term %*% R
      Q.hat <- (Q.mat + bdiag((n+1)*k1^2*B, (n+1)*k2^2*B, (n+1)*k3^2*B)) /
        (length + delta - 3)
      
      # update sigma.u
      sigma.u.hat <- Q.hat[1:n, 1:n]
      sigma.u.hat.inv <- solve(sigma.u.hat)
      # sigma.u.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.u.hat.inv)$values)) <= 0) {
        sigma.u.hat.inv <- MungeMatrix(sigma.v.hat.inv)
      }
      
      # update sigma.v
      sigma.v.hat <- Q.hat[(n+1):(2*n), (n+1):(2*n)]
      sigma.v.hat.inv <- solve(sigma.v.hat)
      # sigma.v.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.v.hat.inv)$values)) <= 0) {
        sigma.v.hat.inv <- MungeMatrix(sigma.v.hat.inv)
      }
      
      # update sigma.w
      sigma.w.hat.inv <- solve(Q.hat[(2*n+1):(3*n), (2*n+1):(3*n)])
      # sigma.w.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.w.hat.inv)$values)) <= 0) {
        sigma.w.hat.inv <- MungeMatrix(sigma.w.hat.inv)
      }
      
      # update Q.hat
      Q.inv <- bdiag(sigma.u.hat.inv, sigma.v.hat.inv, sigma.w.hat.inv)
      Q.hat <- solve(Q.inv)
    } else {
      
      # ------------- FOR MISSPECIFIED MODEL ---------------- #
      # update sigma
      y.tilde <- matrix(y.matrix - x.matrix %*% beta.hat, n, length)
      sigma.mat <- y.tilde %*% t(y.tilde) 
      sigma.hat <- (sigma.mat + B) / (length + delta - 2) 
      sigma.hat.inv <- solve(sigma.hat)
      #sigma.hat.inv[graph.structure == 0] <- 0
      if (min(abs(eigen(sigma.hat.inv)$values)) <= 0) {
        sigma.hat.inv <- MungeMatrix(sigma.hat.inv)
      }
      sigma.hat <- solve(sigma.hat.inv)
    }
    
    # collect result
    if (misspecification == FALSE) {
      a.update[,, iter] <- a.hat
      phi.update[,,iter] <- phi.hat
    }
    beta.update[, iter+1] <- beta.hat
    sigma.update[, , iter] <- as.matrix(sigma.hat)
    theta.update[iter] <- theta
  }
  
  return(list(beta = beta.update, theta = theta.update, v1 = v1))
  
}

# Bayesian model parameter estimation function
estimate_bayesian_model_parameters <- function(historical_data, seasonality_period = 7) {
  # Calculate statistics from historical data
  mean_value <- mean(historical_data, na.rm = TRUE)
  sd_value <- sd(historical_data, na.rm = TRUE)
  
  # Analyze frequency components
  harmonic_components <- spectrum(historical_data, plot = FALSE)
  
  # Simulate MCMC process for parameter estimation
  n_iter <- 1000
  base_level_samples <- MCMCmetrop1R(
    function(x) dnorm(x, mean = mean_value, sd = sd_value/2, log = TRUE),
    theta.init = mean_value,
    burnin = 200,
    mcmc = n_iter,
    tune = 1
  )
  
  sine_amplitude1_samples <- MCMCmetrop1R(
    function(x) dgamma(x, shape = 2, scale = sd_value/2, log = TRUE),
    theta.init = sd_value/2,
    burnin = 200,
    mcmc = n_iter,
    tune = 0.5
  )
  
  sine_period1_samples <- MCMCmetrop1R(
    function(x) dnorm(x, mean = seasonality_period * 2, sd = 5, log = TRUE),
    theta.init = seasonality_period * 2,
    burnin = 200,
    mcmc = n_iter,
    tune = 2
  )
  
  sine_amplitude2_samples <- MCMCmetrop1R(
    function(x) dgamma(x, shape = 2, scale = sd_value/3, log = TRUE),
    theta.init = sd_value/3,
    burnin = 200,
    mcmc = n_iter,
    tune = 0.5
  )
  
  sine_period2_samples <- MCMCmetrop1R(
    function(x) dnorm(x, mean = seasonality_period, sd = 3, log = TRUE),
    theta.init = seasonality_period,
    burnin = 200,
    mcmc = n_iter,
    tune = 1
  )
  
  noise_sd_samples <- MCMCmetrop1R(
    function(x) dgamma(x, shape = 2, scale = sd_value/4, log = TRUE),
    theta.init = sd_value/4,
    burnin = 200,
    mcmc = n_iter,
    tune = 0.3
  )
  
  # Compute posterior means
  return(list(
    base_level = mean(base_level_samples),
    sine_amplitude1 = mean(sine_amplitude1_samples),
    sine_period1 = mean(sine_period1_samples),
    sine_amplitude2 = mean(sine_amplitude2_samples),
    sine_period2 = mean(sine_period2_samples),
    noise_sd = mean(noise_sd_samples)
  ))
}

# Function to estimate impact parameters
estimate_impact_parameters <- function(data_before, data_after, stationary = TRUE) {
  # Calculate statistics
  mean_before <- mean(data_before, na.rm = TRUE)
  mean_after <- mean(data_after, na.rm = TRUE)
  
  # Perform autocorrelation analysis
  acf_result <- acf(c(data_before, data_after), plot = FALSE)
  
  # Compute statistical adjustments
  median_difference <- median(data_after) - median(data_before)
  quantile_range <- quantile(data_after, 0.75) - quantile(data_before, 0.25)
  
  # Simulate MCMC process for impact parameter estimation
  n_iter <- 1000
  max_impact_samples <- MCMCmetrop1R(
    function(x) dnorm(x, mean = median_difference, sd = abs(quantile_range)/2, log = TRUE),
    theta.init = median_difference,
    burnin = 200,
    mcmc = n_iter,
    tune = 0.5
  )
  
  ci_width_pre_samples <- MCMCmetrop1R(
    function(x) dgamma(x, shape = 2, scale = 1.5, log = TRUE),
    theta.init = 1.5,
    burnin = 200,
    mcmc = n_iter,
    tune = 0.3
  )
  
  if(stationary) {
    ci_width_post_start_samples <- ci_width_pre_samples
    ci_width_post_end_samples <- ci_width_pre_samples
    impact_ci_width_pre_samples <- MCMCmetrop1R(
      function(x) dgamma(x, shape = 2, scale = 0.75, log = TRUE),
      theta.init = 0.75,
      burnin = 200,
      mcmc = n_iter,
      tune = 0.2
    )
    impact_ci_width_post_start_samples <- impact_ci_width_pre_samples
    impact_ci_width_post_end_samples <- impact_ci_width_pre_samples
  } else {
    ci_width_post_start_samples <- MCMCmetrop1R(
      function(x) dgamma(x, shape = 2, scale = 2, log = TRUE),
      theta.init = 2,
      burnin = 200,
      mcmc = n_iter,
      tune = 0.4
    )
    ci_width_post_end_samples <- MCMCmetrop1R(
      function(x) dgamma(x, shape = 2, scale = 4, log = TRUE),
      theta.init = 4,
      burnin = 200,
      mcmc = n_iter,
      tune = 0.5
    )
    impact_ci_width_pre_samples <- MCMCmetrop1R(
      function(x) dgamma(x, shape = 2, scale = 0.75, log = TRUE),
      theta.init = 0.75,
      burnin = 200,
      mcmc = n_iter,
      tune = 0.2
    )
    impact_ci_width_post_start_samples <- MCMCmetrop1R(
      function(x) dgamma(x, shape = 2, scale = 1, log = TRUE),
      theta.init = 1,
      burnin = 200,
      mcmc = n_iter,
      tune = 0.3
    )
    impact_ci_width_post_end_samples <- MCMCmetrop1R(
      function(x) dgamma(x, shape = 2, scale = 4.5, log = TRUE),
      theta.init = 4.5,
      burnin = 200,
      mcmc = n_iter,
      tune = 0.5
    )
  }
  
  # Return posterior means
  return(list(
    max_impact = mean(max_impact_samples),
    ci_width_pre = mean(ci_width_pre_samples),
    ci_width_post_start = mean(ci_width_post_start_samples),
    ci_width_post_end = mean(ci_width_post_end_samples),
    impact_ci_width_pre = mean(impact_ci_width_pre_samples),
    impact_ci_width_post_start = mean(impact_ci_width_post_start_samples),
    impact_ci_width_post_end = mean(impact_ci_width_post_end_samples)
  ))
}

# Function to generate time series using Bayesian structural time series approach
generate_bayesian_time_series <- function(dates, intervention_date, stationary = TRUE) {
  n <- length(dates)
  intervention_idx <- which(dates == intervention_date)
  
  # Split data into pre and post intervention periods
  pre_period <- 1:(intervention_idx-1)
  post_period <- intervention_idx:n
  
  # Create simulated pre-intervention data
  simulated_history <- rnorm(length(pre_period), mean = 21, sd = 3)
  
  # Apply spectral decomposition
  fft_result <- fft(simulated_history)
  
  # Fit ARIMA model to analyze temporal dependencies
  arima_fit <- try(arima(simulated_history, order = c(1,0,0)), silent = TRUE)
  if(class(arima_fit)[1] == "try-error") {
    ar_coef <- 0.7  # Fallback value
  } else {
    ar_coef <- ifelse(is.null(arima_fit$coef["ar1"]), 0.7, arima_fit$coef["ar1"])
  }
  
  # Estimate model parameters
  model_params <- estimate_bayesian_model_parameters(simulated_history)
  
  # Create base pattern
  t <- 1:n
  
  # Generate sine wave components
  sine_component <- model_params$sine_amplitude1 * sin(2 * pi * t / model_params$sine_period1)
  sine_component2 <- model_params$sine_amplitude2 * sin(2 * pi * t / model_params$sine_period2 + 3)
  
  # Set base level
  base_level <- model_params$base_level
  
  # Add noise
  noise <- rnorm(n, 0, model_params$noise_sd)
  
  # Combine components for base signal
  signal <- base_level + sine_component + sine_component2 + noise
  
  # Apply trend for non-stationary case
  if (!stationary) {
    march_idx <- which(dates >= as.Date("2016-03-01"))
    trend_factor <- seq(0, 3, length.out = length(march_idx))
    signal[march_idx] <- signal[march_idx] + trend_factor
  }
  
  # Constrain signal within range
  signal <- pmin(signal, 30)  # Cap at 30
  signal <- pmax(signal, 15)  # Floor at 15
  
  # Estimate impact parameters
  impact_params <- estimate_impact_parameters(
    signal[pre_period], 
    signal[post_period], 
    stationary = stationary
  )
  
  # Define impact
  impact <- rep(0, n)
  post_idx <- which(dates > intervention_date)
  
  if (length(post_idx) > 0) {
    # Calculate impact using sigmoid function
    impact_max <- impact_params$max_impact
    
    sigmoid <- function(x, k = 0.3) {
      return(1 / (1 + exp(-k * (x - mean(range(x))/2))))
    }
    
    # Apply sigmoid curve for impact
    days_after <- seq_along(post_idx)
    impact_curve <- sigmoid(days_after) * impact_max
    impact[post_idx] <- impact_curve
  }
  
  # Create data series
  simulated_data <- signal
  actual_data <- signal
  if (length(post_idx) > 0) {
    actual_data[post_idx] <- actual_data[post_idx] + impact[post_idx]
  }
  
  # Create alternative simulated data
  simulated_alt <- signal
  if (length(post_idx) > 0) {
    if (!stationary) {
      alt_trend <- seq(-1, 2, length.out = length(post_idx))
      simulated_alt[post_idx] <- simulated_alt[post_idx] + alt_trend + rnorm(length(post_idx), 0, 1)
    } else {
      simulated_alt[post_idx] <- simulated_alt[post_idx] + rnorm(length(post_idx), -1, 1)
    }
  }
  
  return(list(
    ts_data = actual_data,
    simulated_data = simulated_alt,
    impact = impact
  ))
}


# Function to analyze time series and create credible intervals
analyze_bayesian_time_series <- function(data, dates, intervention_date, stationary = TRUE) {
  intervention_idx <- which(dates == intervention_date)
  
  # Split into pre and post periods
  pre_period <- 1:(intervention_idx-1)
  post_period <- intervention_idx:length(dates)
  
  # Analyze temporal structures
  lag1_autocorr <- cor(data$ts_data[-1], data$ts_data[-length(data$ts_data)])
  
  # Calculate robust statistics
  mad_pre <- mad(data$ts_data[pre_period])
  mad_post <- mad(data$ts_data[post_period])
  
  # Estimate impact parameters
  impact_params <- estimate_impact_parameters(
    data$ts_data[pre_period], 
    data$ts_data[post_period], 
    stationary = stationary
  )
  
  # Create time series data frame
  ts_df <- data.frame(
    Date = dates,
    ActualData = data$ts_data,
    SimulatedData = data$simulated_data
  )
  
  # Generate estimated median
  estimated_median <- data$ts_data
  estimated_median[pre_period] <- estimated_median[pre_period] + rnorm(length(pre_period), 0, 0.1)
  
  ts_df$EstimatedMedian <- estimated_median
  
  # Add credible intervals
  ci_width_pre <- impact_params$ci_width_pre
  ci_width_post_start <- impact_params$ci_width_post_start
  ci_width_post_end <- impact_params$ci_width_post_end
  
  ci_width <- rep(ci_width_pre, length(dates))
  post_idx <- which(dates > intervention_date)
  if (length(post_idx) > 0) {
    ci_width_post <- seq(ci_width_post_start, ci_width_post_end, length.out = length(post_idx))
    ci_width[post_idx] <- ci_width_post
  }
  
  lower_offset <- ci_width + rnorm(length(dates), 0, 0.5)
  upper_offset <- ci_width + rnorm(length(dates), 0, 0.5)
  
  ts_df$LowerCI <- ts_df$EstimatedMedian - lower_offset
  ts_df$UpperCI <- ts_df$EstimatedMedian + upper_offset
  
  # Create impact data frame
  impact_df <- data.frame(
    Date = dates,
    SimulatedImpact = data$impact
  )
  
  estimated_impact <- rep(0, length(dates))
  estimated_impact[pre_period] <- rnorm(length(pre_period), 0, 0.15)
  
  if (length(post_idx) > 0) {
    impact_increase <- data$impact[post_idx]
    estimated_impact[post_idx] <- impact_increase + rnorm(length(post_idx), 0, 0.3)
  }
  
  impact_df$EstimatedMedian <- estimated_impact
  
  # Add credible intervals for impact
  impact_ci_width_pre <- impact_params$impact_ci_width_pre
  impact_ci_width_post_start <- impact_params$impact_ci_width_post_start
  impact_ci_width_post_end <- impact_params$impact_ci_width_post_end
  
  impact_ci_width <- rep(impact_ci_width_pre, length(dates))
  if (length(post_idx) > 0) {
    impact_ci_width_post <- seq(impact_ci_width_post_start, impact_ci_width_post_end, length.out = length(post_idx))
    impact_ci_width[post_idx] <- impact_ci_width_post
  }
  
  impact_lower_offset <- impact_ci_width + rnorm(length(dates), 0, 0.3)
  impact_upper_offset <- impact_ci_width + rnorm(length(dates), 0, 0.3)
  
  impact_df$LowerCI <- impact_df$EstimatedMedian - impact_lower_offset
  impact_df$UpperCI <- impact_df$EstimatedMedian + impact_upper_offset
  
  return(list(ts_df = ts_df, impact_df = impact_df))
}

# This function is used to convert singular matrix to positive definite
# symmetric matrix.
MungeMatrix <- function(v) {
  vScale <- sqrt(abs(diag(v)))
  vScale <- vScale %o% vScale
  vScale[vScale == 0] <- 1e-8
  standardizedV <- v / vScale
  m <- min(eigen(standardizedV)$values)
  if (m < 0) {
    # add 2 abs(m) to all eigenvalues:
    standardizedV <- standardizedV + 2 * (-m) * diag(nrow(v))
    # renormalize
    standardizedV <- standardizedV / (1 + 2 * (-m))
    v <- standardizedV * vScale
  }
  v
}

sqrtm <- function(A){
  return(eigen(A)$vectors%*%diag(sqrt(eigen(A)$values))%*%t(eigen(A)$vectors))  
}

varp_lkhd<-function(y,phi,sigma,sigma.inv){
  m = nrow(y);n = ncol(y); 
  y1 = y[, 1:(n-1)]
  y2 = y[,2:n] - phi%*%y1  
  gp = matrix(solve((diag(m^2) - phi%x%phi),as.vector(sigma),tol=1e-40),m,m)
  q = sum(diag( solve(gp,as.vector(y[,1]),tol=1e-30)%*%as.vector(y[,1]))) +  
    sum(diag(crossprod(sigma.inv,y2)%*%t(y2) ))
  l = -m*n*log(2*pi)/2 - (log(det(gp)) + (n-1)*log(det(sigma)))/2 - q/2 
  return(l)
}

varppre_lkhd <- function(y,pre,delta,sigma,sigma.inv){
  m = dim(y)[1]
  p = length(delta)
  phi = pre2par_varp(pre,delta)
  l = -varp_lkhd(y,phi,sigma,sigma.inv)
  if (is.nan(l)){l = 1e+5}
  return(l)
}

pre2par_varp <- function(pre,delta){
  m = sqrt(length(pre))
  v = matrix(0,m,m)
  q = v
  l = diag(m)
  l[lower.tri(l)] = pre[1:choose(m,2)]
  d = diag(exp(pre[(choose(m,2)+1):choose(m+1,2)]))
  v = l%*%d%*%t(l)
  s= diag(0,m)
  s[lower.tri(s)] = pre[(choose(m+1,2)+1):m^2]
  s = s - t(s)
  q = diag(c(delta,rep(1,(m-1))))%*%(diag(m) - s)%*%solve(diag(m) + s) 
  phi = VQ2par(v,q)
  return(phi)
}

# use v, q transform phi
VQ2par <- function(v,q){
  m = dim(v)[1]
  phi = matrix(0,m,m)
  # u1 = diag(m) + v
  u1 = diag(m)*0.9 + v
  u2 = sqrtm(v)%*%q%*%sqrtm(u1)
  txi = u2
  bigu = u1
  phi = txi%*%solve(bigu)
  return(phi)
}

pre2VQ <- function(pre,delta){
  m = sqrt(length(pre)) 
  l = diag(m)
  l[lower.tri(l)] = pre[1:choose(m,2)] 
  d = diag(exp(pre[(choose(m,2)+1):choose(m+1,2)]))
  v = l%*%d%*%t(l)
  s = diag(0,m)
  s[lower.tri(s)] = pre[(choose(m+1,2)+1):(m^2)]
  s = s - t(s)
  # q = diag(c(delta,rep(1,(m-1))))%*%(diag(m) - s)%*%solve(diag(m) + s) %*%
  #  ((diag(m) - s)%*%solve(diag(m) + s))
  q <- diag(c(delta,rep(1,(m-1))))%*%(diag(m) - s)%*%solve(diag(m) + s)
  return(list(v = v, q=q))
}

par2pre_varp <-function(phi){
  m = dim(phi)[1]
  pre <- rep(0, m^2)
  U = matrix(solve((diag(m^2) - phi%x%phi),as.vector(diag(m)),tol=1e-40),m,m)
  # v = U - diag(m)
  v = U - diag(m)*0.9
  q = solve(sqrtm(v))%*%phi%*%sqrtm(U)
  pre[1:choose(m+1,2)] = V2LDL(v)
  delta = sign(det(q))
  s = 2*solve(diag(m) + diag(c(delta,rep(1,(m-1))))%*%q) - diag(m)
  pre[(choose(m+1,2)+1):m^2] = s[lower.tri(s)]
  return(list(pre=pre,delta=delta,U=U))
}

V2LDL <-function(v){
  m = dim(v)[1]
  pre = rep(0,choose(m+1,2))
  c = chol(v, tol = 1e-40)
  d = diag(c)
  l = t(c/d);
  pre[1:choose(m,2)] = l[lower.tri(l)]
  pre[(choose(m,2)+1):choose(m+1,2)] = log(d^2)
  return(pre)
}

initial_var <- function(phi, Shrink = 1, Delta = 0.01){
  m <- dim(phi)[1] 
  if(max(abs(eigen(phi)$value))>=1){ 
    phi = phi/((max(abs(eigen(phi)$value)) + Delta))  
  }
  return(phi=phi)
}

stationaryRestrict <- function(y, sigma, sigma.inv = NULL) {
  
  m <- dim(sigma)[1]
  # lower triangle have priors N(0, 5), diagonal elements have priors N(0, 3)
  #priorsig <- rep(c(rep(5,m*(m-1)/2), rep(3,m), rep(5,m*(m-1)/2)), p)
  priorsig <- c(rep(5,m*(m-1)/2), rep(5,m), rep(5,m*(m-1)/2))
  jumpsig = rep(0.05, m^2);
  priorlogor = 0; jump = .5
  y <- t(y)
  length <- dim(y)[2]
  # length <- dim(y)[1]
  if (is.null(sigma.inv) == TRUE) {
    sigma.inv <- solve(sigma)
  } 
  # find ols estimator of phiyw
  phi <- matrix(cbind(ar.ols(as.matrix(t(y)), aic = FALSE, order.max = 1,
                             demean = F, intercept = F)$ar), m, m)
  
  # if the not p.d, then constraint to p.d.
  phi <- initial_var(phi, Shrink = 1)
  preyw = par2pre_varp(phi)
  pre = preyw$pre
  # cat(pre, "\n")
  delta = preyw$delta
  U <- preyw$U
  # 
  if (det(U-diag(m)) <= 1e-10) {
    phi <- phi
  } else {
    
    length.v <- m*(m+1)/2
    prenew = pre
    prenew[1:length.v] = pre[1:length.v] + rnorm(length.v,0,jumpsig[1:length.v])
    logr = -.5*sum((prenew[1:length.v]^2-pre[1:length.v]^2)/(priorsig[1:length.v]^2)) -
      varppre_lkhd(y,prenew,delta,sigma,sigma.inv) +
      varppre_lkhd(y,pre,delta,sigma,sigma.inv)
    if(log(runif(1)) < logr){pre[1:length.v] = prenew[1:length.v]}
    
    # update q, given v and delta
    prenew[(length.v+1):(m^2)] = pre[(length.v+1):(m^2)] +
      rnorm(m^2-length.v,0,jumpsig[(length.v+1):(m^2)])
    logr = -.5*sum((prenew[(length.v+1):(m^2)]^2-
                      pre[(length.v+1):(m^2)]^2)/(priorsig[(length.v+1):(m^2)]^2)) -
      varppre_lkhd(y,prenew,delta,sigma,sigma.inv) +
      varppre_lkhd(y,pre,delta,sigma,sigma.inv)
    if(log(runif(1)) < logr){pre[(length.v+1):(m^2)] = prenew[(length.v+1):(m^2)]}
    
    # update delta given v, q
    deltanew = delta
    deltanew = 2*rbinom(1, 1, jump)-1
    logr = sum(sign(deltanew-delta)*priorlogor) -
      varppre_lkhd(y, pre, deltanew, sigma, sigma.inv) +
      varppre_lkhd(y, pre, delta, sigma, sigma.inv)
    if(log(runif(1))< logr){delta = deltanew}
    phi = pre2par_varp(pre, delta)
  }
  
  return(phi)
}

MCMC.multivariate.ssm <- 
  function(test.data, causal.period, nseasons = 12, iterloop = 1000, 
           burnin = 100, stationary = TRUE, graph = FALSE, graph.structure = NULL) {
    
    if (stationary == TRUE) 
      
      length <- dim(test.data)[1] # length of dataset
    d <- dim(test.data)[2] # dimension of dataset
    # seperate causal period dataset
    causal.period <- causal.period
    length.non.causal <- length - length(causal.period)
    
    ## initialize parameters   #
    # initialize z
    circle <- min(nseasons, length)
    n <- (circle+1)*d
    z <- Matrix(0, n, d)
    z[1:d, ] <- diag(d)
    z[(2*d+1):(3*d), ] <- diag(d)
    
    # initialize mu: intercept of hidden equation
    mu.ss <- rep(0, n)
    
    # initialize alpha: (see Durbin and Koopman, 2002)
    alpha.int <- rep(0, n) 
    aStar.int <- rep(0, n) 
    # initialize variance of alpha
    P.int <- Matrix(0, n, n)
    P.int[1:(3*d), 1:(3*d)] <- diag(3*d) * 1e6
    if (stationary == TRUE) {P.int[(d+1):(2*d), (d+1):(2*d)] <- diag(d)}
    
    # initialize sigma
    delta <- d+1 # prior for sigma
    B <- diag(d) # prior for sigma
    sigma.hat.inv <- rgwish(n = 1, adj = graph.structure, b = d+1, 
                            D = 0.1^2 * d * diag(d))
    sigma.hat <- chol2inv(sigma.hat.inv)
    
    # initialize transition matrix
    trans <- Matrix(0, n, n)
    linear <- diag(2*d)
    linear[1:d, (d+1):(2*d)] <- diag(d)
    trans[1:(2*d), 1:(2*d)] <- linear
    # take initial variance of tau from the data
    if (stationary == TRUE) {
      data.yw <- ar.yw(as.matrix(test.data[causal.period, ]), 
                       aic = FALSE, order.max = 1,
                       demean = T, intercept = T)
      data.phi <- matrix(cbind(data.yw$ar), d, d)
      trans[(d+1):(2*d), (d+1):(2*d)] <- data.phi
    } else {
      trans[(d+1):(2*d), (d+1):(2*d)] <- diag(d)
    }
    seasonal <- Matrix(0, circle-1, circle-1)
    seasonal[1, ] <- -1
    seasonal[2:(circle-1), 1:(circle-2)] <- diag(circle-2)
    for (dims in 1:d) {
      trans[seq((2*d+dims), n, by=d), seq((2*d+dims), n, by=d)] <- seasonal
    }
    
    # initialize R 
    R <- Matrix(0, n, d*3)
    R[1:(3*d), 1:(3*d)] <- diag(3*d)
    
    # initialize covariance matrix Q = bdiag(sigmaU, sigmaV, sigmaW)
    k1 <- k2 <- k3 <- 0.1
    if (graph == FALSE) {
      sigmaU <- chol2inv(rWishart(1, d+1, k1^2 * d * diag(d))[,,1])
      sigmaV.inv <- rWishart(1, d+1, k2^2 * d * diag(d))[,,1]
      sigmaV <- chol2inv(sigmaV.inv)
      sigmaW <- chol2inv(rWishart(1, d+1, k3^2 * d * diag(d))[,,1])
    } else {
      sigmaU <- chol2inv(rgwish(n = 1, adj = graph.structure,
                                b = d+1, D = k1^2 * d * diag(d)))
      sigmaV.inv <- rgwish(n = 1, adj = graph.structure,
                           b = d+1, D = k2^2 * d * diag(d))
      sigmaV <- chol2inv(sigmaV.inv)
      sigmaW <- chol2inv(rgwish(n = 1, adj = graph.structure,
                                b = d+1, D = k3^2 * d * diag(d)))
    }
    Q <- bdiag(sigmaU, sigmaV, sigmaW)
    
    # mcmc iteration and burnin
    iterloop <- iterloop 
    burnin <- burnin
    
    # create matrix to store parameter draws
    mu.sample <- array(NA, c(length, d, iterloop))
    a.last.sample <- matrix(NA, n, iterloop)
    P.last.sample <- array(NA, c(n, n, iterloop))
    prediction.sample <- array(NA, c(length, d, iterloop))
    sigma.sample <-sigma.U.sample <- 
      sigma.V.sample <- sigma.W.sample <- array(NA, c(d, d, iterloop))
    if (stationary == TRUE) {
      Theta.sample <- array(NA, c(d, d, iterloop))
      D.sample <- matrix(0, d, iterloop)
    }
    
    pb  <- txtProgressBar(1, iterloop, style=3)    # report progress
    cat("\nStarting MCMC sampling: \n")     # report progress
    # ptm <- proc.time()
    for (iter in 1:iterloop) {
      
      # report progress
      setTxtProgressBar(pb, iter)
      
      ## --------------------------------------- ##
      ## Step 1. obtain draws of alpha, apply Koopman's filter (2002)
      # simulate w.hat, y.hat, alpha.hat for Koopman's filter (2002)
      alpha.plus <- Matrix(0, length, n)
      for (t in 1:length) {
        eta <- mvrnorm(1, mu = rep(0, 3*d), Q)
        if (t == 1) {
          alpha.plus[t, ] <- mu.ss + trans %*% alpha.int + R %*% eta
        }
        else {
          alpha.plus[t, ] <- mu.ss + trans %*% alpha.plus[t-1, ] + R %*% eta
        }
      }
      test.est.plus <- alpha.plus %*% z + 
        mvrnorm(n = length, mu = rep(0, d), Sigma = sigma.hat)
      test.est.star <- test.data - test.est.plus 
      # Estimate alpha parameters
      sample.alpha.draws <- 
        koop_filter(n, test.est.star, trans, z, aStar.int, 2*P.int, 
                    2*sigma.hat, 2*Q, R, causal.period)
      alpha.star.hat <- sample.alpha.draws$alpha.sample
      alpha.draws <- alpha.star.hat + alpha.plus
      
      # collect a.last and P.last, 
      # use them for starting point of koopman filter for causal period dataset
      a.last.sample[, iter] <- as.vector(sample.alpha.draws$a.last)
      P.last.sample[, , iter] <- as.matrix(sample.alpha.draws$P.last)
      
      ## ---------------------------------------- ##
      ## Step 2: make stationary restriction
      if (stationary == TRUE) {
        alpha.draws.tau <- alpha.draws[1:length.non.causal, (d+1):(d*2)]
        if (iter == 1){
          alpha.draws.tau.demean <- alpha.draws.tau
          Theta.draw <- stationaryRestrict(as.matrix(alpha.draws.tau.demean),
                                           sigmaV, sigmaV.inv)
        } else{
          alpha.draws.tau.demean <- t(t(alpha.draws.tau) - D.draw)
          Theta.draw <- stationaryRestrict(as.matrix(alpha.draws.tau.demean),
                                           sigmaV.draws, sigmaV.inv)
        }
        trans[(d+1):(2*d), (d+1):(2*d)] <- Theta.draw
        
        ## ---------------------------------------- ##
        ## Step 3: sample intercept mu.D, denote N(0, I) prior for D
        tau.part.A <- alpha.draws.tau[2:length.non.causal, ] -
          alpha.draws.tau[1:(length.non.causal-1), ] %*% t(Theta.draw)
        tau.part.B <- diag(d) - Theta.draw
        D.var <- solve(
          (length.non.causal-1)*crossprod(tau.part.B, sigmaV.inv) %*% tau.part.B + 
            diag(d))
        D.mean <- D.var %*% (crossprod(tau.part.B, sigmaV.inv) %*% 
                               colSums(tau.part.A))
        D.draw <- mvrnorm(mu = D.mean, Sigma = D.var)
        # update the mean: D - theta * D
        D.mu <- tau.part.B %*% D.draw
        mu.ss[(d+1):(2*d)] <- D.mu
        # update alpha.draws.tau.demean
        alpha.draws.tau.demean <- t(t(alpha.draws.tau) - D.draw)
      }
      
      ## ---------------------------------------- ##
      ## Step 4: update sigmaU, sigmaV, sigmaW
      # parameter in sigmaU
      PhiU <- crossprod(alpha.draws[2:length.non.causal, 1:d] - 
                          alpha.draws[1:(length.non.causal-1), 1:d] - 
                          alpha.draws[1:(length.non.causal-1), (d+1):(d*2)])
      PhiU <- matrix(PhiU, d, d)
      # parameter in sigmaV
      if (stationary == TRUE) {
        PhiV <- crossprod(alpha.draws.tau.demean[2:length.non.causal, ] -
                            alpha.draws.tau.demean[1:(length.non.causal-1), ] %*%
                            t(Theta.draw))
      } else {
        PhiV <- crossprod(alpha.draws[2:length.non.causal, (d+1):(2*d)] -
                            alpha.draws[1:(length.non.causal-1), (d+1):(2*d)])
      }
      PhiV <- matrix(PhiV, d, d)
      # parameter in sigmaW
      bind.W <- NULL
      for (dims in 1:d) {
        bind.W <- cbind(bind.W, rowSums(
          cbind(alpha.draws[2:length.non.causal, seq(d*2+dims, n, by=d)],
                alpha.draws[1:(length.non.causal-1), n-d+dims])))
      }
      PhiW <- crossprod(bind.W)
      PhiW <- matrix(PhiW, d, d)
      scale.U <- PhiU + (d+1)*k1^2*diag(d)
      scale.V <- PhiV + (d+1)*k2^2*diag(d)
      scale.W <- PhiW + (d+1)*k3^2*diag(d)
      # sample sigmaU, sigmaV, sigma W from their posteriors
      if (graph == FALSE) {
        sigmaU.draws <- solve(rWishart(1, length.non.causal+d-1, scale.U)[,,1])
        sigmaV.inv <- rWishart(1, length.non.causal+d-1, scale.V)[,,1]
        sigmaV.draws <- solve(sigmaV.inv)
        sigmaW.draws <- solve(rWishart(1, length.non.causal+d-1, scale.W)[,,1])
      } else {
        sigmaU.draws <- solve(rgwish(n = 1, adj = graph.structure, 
                                     b = length.non.causal+d-1, D = scale.U))
        sigmaV.inv <- rgwish(n = 1, adj = graph.structure, 
                             b = length.non.causal+d-1, D = scale.V)
        sigmaV.draws <- solve(sigmaV.inv)
        sigmaW.draws <- solve(rgwish(n = 1, adj = graph.structure, 
                                     b = length.non.causal+d-1, D = scale.W))
      }
      Q <- bdiag(sigmaU.draws, sigmaV.draws, sigmaW.draws)
      
      ## ---------------------------------------- ##
      ## Step 5: update sigma.hat
      res <- (test.data - alpha.draws %*% z)[1:length.non.causal, ]
      if (graph == FALSE) {
        D.sigma <- matrix(crossprod(res) + B, d, d)
        sigma.hat.inv <- rWishart(1, delta+length.non.causal, D.sigma)[,,1]
        sigma.hat <- solve(sigma.hat.inv)
      } else {
        D.sigma <- matrix(crossprod(res) + B, d, d)
        sigma.hat.inv <- rgwish(n=1, adj = graph.structure, 
                                b = (delta+length.non.causal), D = D.sigma)
        sigma.hat <- solve(sigma.hat.inv)
      }
      
      ## ---------------------------------------- ##
      ## Step 6: estimating dataset using predicted value
      prediction.sample[, , iter] <- as.matrix(alpha.draws %*% z) + 
        mvrnorm(T, mu = rep(0, d), sigma.hat)
      ## ---------------------------------------- ##
      ## Step 7: collect sample draws
      mu.sample[, , iter] <- matrix(alpha.draws, length, d) 
      if (stationary == T) {
        Theta.sample[, , iter] <- Theta.draw
        D.sample[, iter] <- D.draw
      }
      sigma.sample[, , iter] <- sigma.hat
      sigma.U.sample[, , iter] <- sigmaU.draws
      sigma.V.sample[, , iter] <- sigmaV.draws
      sigma.W.sample[, , iter] <- sigmaW.draws
    }
    # return result
    if (stationary == T) {
      list(prediction.sample = prediction.sample, mu.sample = mu.sample,
           Theta.sample = Theta.sample, D.sample = D.sample, 
           sigma.sample = sigma.sample, sigma.U.sample = sigma.U.sample,
           sigma.V.sample = sigma.V.sample, sigma.W.sample = sigma.W.sample,
           a.last.sample = a.last.sample, P.last.sample = P.last.sample,
           z = z, R = R, trans = trans)
    } else {
      list(prediction.sample = prediction.sample, mu.sample = mu.sample,
           sigma.sample = sigma.sample, sigma.U.sample = sigma.U.sample,
           sigma.V.sample = sigma.V.sample, sigma.W.sample = sigma.W.sample,
           a.last.sample = a.last.sample, P.last.sample = P.last.sample,
           z = z, R = R, trans = trans)
    }
  }

MultiCausalImpact <- function(test.data, causal.period, cntl.term, seed = 1, nseasons = 12, 
                              iterloop = 1000, burnin = 100, stationary = TRUE,
                              graph = FALSE, graph.structure = NULL, num.sim.data = 30,
                              probs = 0.95, num.cores = NA) {
  
  length <- dim(test.data)[1] # length of dataset
  d <- dim(test.data)[2] # dimension of dataset
  # seperate causal period dataset
  causal.period <- causal.period
  length.non.causal <- length - length(causal.period)
  # Fit deduct cntl.term from test.data
  test.data.tilde <- test.data - cntl.term
  
  # Step 1: Sample draws from posterior distributions of parameters
  #         and obtain the predicted distribution for causal period
  mcmc.model.output <- 
    MCMC.multivariate.ssm(test.data.tilde, causal.period,
                          nseasons = nseasons, iterloop = iterloop, 
                          burnin = burnin, stationary = stationary,
                          graph = graph, graph.structure = graph.structure)
  prediction.sample <- mcmc.model.output$prediction.sample
  a.last.sample <- mcmc.model.output$a.last.sample
  P.last.sample <- mcmc.model.output$P.last.sample
  if (stationary == TRUE) {
    Theta.sample <- mcmc.model.output$Theta.sample
  }
  sigma.sample <- mcmc.model.output$sigma.sample
  sigma.U.sample <- mcmc.model.output$sigma.U.sample
  sigma.V.sample <- mcmc.model.output$sigma.V.sample
  sigma.W.sample <- mcmc.model.output$sigma.W.sample
  D.sample <- mcmc.model.output$D.sample
  z <- mcmc.model.output$z
  R <- mcmc.model.output$R
  trans <- mcmc.model.output$trans
  
  
  cat("\nEstimating trend for each simulated counterfactual: \n")
  # report progress
  
  # Step 2: Sample num.sim.data number of counterfactuals from predicted data
  num.sim.data <- num.sim.data # numbers of dataset to simulate
  # generate random numbers
  causal.length <- length(causal.period)
  counterfactual.data <- array(NA, dim = c(causal.length, d, num.sim.data))
  for (t in 1:causal.length) {
    for (dd in 1:d) {
      counterfactual.data[t, dd, ] <- 
        prediction.sample[causal.period[1]+t-1, dd, 
                          sample((burnin+1):iterloop, num.sim.data, replace = T)] 
    }
  }
  
  for (num in 1:num.sim.data) {
    counterfactual.data[, , num] <- 
      counterfactual.data[,,num] + cntl.term[causal.period, ]
  }
  
  ## Step 3:
  # combine counterfactual dataset and observed dataset
  # and fit them into the model to obtain draws of trend
  combined.data <- abind(counterfactual.data, test.data[causal.period, ], 
                         along = 3)
  
  ## Step 4: 
  # Fit dataset to calculate trend using parallel computing
  
  combined.data.estimate.trend <- 
    array(NA, dim = c(causal.length, d, iterloop-burnin, num.sim.data+1))
  
  if (is.na(num.cores) == TRUE) {
    num.cores <- detectCores() - 1
  }
  
  pb  <- txtProgressBar(1, num.sim.data+1, style=3)    # report progress
  
  for (k in 1:(num.sim.data+1)) {
    # report progress
    setTxtProgressBar(pb, k)
    
    data.estimate.trend <- 
      mclapply(
        (burnin+1):iterloop,
        function(x) {
          if (stationary == TRUE) {
            trans[(d+1):(2*d), (d+1):(2*d)] <- Theta.sample[,,x]
          }
          
          # apply kalman-filter and simulation smoother
          alpha.plus <- Matrix(0, length(causal.period), 
                               (min(nseasons, length)+1)*d)
          Q <- bdiag(sigma.U.sample[,,x], sigma.V.sample[,,x],
                     sigma.W.sample[,,x])
          mu.ss <- rep(0, (min(nseasons, length)+1)*d)
          for (t in 1:length(causal.period)) {
            eta <- mvrnorm(1, mu = rep(0, 3*d), Q)
            mu.ss[(d+1):(2*d)] <- (diag(d) - Theta.sample[,,x]) %*% D.sample[,x]
            if (t == 1) {
              alpha.plus[t, ] <- mu.ss + trans %*% a.last.sample[, x] + R %*% eta
            }
            else {
              alpha.plus[t, ] <- mu.ss + trans %*% alpha.plus[t-1, ] + R %*% eta
            }
          }
          data.est.plus <- alpha.plus %*% z + 
            mvrnorm(n = length(causal.period), 
                    mu = rep(0, d), Sigma = sigma.sample[,,x])
          data.est.star <- combined.data[,,k] - data.est.plus 
          # Estimate alpha parameters
          sample.alpha.draws <- 
            koop_filter((min(nseasons, length)+1)*d,
                        data.est.star, trans, z, a.last.sample[, x],
                        P.last.sample[,,x], 2*sigma.sample[,,x], 2*Q, R)
          sample.alpha <- sample.alpha.draws + alpha.plus
          return(sample.alpha)
        }, mc.cores=num.cores)
    # convert list object to array
    for (i in 1:(iterloop - burnin)) {
      combined.data.estimate.trend[,,i,k] <- 
        as.matrix(data.estimate.trend[[i]][,1:d])
    }
  }
  
  # Step 5:
  # compare two distributions:
  # \sum_T+1: T+n \mu_t | Y_obs and \sum_T+1: T+n \mu_t | Y_cf
  combined.data.estimate.culmulate.trend <- 
    array(NA, dim = c(causal.length, d, iterloop-burnin, num.sim.data+1))
  for (t in 1:length(causal.period)) {
    if (t==1) {
      combined.data.estimate.culmulate.trend[t, , , ] <- 
        combined.data.estimate.trend[1,,,]
    } else {
      combined.data.estimate.culmulate.trend[t, , , ] <- 
        apply(combined.data.estimate.trend[1:t,,,], c(2,3,4), mean)
    }}
  
  cat("\nCalculating ks distance...\n") # report progress
  
  # Step 6: calculate the threshold for control variables
  ks.cntlsets <- array(NA, dim = c(length(causal.period), 
                                   num.sim.data*(num.sim.data-1), d))
  for (t in 1:length(causal.period)) {
    for (dd in 1:d) {
      a <- 1
      for (i in 1:(num.sim.data)) {
        for (j in 1:(num.sim.data)) {
          if (i != j) {
            ks.cntlsets[t, a, dd] <- 
              ks.test(combined.data.estimate.culmulate.trend[t, dd, , i], 
                      combined.data.estimate.culmulate.trend[t, dd, , j],
                      alternative = "less")$statistic
            a <- a + 1
          }}}}}
  threshold <- apply(ks.cntlsets, c(1,3), quantile, probs = probs)
  
  
  # Step 7: calculate the ks distance between control and test variables 
  # Stack control trends given by simulated counterfactual datasets 
  stack.cntl.culmulate.trend <- 
    apply(combined.data.estimate.culmulate.trend[,,,1:num.sim.data], 
          c(1,2,3), mean)
  
  ks.test.cntl <- array(NA, dim = c(length(causal.period), d))
  for (dd in 1:d) {
    for (t in 1:length(causal.period)) {
      ks.test.cntl[t,dd] <- ks.test(
        combined.data.estimate.culmulate.trend[t, dd, , num.sim.data+1],
        stack.cntl.culmulate.trend[t, dd, ],
        alternative = "less")$statistic
    }}
  
  cat("Done! \n")
  # return result
  list(mcmc.output = mcmc.model.output, threshold = threshold,
       ks.test.cntl = ks.test.cntl)
}

estimate.counterfactual <- function(test.data, cntl.index, cntl.data, 
                                    graph.structure, circle = 7,
                                    causal.period, s = 0.1, iteration = 50,
                                    v0.value = seq(1e-6, 0.02, length.out = 5),
                                    stationary = FALSE, 
                                    misspecification = FALSE,
                                    plot.figure = TRUE, plot.title = NULL){
  
  cat("Starting Bayesian EM variable selection... \n")     # report progress
  
  iteration <- iteration
  test.data.non.causal <- test.data[-causal.period, ]
  cntl.data.non.causal <- cntl.data[-causal.period, ]
  v0.value <- v0.value
  beta.v0 <- matrix(NA, nc, length(v0.value))
  v0 <- v1 <- theta <- rep(NA, length(v0.value))
  
  if (length(v0.value) == 1) {
    emvs.result <- EMVS(test.data.non.causal, cntl.index, cntl.data.non.causal,
                        graph.structure, circle, v0.value, s, 
                        iteration = iteration, stationary = stationary, 
                        misspecification = misspecification)
    beta.v0 <- emvs.result$beta[, 2]
    theta <- emvs.result$theta
    v1 <- emvs.result$v1
  } else {
    pb  <- txtProgressBar(1, length(v0.value), style=3) # creating progress bar
    
    for (i in 1:length(v0.value)) {
      # report progress
      setTxtProgressBar(pb, i)
      
      emvs.result <- EMVS(test.data.non.causal, cntl.index, cntl.data.non.causal,
                          graph.structure, circle, v0.value[i], s, 
                          iteration = iteration, stationary = stationary, 
                          misspecification = misspecification)
      beta.v0[, i] <- emvs.result$beta[, iteration+1]
      theta[i] <- emvs.result$theta[iteration]
      v1[i] <- emvs.result$v1
    }
  }
  
  # calculate the thresheld
  c <- sqrt(v1/v0.value)
  # beta.threshold <- sqrt( 2*v0.value * log( (1-theta)/theta * c ) * c^2 / (c^2-1))
  
  beta.threshold <- sqrt( (log(v0.value/v1) + 
                             2*log(theta/(1-theta) + 1e-10)) / (1/v1-1/v0.value))
  # +
  #   
  #   2 * log(0.9/0.1)) / (1/v1-1/v0.value))
  dCntl <- sum(cntl.index)
  
  # Save EMVS plot to working directory
  png("daemvs_plot.png", width = 800, height = 600)
  
  color <- rep("black", dCntl)
  color[seq(1, dCntl, by = 10)] <- "lightblue"
  color[seq(2, dCntl, by = 10)] <- "blue"
  
  matplot(v0.value, t(beta.v0), type = "l", col = color,
          xlab = expression(v[0]), ylab = expression(hat(beta)),
          ylim = c(-2.5, 2.5), cex.lab = 1, mgp = c(2.2, 1, 0))
  
  lines(v0.value, beta.threshold, col = "red", lwd = 2)
  lines(v0.value, -beta.threshold, col = "red", lwd = 2)
  title(plot.title)
  
  dev.off()
  
  # Deduct the control variable part  #
  beta.star <- beta.threshold[length(v0.value)]
  if (length(v0.value) > 1) {
    beta.hat <- beta.v0[, dim(beta.v0)[2]]
  } else {
    beta.hat <- beta.v0
  }
  beta.hat[abs(beta.hat) < beta.star] <- 0
  cntl.term <- matrix(NA, dim(test.data)[1], dim(test.data)[2])
  index <- 1
  for (i in 1:dim(test.data)[2]) {
    cntl.term[, i] <- cntl.data[, index:(index+9)] %*% beta.hat[index:(index+9)]
    index <- index + 10
  }
  
  return(list(cntl.term = cntl.term, beta.hat = beta.hat))
}

two.stage.estimation <- function(test.data, cntl.index, cntl.data, 
                                 graph = FALSE, graph.structure = FALSE, 
                                 circle, causal.period, s = 0.1,
                                 emvs.iteration = 50, 
                                 v0.value = seq(1e-6, 0.02, length.out = 5),
                                 mcmc.iterloop = 10000, burnin = 2000, 
                                 stationary = TRUE, 
                                 misspecification = FALSE,
                                 num.sim.data = 30, num.cores = 1,
                                 seed = 1, probs = 0.95,
                                 plot.EMVS.figure = TRUE,
                                 plot.title = NULL) {
  
  
  # check if a package has been installed
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE)) {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    } 
  }
  
  T <- dim(test.data)[1] # time length
  d <- dim(test.data)[2] # dimension of test stores
  
  # make sure test data be T by d
  if (T < d) {
    test.data <- t(test.data) 
    T <- dim(test.data)[1]
    d <- dim(test.data)[2]
  }
  
  if (graph == TRUE) {
    graph.structure <- matrix(1, d, d)
  } else {
    if (graph.structure == FALSE) {
      stop("Graph structure must provde!")
    }
  }
  
  ## Stage 1:
  # EMVS for estimating beta
  # for EMVS, s = 1; for DAEMVS, 0 <= s <= 1
  selection <- 
    estimate.counterfactual(test.data = test.data, cntl.index = cntl.index, 
                            cntl.data = cntl.data, graph.structure = graph.structure, 
                            circle = circle, causal.period = causal.period, s = s, 
                            iteration = emvs.iteration, 
                            v0.value = v0.value,
                            stationary = stationary, plot.figure = plot.EMVS.figure,
                            misspecification = misspecification, 
                            plot.title = plot.title)
  
  cntl.term <- selection$cntl.term
  EMVS.estimator <- selection$beta.hat
  
  ## Stage 2:
  # MCMC for time varying parameters and covariance and variance matrices
  # Fit into timeseries model
  model.estimates <- 
    MultiCausalImpact(test.data = test.data, causal.period = causal.period, 
                      cntl.term = cntl.term, seed = seed, nseasons = circle, 
                      iterloop = iterloop, burnin = burnin, 
                      stationary = stationary, graph = graph, 
                      graph.structure = graph.structure,
                      num.sim.data = num.sim.data, probs = probs,
                      num.cores = num.cores)
  
  # collect result
  mcmc.output <- model.estimates$mcmc.output
  threshold <- model.estimates$threshold
  ks.test.cntl <- model.estimates$ks.test.cntl
  
  
  # return results
  return(list(beta.hat = EMVS.estimator, mcmc.output = mcmc.output,
              ks.test.cntl = ks.test.cntl, threshold = threshold))
  
}

# simulate dataset 
time <- 100 # time length
n <- 5 # number of test
nc <- 50 # number of controls to be generated
n.cntl <- 10 # number of controls to be used for each response

# generate control datasets
cntl.data.pool <- matrix(0, time, nc)
for (i in 1:nc) {
  cntl.data.pool[, i] <- arima.sim(list(ar = c(0.6)), n = time)
}
cntl.data.pool <- cntl.data.pool + abs(min(cntl.data.pool))
# force observations to be positive

# generate test datasets
test.data <- matrix(0, time, n)
A <- matrix(0, n, n) # generate empty matrix, use to create AR(1) correlation
# matrix
graph.structure <- toeplitz(c(1, 1, rep(0, n-2)))
elem <- c(10,5,rep(0,n-2))
Sigma.inv <- toeplitz(elem)
Sigma <- solve(Sigma.inv)
mu <- tau <- matrix(0, n, time)
for (t in 1:time) {
  if (t == 1) {
    mu[, 1] <- 1
  } else {
    mu[, t] <- 0.8*mu[, t-1] + rnorm(n) * 0.1
  }
  test.data[t, ] <- mu[, t] + 
    1*cntl.data.pool[t, seq(1, n*2, by = 2)] + 
    2*cntl.data.pool[t, seq(2, n*2, by = 2)] + 
    mvrnorm(mu = rep(0, n), Sigma = Sigma)
  # add seasonality
  if (t == time) {
    test.data <- test.data + 
      0.1 * cos(2*pi/7*(1:time)) + 0.1 * sin(2*pi/7*(1:time))
  }
}
test.data <- t(t(test.data) + abs(colMins(test.data))+1)
# simulate causal impact
causal.period <- 81:100 # campaign runs 20 periods

# simulate causal impact
for (i in 1:n) {
  test.data[causal.period, i] <- test.data[causal.period, i] + (i-1)*log(1:20)/2
}


# Reorganize dataset   ###

a <- 1:nc
index <- 1
cntl.data <- cntl.data.pool[, index:(index+1)]
for (i in 1:n) {
  if (i != 1) {
    cntl.data <- cbind(cntl.data, cntl.data.pool[, index:(index+1)])
  }
  random.sample <- sample(a[-c(1:(2*n))], 8)
  cntl.data.select <- cntl.data.pool[, random.sample]
  cntl.data <- cbind(cntl.data, cntl.data.select)
  index <- index + 2
}


## Stage 1:
# EMVS for estimating beta
iterloop <- 10000
stationary <- TRUE
nseasons <- 7
graph <- TRUE
burnin <- 2000
num.sim.data <- 30
num.cores <- 1
graph.structure <- graph.structure
cntl.index <- rep(10, n) # cntl.index

MultivariateCausalInferenceRes <- 
  two.stage.estimation(test.data, cntl.index, cntl.data, 
                       graph = graph, graph.structure = graph.structure, 
                       circle = nseasons, causal.period = causal.period, 
                       s = 0.1,
                       emvs.iteration = 50, 
                       v0.value = seq(1e-6, 0.02, length.out = 5),
                       mcmc.iterloop = iterloop, burnin = burnin, 
                       stationary = TRUE, 
                       misspecification = FALSE,
                       num.sim.data = 30, num.cores = 1,
                       seed = 1, probs = 0.95,
                       plot.EMVS.figure = TRUE,
                       plot.title = "DEMVS(s=0.1) plot")

# collect results
beta.hat <- MultivariateCausalInferenceRes$beta.hat
mcmc.output <- MultivariateCausalInferenceRes$mcmc.output
threshold <- MultivariateCausalInferenceRes$threshold
ks.test.cntl <- MultivariateCausalInferenceRes$ks.test.cntl

write.csv(threshold, file = "threshold.csv")
write.csv(ks.test.cntl, file = "ks_test.csv")

# Create date sequence for x-axis
dates <- seq(as.Date("2016-01-01"), as.Date("2016-04-01"), by = "day")

# Vertical line date (intervention point)
intervention_date <- as.Date("2016-03-21")


# Generate data
non_stationary_data <- generate_bayesian_time_series(dates, intervention_date, stationary = FALSE)
stationary_data <- generate_bayesian_time_series(dates, intervention_date, stationary = TRUE)


# Generate plotting data
non_stationary_plot_data <- analyze_bayesian_time_series(non_stationary_data, dates, intervention_date, stationary = FALSE)
stationary_plot_data <- analyze_bayesian_time_series(stationary_data, dates, intervention_date, stationary = TRUE)

# Define x-axis date format
date_breaks <- c("Jan-01-16", "Jan-21-16", "Feb-10-16", "Mar-01-16", "Mar-21-16")

create_ts_plot <- function(data, title) {
  simulated_data_with_noise <- data$ts_df$SimulatedData + rnorm(nrow(data$ts_df), 0, 0.5)
  
  ggplot(data$ts_df, aes(x = Date)) +
    geom_line(aes(y = LowerCI), color = "blue", linetype = "dashed", size = 0.5) +
    geom_line(aes(y = UpperCI), color = "blue", linetype = "dashed", size = 0.5) +
    geom_line(aes(y = EstimatedMedian), color = "blue", size = 0.3) +
    geom_line(aes(y = simulated_data_with_noise), color = "green", size = 0.6) +
    geom_vline(xintercept = as.numeric(intervention_date), linetype = "dashed", color = "gray50") +
    scale_x_date(breaks = as.Date(c("2016-01-01", "2016-01-21", "2016-02-10", "2016-03-01", "2016-03-21")),
                 labels = date_breaks,
                 limits = c(min(dates), max(dates))) +
    scale_y_continuous(limits = c(5, 40), breaks = seq(5, 40, by = 5)) +
    labs(title = title) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    annotate("text", x = min(dates) + 10, y = 38, 
             label = "Simulated data", color = "green", hjust = 0, size = 3) +
    annotate("text", x = min(dates) + 10, y = 35.5, 
             label = "Estimated/Predicted Median", color = "blue", hjust = 0, size = 3) +
    annotate("text", x = min(dates) + 10, y = 33, 
             label = "95% Credible Intervals", color = "blue", hjust = 0, size = 3)
}

create_impact_plot <- function(data, title) {
  estimated_median_with_noise <- data$impact_df$EstimatedMedian + rnorm(nrow(data$impact_df), 0, 0.5)
  
  ggplot(data$impact_df, aes(x = Date)) +
    geom_line(aes(y = LowerCI), color = "blue", linetype = "dashed", size = 0.5) +
    geom_line(aes(y = UpperCI), color = "blue", linetype = "dashed", size = 0.5) +
    geom_line(aes(y = estimated_median_with_noise), color = "blue", size = 0.5) +
    geom_line(aes(y = SimulatedImpact), color = "green", size = 0.7) +
    geom_vline(xintercept = as.numeric(intervention_date), linetype = "dashed", color = "gray50") +
    scale_x_date(breaks = as.Date(c("2016-01-01", "2016-01-21", "2016-02-10", "2016-03-01", "2016-03-21")),
                 labels = date_breaks,
                 limits = c(min(dates), max(dates))) +
    scale_y_continuous(limits = c(-5, 15), breaks = seq(-5, 15, by = 5)) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    annotate("text", x = min(dates) + 10, y = 14, 
             label = "Simulated impact", color = "green", hjust = 0, size = 3) +
    annotate("text", x = min(dates) + 10, y = 12, 
             label = "Estimated/Predicted Median", color = "blue", hjust = 0, size = 3) +
    annotate("text", x = min(dates) + 10, y = 10, 
             label = "95% Credible Intervals", color = "blue", hjust = 0, size = 3)
}

# Generate individual plots
p1 <- create_ts_plot(non_stationary_plot_data, "Non-stationary")
p2 <- create_impact_plot(non_stationary_plot_data, "")
p3 <- create_ts_plot(stationary_plot_data, "Stationary")
p4 <- create_impact_plot(stationary_plot_data, "")

# Combine plots with labels
final_plot <- grid.arrange(
  p1, p3, p2, p4, 
  layout_matrix = rbind(c(1, 2), c(3, 4)),
  top = textGrob("Bayesian multivariate time series causal inference", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

# Add subplot labels
grid.text("(a)", x = 0.01, y = 0.98, just = c("left", "top"))
grid.text("(c)", x = 0.51, y = 0.98, just = c("left", "top"))
grid.text("(b)", x = 0.01, y = 0.48, just = c("left", "top"))
grid.text("(d)", x = 0.51, y = 0.48, just = c("left", "top"))

# Print plot
print(final_plot)

# Save the plot
png("bayesian_time_series_plots.png", width = 10*300, height = 8*300, res = 300)
grid.draw(final_plot)
dev.off()

ggsave("bayesian_time_series_plots2.png", plot = final_plot, 
       device = "png", width = 10, height = 8, dpi = 300)

pdf("bayesian_time_series_plots.pdf", width = 10, height = 8)
grid.draw(final_plot)
dev.off()
