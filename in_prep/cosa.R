simmvnrom <- function(n_sp = 10,
                      years = 25,
                      tot_abu = 200 * n_sp,
                      power = 1.6, # values closer to 1 indicate more stable dominants and thus more dominance effect
                      bound_pos = TRUE,
                      corr = 0.5, # correlation between time series
                      p = 0.8, # Higher values mean less dominance
                      trend = FALSE,
                      trend_sd = 0.01)  {
  
  # Vector of mean abundances. First get a vector of relative abundances 
  # that add up to 1 using the dirchlet distribtion. Parameter alpha (p)
  # controls the spread of the values, with higher values leading to 
  # more even relative abundances and thus less dominance
  mean_abu <- sort(tot_abu * gtools::rdirichlet(1, alpha = rep(p, n_sp))[1,], decreasing = TRUE)
  
  # Get Sd of abundances from TPL
  sd_abu <- sqrt(mean_abu ** power)
  sd_abu = sd_abu / max(sd_abu)
  
  # Variance-covariance matrix for MVN distribution
  k <- 3
  # This makes positive definitive matrix (necessary for MVN) 
  # with correlation between species
  Lambda <- matrix(rnorm(n_sp * k, sd = sqrt(corr)), n_sp, k)
  Psi <- diag(sd_abu^2 * (1 - corr)) 
  ss <- Lambda %*% t(Lambda) + Psi
  
  # Simulate random variation around mean cover for each species
  # drawn from multivariate normal so species can correlate
  simcom <- matrix(0, years, n_sp)
  simcom[1,] <- mean_abu
  
  for (t in 2:years) {
    r <- MASS::mvrnorm(n = 1, mu = rep(0, n_sp), Sigma = ss)
    simcom[t,] <- simcom[t-1,] * exp(r)
    
    # Force positive values
    if (bound_pos) {
      abi[abi < 0] <- 0
    }
  }
  plot_com(simcom)
  res <- list(sim_data = as.data.frame(simcom),
              av_trend = colMeans(apply(log(simcom),1,diff)))
  return(res)
} 