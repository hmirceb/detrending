#' Simulate community from a multivariate normal
#'
#' @param n_sp Numeric. Number of species in the community.
#' @param years Numeric. NUmber of years to simulate.
#' @param tot_abu Numeric. Total abundance of the community.
#' @param power Numeric. Exponent of the Taylor's Power Law to estimate variance from mean abundance.
#' @param bound_pos Boolean. Bound abundance values to be positive. Default TRUE.
#' @param corr Numeric. Correlation between populations.
#' @param p Numeric. Dispersion value for a Dirichlet distribution bound between 0-1. Controls the evenness of the community with lower values indicating less dominance.
#' @param trend Boolean. Apply a trend to the data. Default FALSE.
#' @param trend_mean Numeric. Mean of the trend.
#' @param trend_sd Numeric. Standard deviation of the trend.
#' @param bimodal_trend 
#'
#' @returns A named list containing a data.frame with the simulated data, the true mean trends of each simulated species and the parameters used to simulate the data.
#' 
#' @export
sim_mvcomm <- function(n_sp = 10,
                    years = 25,
                    tot_abu = 200 * n_sp,
                    power = 1.6,
                    bound_pos = TRUE,
                    corr = 0.5,
                    p = 0.8,
                    switch_trend = FALSE,
                    trend_mean = 1,
                    trend_sd = 0.01,
                    bimodal_trend = FALSE){
  
  # Vector of mean abundances. First get a vector of relative abundances 
  # that add up to 1 using the dirchlet distribtion. Parameter alpha (p)
  # controls the spread of the values, with higher values leading to 
  # more even relative abundances and thus less dominance
  mean_abu <- sort(tot_abu * gtools::rdirichlet(1, alpha = rep(p, n_sp))[1,], decreasing = TRUE)
  
  # Create a matrix of abundances
  abu_matrix <- matrix(rep(mean_abu, times = years),
                       nrow = years, ncol = n_sp, byrow = TRUE)
  
  trend <- seq(-1, 1, length.out = years)
  trend_resp <- response(state = switch_trend,
                         n_sp = n_sp,
                         mean = trend_mean,
                         sd = trend_sd,
                         bimodal = bimodal_trend,
                         comp = FALSE)
  
  # Get Sd of abundances from TPL
  sd_abu <- sqrt(mean_abu ** power)
  
  # Variance-covariance matrix for MVN distribution
  k <- n_sp-1
  # This makes positive definitive matrix (necessary for MVN) 
  # with correlation between species
  Lambda <- matrix(rnorm(n_sp * k, sd = 1), n_sp, k)
  Psi <- diag(sd_abu^2 * (1 - corr)) 
  ss <- Lambda %*% t(Lambda) + Psi
  
  # Simulate random variation around mean cover for each species
  # drawn from multivariate normal so species can correlate
  simcom <- matrix(0, years, n_sp)
  for (j in 1:years) {
    abi <- MASS::mvrnorm(n = 1, 
                         mu = abu_matrix[j,] * (1 + trend[j] * trend_resp), 
                         Sigma = ss)
    # Force positive values if necessary
    if (bound_pos) {
      abi[abi < 0] <- 0
    }
    simcom[j, ] <- abi
  }
  
  # Add a small offset (1% of the mean abundance of each species) to avoid having 0s
  offset <- colMeans(simcom)*0.01 
  # Add vector to matrix rowwise
  simcom <- sweep(x = simcom, MARGIN = 2, STATS = offset, FUN = "+")
  
  # Results into list
  res <- list(sim_data = as.data.frame(simcom),
              true_trend = exp(colMeans(apply(log(simcom), 2, diff))),
              params = c(n_sp = n_sp,
                         years = years,
                         tot_abu = tot_abu,
                         power = power,
                         bound_pos = bound_pos,
                         corr = corr,
                         p = p,
                         switch_trend = switch_trend,
                         trend_mean = trend_mean,
                         trend_sd = trend_sd,
                         bimodal_trend = bimodal_trend))
  
  return(res)
}


#' Simulate fluctuations in species abundances across time
#'
#' The arguments in \code{syngenr} are set to predefined default values, which can of course
#' be customised, or sampled from parameter spaces.
#'
#' @param years the length of the timeseries in years.
#' @param n_sp number of species in the community.
#' @param tot_abu total abundance of the community, representing e.g. the number of
#' individuals or the total amount of biomass.
#' @param power the slope of the relationship between log(mean) and log(variance) of
#' the abundances of the species.
#' @param switch_env character vector of length 1, either "on" or "off". Defines if the
#' species abundances repond to a hypothetical environmental cue.
#' @param mean_env_resp the mean of the normal distribution from which each of the
#' species responses to the environemtnal cue is drawn.
#' @param sd_env_resp the standard deviation around the mean of the normal distribution
#' from which each of the species responses to the environemtnal cue is drawn.
#' @param bimodal_env Logical. Making the reponse to the environmental cue either
#' uniform among the species (if FALSE), or making the majority of half of the
#' species respond positively, and the other half negatively.
#' @param comp Logical. If TRUE, species exhibit compensatory dynamics, i.e. the gain
#' of abundance in a species from one year to the next, is compensated by the
#' loss of abundance in another species, where the latter has a similar mean
#' abundance value.
#' @param switch_trend Either "on" or "off", defining if there is a general monotonic
#' trend in abundances of species across the timeseries.
#' @param mean_trend_resp the mean of the normal distribution from which each of the
#' species responses to the longterm trend is drawn. Default is 1.
#' @param sd_trend_resp the standard deviation around the mean of the normal
#' distribution from which each of the species responses to the longterm trend is
#' drawn. Default is 1.
#' @param bimodal_trend Logical. If TRUE, the majority of half of the species exhibit a
#' positive long term trend of abundances, and the other half expresses a
#' negative long term trend. If FALSE, most of the species exhibit a uniform long
#' term trend in abundances, depending on the value set for mean_trend and
#' mean_sd.
#' @param bound_pos Logical. If true, abundance values that are simulated to be
#' negative, will be set to zero.
#'
#' @return The output is a list with four elements
#'
#' time_species_matrix: the simulated temporal community data where species are
#' columns and years are rows.
#'
#' param_years: values for the evnironmental cue and the trend throughout the
#' years. Note that these contain values even if the environment or trend are
#' switched off.
#'
#' param_species: a data frame containing the responses to the environment and
#' the long term trend, as well as the mean abundance and its standard deviation
#' for each species in the community.
#'
#' param_general: a data frame with only one row, containing all the parameter
#' settings from the function call.
#' @export syngenr
syngenr <- function(years = 100,
                    n_sp = 16,
                    max_rel_abu = 0.6,
                    tot_abu = 300,
                    power = 1.8,
                    switch_env = FALSE,
                    mean_env_resp = 1,
                    sd_env_resp = 1,
                    bimodal_env = FALSE,
                    comp = FALSE,
                    switch_trend = FALSE,
                    mean_trend_resp = 1,
                    sd_trend_resp = 1,
                    bimodal_trend = FALSE,
                    bound_pos = TRUE) {
  
  mean_abu <- tot_abu * geom_seq(max_rel_abu, n_sp)
  sd_abu <- sqrt(mean_abu ** power)
  env <- sample(seq(-0.5, 0.5, length.out = 3), years, replace = T)
  env_resp <- response(
      state = switch_env,
      n_sp = n_sp,
      mean = mean_env_resp,
      sd = sd_env_resp,
      bimodal = bimodal_env,
      comp = comp
    )
  trend = seq(-1, 1, length.out = years)
  trend_resp <- response(state = switch_trend,
                         n_sp = n_sp,
                         mean = mean_trend_resp,
                         sd = sd_trend_resp,
                         bimodal = bimodal_trend,
                         comp = FALSE
                         )
  
  simcom <- matrix(0, years, n_sp)
  er <- env_resp
  tr <- trend_resp
  for (i in 1:n_sp) {
    abi <- vector("numeric", years)
    for (j in 1:years) {
      abi[j] <-
        rnorm(1, mean_abu[i] * (1 + env[j] * er[i]) * (1 + trend[j] * tr[i]),
              sd_abu[i])
      if (bound_pos)
        abi[abi < 0] <- 0
    }
    simcom[, i] <- abi
  }
  param_species <-
    data.frame(env_resp, trend_resp, mean_abu, sd_abu)
  param_years <- data.frame(env, trend)
  param_general <- data.frame(
    n_sp,
    max_rel_abu,
    tot_abu,
    power,
    switch_env,
    bimodal_env,
    mean_env_resp,
    sd_env_resp,
    switch_trend,
    bimodal_trend,
    mean_trend_resp,
    sd_trend_resp,
    comp,
    bound_pos
  )
  res <- list(
    time_species_matrix = simcom,
    param_years = param_years,
    param_species = param_species,
    param_general = param_general
  )
  return(res)
}

#' Function to generate a geometric distribution of relative abundances
#'
#' @param max_rel_abu a single numeric value that represents the realtive
#' abundance of the most abundant species in the community.
#' @param n_sp number of species in the community.
#'
#' @return A vector with the length of the number of species containing the
#' relative abundances of the species.
#' @export geom_seq
geom_seq <- function(max_rel_abu, n_sp) {
  rel.abu <- max_rel_abu
  remaining <- 1 - rel.abu
  for (i in 2:n_sp) {
    rel.abu[i] <- max_rel_abu * remaining
    remaining <- 1 - sum(rel.abu)
  }
  geom <- rel.abu / sum(rel.abu)
  return(geom)
}


#' Function to simulate response of species abundances to environmental cues
#'
#' \code{response} simulates the strength and direction with which the change in
#' species abundance across time responds to a stochastic environmental cue, or,
#' a general directional trend in abundance. This is only a helper function that
#' is used internally in the main function to simulate temporal communities,
#' syngenr.
#'
#' @param state Boolean. Defines if the species have some kind of trend. Default TRUE.
#' @param bimodal Boolean. If TRUE, half of the species respond opposite to the other
#' half, if FALSE, all species respond in the same direction. How many species
#' respond in a given direction also depends on the settings of mean and sd.
#' @param mean Numeric. The mean of the normal distribution from which each of the species
#' responses is drawn.
#' @param sd Numeric. The standard deviation around the mean of the normal distribution from
#' which each of the species responses is drawn.
#' @param n_sp Numeric. Number of species in the community.
#' @param comp Boolean. If TRUE, species will exhibit compensatory dynamics. This is
#' simulated by having species of similar abundance to respond in opposite
#' directions. This argument is therefore only meaning full when having set
#' bimodal = TRUE.
#'
#' @return The output is a vector with the length of the number of species, each
#' value representing the response of a species to hypotehtical cue (i.e. in the
#' context of the simulation, either an environmental signal or a longterm
#' monotonic abundance trend).
#' @export response
response <- function(state = TRUE,
                     bimodal = FALSE,
                     mean = 1,
                     sd = 1,
                     n_sp,
                     comp = FALSE) {
  if ( !isTRUE(state) ) {
    # species do not respond
    resp <- rep(0, n_sp)
  } else{
    if (bimodal) {
      # species split into negative and positive responders
      n_half <- ceiling(n_sp / 2)
      n_otherhalf <- n_sp - n_half
      resp_pos <- rnorm(n_half, abs(mean), sd)
      resp_neg <- rnorm(n_otherhalf,-abs(mean), sd)
      
      if (comp) {
        # model compensatory dynamics
        resp <- vector("numeric", n_sp)
        resp[!is_even(1:n_sp)] <- resp_pos
        resp[is_even(1:n_sp)] <- resp_neg
      } else{
        index <- sample(1:n_sp, n_half)
        resp <- vector("numeric", n_sp)
        resp[index] <- resp_pos
        resp[-index] <- resp_neg
      }
    } else{
      # majority of species respond either positive or negative
      resp <- rnorm(n_sp, mean, sd)
    }
  }
  return(resp)
}