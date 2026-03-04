#' Estimate mean annual growth rate following Dennis et al. (2001). 
#'
#' @param x Numeric. A vector of abundances.
#' @param time Numeric. A vector with the time steps corresponding to each value in x.
#'
#' @returns A named list with the mean annual growth rate in the natural logarithm scale and its variance.
#' 
#' @references
#' - Dennis, B., Munholland, P. L., & Scott, J. M. (1991). Estimation of Growth and Extinction Parameters for Endangered Species. Ecological Monographs, 61(2), 115–143.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @export
trend_dennis  <- function(x, time = NULL){
  # Make time variable if missing
  if ( is.null(time) ) {
    time <- seq_along(x)
  }
  # Get good years
  g_y <- ( x >= 0 | !is.na(x) )
  # Keep only good data
  x <- x[g_y]
  # Keep only good times
  time <- time[g_y]
  # Square root of difference between time steps
  d_t <- sqrt(diff(time))
  # Log ratios between time steps
  d_n <- diff(log(x))
  # Average log ratios
  y <- d_n / d_t
  # Lm for log ratios and corrected time
  t_m <- lm(y~0+d_t)
  # Stochastic variance of growth rates
  stoc_v <- sigma(t_m)^2
  # Return mean growth rate, variance and confidence interval
  res <- c(coef(t_m), stoc_v, confint(t_m))
  names(res) <- c("trend", "var", "l95", "u95")
  return(res)
}

#' Estimate mean annual growth rate using linear regression 
#'
#' @param x Numeric. A vector of abundances.
#' @param time Numeric. A vector with the time steps corresponding to each value in x.
#'
#' @returns A numeric value with the slope of a linear regression of the log transformed abundances.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @export
trend_loglinear  <- function(x, time = NULL){
  # Make time variable if missing
  if ( is.null(time) ) {
    time <- seq_len(length(x))
  }
  # Get good years
  g_y <- ( x >= 0 | !is.na(x) )
  # Keep only good years
  x <- x[g_y]
  # Keep only good years
  time <- time[g_y]
  # Set time as increasing from 0
  d_t <- time[g_y]-min(time[g_y])
  # log transform
  d_n <- log(x)
  # linear regression
  t_m <- lm(d_n~d_t)
  # Return slope and confidence interval
  res <- c(coef(t_m)[2], confint(t_m)[2,])
  names(res) <- c("trend", "l95", "u95")
  
  return(res)
}


