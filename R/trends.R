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
  t_m <- stats::lm(y~0+d_t)
  # Stochastic variance of growth rates
  stoc_v <- stats::sigma(t_m)^2
  # Return mean growth rate, variance and confidence interval
  res <- c(stats::coef(t_m), stoc_v, stats::confint(t_m))
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
  t_m <- stats::lm(d_n~d_t)
  # Return slope and confidence interval
  res <- c(stats::coef(t_m)[2], stats::confint(t_m)[2,])
  names(res) <- c("trend", "l95", "u95")
  
  return(res)
}


#' Title
#'
#' @param x A data.frame. A community matrix of species abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' @param method Character. Method to estimate the trends, one of "dennis" or "loglinear". Default "dennis". 
#'
#' @returns A data.frame with one row per species in the community.
#' 
#' @export
comm_trend <- function(x, time_col = "time", method = "dennis", plot = FALSE){
  # Match variance function
  trend_func <- switch(
    method,
    dennis = trend_dennis,
    loglinear = trend_loglinear
  )
  
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = "var", rm = TRUE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  trends <- as.data.frame(
    cbind(
      taxa = colnames(x), 
      as.data.frame(
        do.call("rbind", 
                apply(x, MARGIN = 2, trend_func, simplify = F)
        )
      )
    )
  )
  rownames(trends) <- NULL
  
  # Plot abundances
  if (plot) {
    par(mfrow = c(1,2))
    plot_com(x)
    # for (i in seq_len(nrow(trends))) {
    #   graphics::abline(a = mean(log(x[,i])), 
    #                    b = trends$trend[i],
    #                    col = i)
    # }
    
    plot(x = trends[1,]$trend, y = seq_along(trends$taxa)[1],
         xlim = c(min(trends$l95), max(trends$u95)),
         ylim = c(min(seq_along(trends$taxa)), max(seq_along(trends$taxa))),
         pch = 19,
         col = 1,
         xlab = "trend (log)",
         ylab = "taxa", 
         yaxt = "n")
    arrows(x0 = trends$l95, x1 = trends$u95, y0 = seq_along(trends$taxa),
           code = 3, length = 0.05, angle = 90)
    for (i in 2:nrow(trends)) {
      graphics::points(x = trends[i,]$trend, y = seq_along(trends$taxa)[i], 
                       pch = 19,
                       col = i)
    }
    
    axis(2, at = seq_along(trends$taxa), labels = trends$taxa, las = 2)
    abline(v = 0, lty = "dashed")
    
    par(mfrow = c(1,1), 
        xpd=FALSE, 
        mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  return(trends)
}