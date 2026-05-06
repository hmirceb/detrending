#' Estimate mean annual growth rate following Dennis et al. (2001). 
#'
#' This function estimates the mean annual growth rate of a time series of abundance data following Dennis et al. (2001) (see Details).
#'
#' @param x Numeric. A vector of abundances.
#' @param time Numeric. A vector with the time steps corresponding to each value in x. If not provided data are assumed to be in order.
#'
#' @details For a given a time series of abundances \eqn{n_{t}} the function estimates its mean growth rate \eqn{\mu} and variance \eqn{\sigma^2} using a linear regression model without intercept such as:
#' \deqn{y_{i} \sim \mu t_{i} + \epsilon_{i} }
#' \deqn{\epsilon_{i} \sim Normal(0, \sigma^2)}
#'  where \eqn{ y_{i} = \dfrac{ \ln{ ( n_{t}/n_{t-1}) } }{ t_{i} } } and \eqn{ t_{i} = \sqrt{ t_{t}-t_{t-1} } }
#'  
#'  Note that the confidence intervals for the estimate of \eqn{\mu} are not reliable.
#'  
#' @returns A named vector with the mean annual growth rate of a population in the natural logarithm scale, its confidence interval and p-value.
#' 
#' @references
#' - Dennis, B., Munholland, P. L., & Scott, J. M. (1991). Estimation of Growth and Extinction Parameters for Endangered Species. Ecological Monographs, 61(2), 115–143.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate time series
#' ts <- 5^seq(1, 2, by = 0.01)
#' mean(diff(log(ts))) # True trend = 0.016 (~1.62%)
#' # Simulate some random noise
#' noise <- rnorm(length(ts))
#' 
#' # Estimate trend
#' trend_dennis(ts+noise)
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
  res <- c(stats::anova(t_m)$F[1], stats::coef(t_m), stats::confint(t_m), stats::anova(t_m)$`Pr(>F)`[1])
  names(res) <- c("F", "trend", "l95", "u95", "p")
  
  return(res)
}

#' Estimate mean annual growth rate using linear regression 
#'
#' @param x Numeric. A vector of abundances.
#' @param time Numeric. A vector with the time steps corresponding to each value in x. If not provided data are assumed to be in order.
#'
#' @details For a given a time series of abundances \eqn{n_{t}} the function estimates a linear regression model with log-transformed abundances as the response variable and time steps (with \eqn{t_{0} = 0, t_{1} = 1...}) as the explanatory variable.
#'
#' @returns A named vector with the slope of a linear regression of the log transformed abundances, its confidence interval and p-value.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate time series
#' ts <- 5^seq(1, 2, by = 0.01)
#' mean(diff(log(ts))) # True trend = 0.016 (~1.62%)
#' # Simulate some random noise
#' noise <- rnorm(length(ts))
#' 
#' # Estimate trend
#' trend_loglinear(ts+noise)
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
  res <- c(stats::anova(t_m)$F[1], stats::coef(t_m)[2], stats::confint(t_m)[2,], stats::anova(t_m)$`Pr(>F)`[1])
  names(res) <- c("F", "trend", "l95", "u95", "p")
  
  return(res)
}

#' Perform Redundancy Analysis on species abundances with time as explanatory variable to assess community level trends
#'
#' @param x A data.frame. A community matrix of species abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' @param community_col Character. Name of column with the community identifier.
#' @param scale Boolean. Scale abundances to have mean 0 and standard deviation 1. Default TRUE.
#' @param perm Numeric. Number of permutations for significance testing. Default 999.
#' @details This function estimates temporal trends in abundance at the community level by conducting Redundanncy Analysis (RDA) on species abundances with time as an explanaroty variable.
#'
#' @returns An object of class `mv_trend`, a named list with three elements:
#'  - `anova`: A data.frame with the *F-value* of the Redundancy Analysis and its p-value based on the number of permutations chosen.
#'  
#'  - `rda`: An object of class `rda` with the scores of the RDA.
#'  
#'  - `time`: A vector with the time variable used in the RDA.
#'
#' @references
#'  - Legendre, P. and Legendre, L. (2012) Numerical Ecology. 3rd English ed. Elsevier.
#'  
#' @author Jan Lepš, \email{suspa@@prf.jcu.cz} 
#' @author Aleš Lisner, \email{lisnea00@@jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#'
#' @examples
#' require(detrending)
#' 
#' # Simulate community data with trends
#' comm_df <- sim_mvcomm(trend_mean = 0.3, bimodal_trend = TRUE)
#' 
#' # Estimate multivariate community trend
#' trend_mv(comm_df$sim_data, time_col = "time")
#' 
#' @export
trend_mv <- function(x, time_col = "time", community_col = "comm", scale = TRUE, perm = 999){
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = "var", rm = FALSE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # prepare data matrices for RDA
  comm <- x[!colnames(x) %in% c(community_col, time_col)]
  t <- x[colnames(x) %in% time_col]
  
  # make formula for RDA
  rda_formula <- stats::as.formula(paste("comm", "~", time_col))
  
  # RDA
  comm_rda <- vegan::rda(formula = rda_formula, data = t, scale = scale)
  rda_sign <- stats::anova(comm_rda, permutations = perm)
  # format anova output
  res <- c(rda_sign$F[1], NA, NA, NA, rda_sign$`Pr(>F)`[1])
  names(res) <- c("F", "trend", "l95", "u95", "p")
  # make mv_trend object
  results <- list(anova = res, rda = comm_rda, time = t[,1])
  class(results) <- "mv_trend"
  return(results)

}

#' Estimate population trends in a community
#'
#' @param x A data.frame. A community matrix of species abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' @param community_col Character. Name of column with the community identifier.
#' @param method Character. Method to estimate the trends, one of "dennis", "loglinear" or "rda". Default "dennis".
#' @param plot Boolean. Plot species abundances and their estimated trends. Default FALSE. 
#' @param title Character. Title for the plot. Default NULL.
#' 
#' @returns A data.frame with the trend (in the natural logarithm scale) for each species in the community along with its variance and 95% confidence interval.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data with trends
#' comm_df <- sim_mvcomm(trend_mean = 0.3, bimodal_trend = TRUE)
#' 
#' # Estimate trend for each species and plot them
#' comm_trend(comm_df$sim_data, method = "loglinear", plot = TRUE)
#' @export
comm_trend <- function(x, time_col = "time", community_col = "comm", method = "loglinear", plot = FALSE, title = NULL){
  # Match trend estimation function
  method_matched <- match.arg(method, choices = c("dennis", "loglinear", "rda"))
  
  trend_func <- switch(
    method_matched,
    dennis = trend_dennis,
    loglinear = trend_loglinear,
    rda = trend_mv
  )
  
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = "var", rm = TRUE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # remove community column
  x <- x[,!colnames(x) %in% c(time_col, community_col)]
  
  z <- apply(x, 2, mean)*0.01
  x <- sweep(x, 2, STATS = z, FUN = "+")
  
  # Estimate trends using apropriate function (this is like this because rda has different output format)
  if (method_matched == "rda") {
    trends <- trend_func(x = x, time_col = time_col, community_col = community_col)
    linear_trends <- as.data.frame(
      cbind(
        taxon = colnames(x), 
        as.data.frame(
          do.call("rbind", 
                  apply(x, MARGIN = 2, trend_loglinear, simplify = F)
          )
        )
      )
    )
    rownames(linear_trends) <- NULL
  }
  
  if (method_matched %in% c("dennis", "loglinear")) {
    trends <- as.data.frame(
      cbind(
        taxon = colnames(x), 
        as.data.frame(
          do.call("rbind", 
                  apply(x, MARGIN = 2, trend_func, simplify = F)
          )
        )
      )
    )
    rownames(trends) <- NULL
  }
  
  # Plot abundances
  if ( isTRUE(plot) ) {
    
    # plot mv trends
    if (method_matched == "rda") {
      # panel layout
      graphics::layout(matrix(c(1,1,2,2), nrow = 1, ncol = 4, byrow = TRUE))
      # plot community
      # plot_com(x, title = title)
      # plot rda and trends
      plot(trends)
      
      # plot species trends
      plot(x = linear_trends[1,]$trend, y = seq_along(linear_trends$taxon)[1],
           xlim = c(min(linear_trends$l95), max(linear_trends$u95)),
           ylim = c(min(seq_along(linear_trends$taxon)), max(seq_along(linear_trends$taxon))),
           pch = 19,
           col = 1,
           xlab = "Trend (log)",
           ylab = "", 
           yaxt = "n")
      graphics::arrows(x0 = linear_trends$l95, x1 = linear_trends$u95, y0 = seq_along(linear_trends$taxon),
                       code = 3, length = 0.05, angle = 90)
      for (i in 2:nrow(linear_trends)) {
        graphics::points(x = linear_trends[i,]$trend, y = seq_along(linear_trends$taxon)[i], 
                         pch = 19,
                         col = i)
      }
      # shorten taxa names
      labs <- short_names(linear_trends$taxon)
      graphics::axis(2, at = seq_along(linear_trends$taxon), 
                     labels = labs, 
                     las = 2)
      graphics::abline(v = 0, lty = "dashed") 
    }
    # plot log trends
    if (method_matched %in% c("dennis", "loglinear")) {
      # setup two panels
      graphics::layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
      # plot community
      plot_com(x, title = title)
      # plot species trends
      plot(x = trends[1,]$trend, y = seq_along(trends$taxon)[1],
           xlim = c(min(trends$l95), max(trends$u95)),
           ylim = c(min(seq_along(trends$taxon)), max(seq_along(trends$taxon))),
           pch = 19,
           col = 1,
           xlab = "Trend (log)",
           ylab = "", 
           yaxt = "n")
      graphics::arrows(x0 = trends$l95, x1 = trends$u95, y0 = seq_along(trends$taxon),
                       code = 3, length = 0.05, angle = 90)
      for (i in 2:nrow(trends)) {
        graphics::points(x = trends[i,]$trend, y = seq_along(trends$taxon)[i], 
                         pch = 19,
                         col = i)
      }
      # shorten taxa names
      labs <- short_names(trends$taxon)
      graphics::axis(2, at = seq_along(trends$taxon), 
                     labels = labs, 
                     las = 2)
      graphics::abline(v = 0, lty = "dashed") 
    }
    # reset graphics
    graphics::par(mfrow = c(1,1), 
                  xpd = FALSE, 
                  mar = c(5.1, 4.1, 4.1, 2.1))
  }
  
  if (method_matched == "rda") {
    return(trends$anova)
  }
  if (method_matched %in% c("dennis", "loglinear")) {
    return(trends)
  }
}

# #' @export
# plot.mv_trend <- function(x, ...) {
#   # get community scores 
#   rda_sites <- vegan::scores(x$rda)$sites
#   # get species scores
#   rda_species <- vegan::scores(x$rda)$species
#   # empty plot
#   plot(x = rda_sites[,1], y = rda_sites[,2],
#        type = "n", 
#        xlab = "Time (RDA1)", ylab = "Community composition (PC1)",
#        xlim = 2 * c(range(rda_sites[,1]))) # adjust limits
#   # add reference lines
#   graphics::abline(h = 0, lty = 3)
#   graphics::abline(v = 0, lty = 3)
#   # arrows from timestep to timestep
#   for (i in 2:nrow(rda_sites)) {
#     # get start and end of each arrow
#     x0 <- rda_sites[i-1,1]
#     y0 <- rda_sites[i-1,2]
#     x1 <- rda_sites[i,1]
#     y1 <- rda_sites[i,2]
#     # only the last one is a proper arrow
#     angle <- ifelse(i == max(nrow(rda_sites)), 20, 0)
#     graphics::arrows(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
#                      angle = angle, length = 0.1)
#   }
#   # add time steps with text and colored points
#   graphics::text(x = rda_sites[,1],
#                  y = rda_sites[,2], 
#                  labels = x$time, 
#                  pos = 3)
#   graphics::points(x = rda_sites[,1],
#                    y = rda_sites[,2], 
#                    pch = 21, 
#                    bg = "grey")
#   # add arrow indicating time direction
#   graphics::arrows(x0 = 0, 
#                    y0 = 0, 
#                    y1 = 0, 
#                    x1 = 1.8 * max(rda_sites[,1]), 
#                    col = "blue", 
#                    length = 0.1)
#   graphics::text(y = 0, 
#                  x = 1.8 * max(rda_sites[,1]), 
#                  labels = "time", 
#                  pos = 3, 
#                  col = "blue")
#   # F statistic and p value from permutation test
#   graphics::mtext(side = 3, 
#                   adj = 1, 
#                   text = paste0("F=", round(x$anova["F"], 3), ", p=", x$anova["p"]))
#   
#   # plot with species trends
#   plot(x = rda_species[,1], 
#        y = seq_along(rownames(rda_species)), 
#        type = "n",
#        xlab = "Trend (RDA1)",
#        ylab = "",
#        yaxt = "n",
#        xlim = 1.1 * range(rda_species[,1]))
#   graphics::abline(v = 0, 
#                    lty = "dashed")
#   # shorten species names
#   labs <- short_names(rownames(rda_species))
#   
#   graphics::axis(2, at = seq_along(rownames(rda_species)), 
#                  labels = labs, 
#                  las = 2)
#   graphics::points(x = rda_species[,1], 
#                    y = seq_along(rownames(rda_species)), 
#                    col = seq_along(rownames(rda_species)), 
#                    pch = 19,
#                    cex = 1.5)
# }

#' @export
plot.mv_trend <- function(x, ...) {
  # get community scores
  rda_sites <- vegan::scores(x$rda)$sites
  # get species scores
  rda_species <- vegan::scores(x$rda)$species
  # empty plot
  plot(x = rda_sites[,1], y = rda_sites[,2],
       type = "n",
       xlab = "Time (RDA1)", ylab = "Community composition (PC1)",
       xlim = 2 * c(range(rda_sites[,1]))) # adjust limits
  # add reference lines
  graphics::abline(h = 0, lty = 3)
  graphics::abline(v = 0, lty = 3)
  # arrows from timestep to timestep
  for (i in 2:nrow(rda_sites)) {
    # get start and end of each arrow
    x0 <- rda_sites[i-1,1]
    y0 <- rda_sites[i-1,2]
    x1 <- rda_sites[i,1]
    y1 <- rda_sites[i,2]
    # only the last one is a proper arrow
    angle <- ifelse(i == max(nrow(rda_sites)), 20, 0)
    graphics::arrows(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
                     angle = angle, length = 0.1)
  }
  # add time steps with text and colored points
  graphics::text(x = rda_sites[,1],
                 y = rda_sites[,2],
                 labels = x$time,
                 pos = 3)
  graphics::points(x = rda_sites[,1],
                   y = rda_sites[,2],
                   pch = 21,
                   bg = "grey")
  # add arrow indicating time direction
  graphics::arrows(x0 = 0,
                   y0 = 0,
                   y1 = 0,
                   x1 = 1.8 * max(rda_sites[,1]),
                   col = "blue",
                   length = 0.1)
  graphics::text(y = 0,
                 x = 1.8 * max(rda_sites[,1]),
                 labels = "time",
                 pos = 3,
                 col = "blue",
                 font = 2)
  # F statistic and p value from permutation test
  graphics::mtext(side = 3,
                  adj = 1,
                  text = paste0("F=", round(x$anova["F"], 3), ", p=", x$anova["p"]))
  # add species info
  graphics::points(x = rda_species[,1],
                   y = rda_species[,2],
                   col = seq_along(rownames(rda_species)),
                   pch = 19,
                   cex = 1.5)
  graphics::text(x = rda_species[,1],
                 y = rda_species[,2],
                 labels = short_names(rownames(rda_species)),
                 pos = 3,
                 font = 2)
  # add arrows indicating species trend
  spec_arrow <- rda_species[,1] + 0.2 * rda_species[,1] * vegan::scores(x$rda)$regression[1]
  graphics::arrows(x0 = rda_species[,1],
                   y0 = rda_species[,2],
                   y1 = rda_species[,2],
                   x1 = spec_arrow,
                   col = seq_along(rownames(rda_species)),
                   length = 0.08)
}