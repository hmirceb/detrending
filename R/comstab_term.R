# Functions to decompose community stability

#' Decompose community stability (internal function)
#'
#' @noRd
comstab_internal <- function(x, 
                        term = "var",
                        time_col = "time") {
  
  # Match variance function
  var_func <- switch(
    term,
    var = stats::var,
    two = var_t2,
    three = var_t3
  )
  
  # Errors if data is not properly formated
  if ( !is.data.frame(x) ) {
    stop("Error: x is not a data.frame")
  }
  # if ( !is.numeric(x) ) {
  #   stop("Error: non-numerical values in x")
  # } 
  if ( any(x < 0) ) {
    stop("Error: negative values in x")
  }
  if ( nrow(x) == 1 ) {
    stop("Only one year provided.")
  }
  
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = term, rm = TRUE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # Remove columns (species) with constant abundance (min abundance == max abundance)
  x <- x[, apply(X = x, MARGIN = 2, FUN = min) != 
           apply(X = x, MARGIN = 2, FUN = max), 
         drop = FALSE]
  
  # Turn into matrix
  x <- as.matrix(x)
  
  # number of species
  n <- ncol(x) 
  
  ## Community metrics ##
  varsum <- var_func(rowSums(x)) # variance of sum of abundances
  meansum <- mean(rowSums(x)) # mean of sum of abundances
  CV <- sqrt(varsum) / meansum # CV of sum of abundances
  
  # Check if community fluctuates
  if (CV == 0) {
    stop("The community CV is zero. This analysis does not apply to \nperfectly stable communities.")
  }
  
  # Stop if there is only one species
  if (ncol(x) == 1) {
    stop("This analysis is not relevant for single-species communities.")
  }
  
  #--------------#
  # Partitioning #
  #--------------#
  
  # Estimate CVe from TPL of all species
  vari <- apply(X = x, MARGIN = 2, FUN = var_func) # Variance of each species
  meani <- colMeans(x) # Mean of each species
  CVi <- sqrt(vari) / meani  # CV of each species
    
  if (any(CVi == 0)) { # Warn if constant species present
    warning("Non-fluctuating species found in the data.")
    }
    
  CV0 <- which(CVi > 0) # Use only species with CV != 0
  # Calculate TPL between CVs and means. 
  # This is equivalent to estimating TPL for variances and means, estimating the variance of the average species using TPL and average mean and then estimating CVe with that variance
  TPL <- tpl(vari = CVi[CV0],  meani = meani[CV0]) # LM of CVs and means on log scale
  CVe <- 10^TPL["alpha"] * (meansum/n)^TPL["beta"] # Predict CVe from mean abundance and TPL coefficients (backtransformed from log scale)

  # Test correlation between individual CVs and mean abundances (if there are more than 5 species with variation)
  if (sum(CV0) > 5) {
    testcor <- stats::cor.test(log10(CVi[CV0]), log10(meani[CV0]))$p.value > 0.05
    if (testcor) {
      warning("No significant power law between species CVs and abundances.")
      }
  } else {
    warning("Low number of species. The power law between species CVs and abundances cannot be tested.")
  }
    
  ## Dominance effect #
  sumsd <- sum(sqrt(vari)) # sum of individual SDs
  CVtilde <- sumsd / meansum # CV tilde. Weighted mean of individual CVs. sum(pi * sdi/mui) = sum(mui/meansum * sdi/mui) = sum(sdi/meansum) = sum(sdi)/meansum
  Delta <- CVtilde / CVe # Ratio CVtilde / CVe "average species"
  if (Delta > 1) {
    warning("Destabilizing effect of dominants. Relative effects cannot be computed.")
    }
  
  ## Compensatory dynamics ##
  sdsum <- sqrt(varsum) # Square root of sum of variances (SD). Equivalent to SD of the sum of yearly abundances (whole community)
  rootPhi <- sdsum / sumsd # Ratio between SD of whole community vs sum of individual SDs
  sumvar <- sum(vari) # Sum of individual variances
  alpha <- log10(1/2) / (log10(sumvar/(sumsd^2))) # Scaling coefficient (eq. 7)
  Psi <- rootPhi^alpha # Asynchrony effect
  omega <- rootPhi / Psi # Diversity effect
  if (omega > 1) {
    warning("Community diversity is lower than the null diversity. Relative effects cannot be computed.")
  }
  
  ## Partitioning ##
  tau <- Delta * Psi * omega
  CVs <- stats::setNames(object = c(CVe, CVtilde, CVtilde * 
                                      Psi, CV), nm = c("CVe", "CVtilde", "CVa", "CVc"))
  Stabilization <- stats::setNames(object = c(tau, Delta, 
                                              Psi, omega), nm = c("tau", "Delta", "Psi", "omega"))
  if (any(Stabilization > 1)) {
    Relative <- stats::setNames(object = rep(NA, 3), 
                                nm = c("Delta_cont", "Psi_cont", "omega_cont"))
  } else { # Return relative importance of each component
    Relative <- stats::setNames(object = c(log10(Delta) / log10(tau), # dominance
                                           log10(Psi) / log10(tau), # asynchrony
                                           log10(omega) / log10(tau)), # averaging 
                                nm = c("Delta_cont", "Psi_cont", "omega_cont"))
  }
  # Results into a list
  res <- list(CVs = CVs, Stabilization = Stabilization, 
              Relative = Relative)
  class(res) <- "comstab"
  return(res)
}

#' Decompose community stability
#'
#' `comstab_term()` partitions the temporal coefficient of variation of a community into the variability of the average species and three stabilizing effects: the dominance, asynchrony and averaging effects. It allows standard estimates of variance and CV as well as their detrended versions using Hill's two and three term local quadratic variance estimates (see Details).
#' 
#' @param x A data.frame. A community matrix of species abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param community_col Character. Name of column with the community identifier.
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' 
#' @details The analytic framework is described in detail in Segrestin *et al.* (2024).
#' In short, the partitioning relies on the following equation: \deqn{CV_{com} = CV_e \Delta \Psi \omega} 
#' where \eqn{CV_{com}} is the community coefficient of variation (reciprocal of community stability), 
#' \eqn{CV_e} is the expected community CV when controlling for the dominance structure and species temporal synchrony,
#' \eqn{ \Delta} is the dominance effect, \eqn{ \Psi} is the asynchrony effect, and \eqn{ \omega} is the averaging effect.
#' 
#' @returns A single object of class `comstab` or multiple of them if data correspond to several communities, each `comstab` object is a list of named vectors containing the following components:
#'  - `CVs`: a named vector of calculated coefficient of variations. `CVe` is the CV of an average species,
#'  `CVtilde` is the mean of species CVs weighted by their relative abundances, `CVa` is the expected community CV if 
#'   the community was stabilized by species asynchrony only, and `CVc` is the observed community CV.
#'   
#'  - `Stabilization`: a named vector of the stabilizing effects. `tau` is the total stabilization, `Delta` is
#'  the dominance effect, `Psi` is the asynchrony effect, and `omega` is the averaging effect.
#'  
#'  - `Relative`: a named vector of the relative contributions of each stabilizing effect to the total stabilization.
#'  `Delta_cont`, `Psi_cont`, and `omega_cont` are the relative contribution of respectively, the dominance, asynchrony, and averaging effects to the total stabilization.
#'  Returns a vector of NAs if any Stabilizing effect is higher than 1.
#'  
#' @references
#'  - Segrestin, J., Götzenberger, L., Valencia, E., de Bello, F., & Lepš, J. (2024). A unified framework for partitioning the drivers of stability of ecological communities. Global Ecology and Biogeography, 33(5), e13828.
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @author Jan Lepš, \email{suspa@@prf.jcu.cz} 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Load and clean data
#' comm_df <- sim_mvcomm(n_sp = 15, years = 30)
#' 
#' # Decompose CV into stability components
#' comstab_term(x = comm_df$sim_data, time_col = "time")
#' @export
comstab_term <- function(x, 
                          term = "var",
                          community_col = "comm",
                          time_col = "time"){
  x <- as.data.frame(x)
  
  # Check community column, if not present create one and assume a single community
  if( !community_col %in% colnames(x) ) {
    warning("Missing community column. Data are assumed to belong to a single community.",
            call. = FALSE)
    community_col <- "comm"
    x <- cbind(comm = as.character(rep(1, times = nrow(x))), x)
  }
  
  # split into communities
  c_list <- split(x, f = as.character(x[, community_col]))
  
  # apply comstab to each community
  comstab <- lapply(c_list, FUN = function(y){
    # remove community column
    y <- y[,!colnames(y) %in% community_col]
    # apply comstab
    comstab <- comstab_internal(x = y,
                            term = term,
                            time_col = time_col)
    return(comstab)
  }
  )
  
  # join and add community id
  # comstab_df <- do.call("rbind", lapply(comstab, as.data.frame))
  # comstab_df <- cbind(comm = unique(x$comm), comstab_df)
  # colnames(comstab_df)[1] <- community_col
  # rownames(comstab_df) <- NULL
  # 
  # print(comstab_df)
  
  # set class and term for plot and print methods
  class(comstab) <- c("comstab", "comstab_list")
  attr(comstab, "term") <- term
  
  return(comstab)
}

#' @export
print.comstab <- function(x, ...){
  
  prt_comstab <- function(y){
    cat("\nPartitionning of the community temporal variability (CV)")
    cat("\n")
    cat(paste0("Community CV = ", round(y$CVs["CVc"], 2),
               "\nTotal stabilization = ", round(y$Stabilization["tau"], 2),
               "\nDominance effect = ", round(y$Stabilization["Delta"], 2),
               "\nAsynchrony effect = ", round(y$Stabilization["Psi"], 2),
               "\nAveraging effect = ", round(y$Stabilization["omega"], 2)))
    cat("\n")
    cat("\nRelatives contributions:")
    cat(paste0("\n% Dominance = ", round(y$Relative["Delta_cont"], 2)),
        paste0("\n% Asynchrony = ", round(y$Relative["Psi_cont"], 2)),
        paste0("\n% Averaging = ", round(y$Relative["omega_cont"], 2)))
  }
  
  if (inherits(x, "comstab_list")) {
    nms <- names(x)
    for (i in seq_along(x)) {
      label <- if (!is.null(nms) && nms[i] != "") nms[i] else paste("Object", i)
      if( i == 1 ){
        cat(paste0("Community: ", label, "\n"))
      } else {
        cat(paste0("\n\nCommunity: ", label, "\n"))
      }
      prt_comstab(x[[i]])
    }
  } else {
    prt_comstab(x)
  }
}

#' @export
as.data.frame.comstab <- function(x, ...){
  # define auxiliary internal function
  comstab_to_df <- function(x) {
    d <- data.frame(CVc = x$CVs[4], 
                    CVe = x$CVs[1],
                    CVtilde = x$CVs[2],
                    CVa = x$CVs[3],
                    tau = x$Stabilization["tau"],
                    delta = x$Stabilization["Delta"],
                    psi = x$Stabilization["Psi"],
                    omega = x$Stabilization["omega"],
                    delta_rel = x$Relative["Delta_cont"],
                    psi_rel = x$Relative["Psi_cont"],
                    omega_rel = x$Relative["omega_cont"])
    rownames(d) <- NULL
    return(d)
  }
  
  # apply based on class
  if (inherits(x, "comstab_list")) {
    dat <- do.call("rbind", lapply(x, comstab_to_df))
  } else {
    dat <- comstab_to_df(x)
  }
  dat <- cbind(comm = rownames(dat), term = attributes(x)$term, dat)
  rownames(dat) <- NULL
  return(dat)
}

#' @export
plot.comstab <- function(x, y = NULL, change = TRUE, relative = TRUE, ...) {
  
  # as data frame based on class
  dat <- as.data.frame(x)
  
  # pot functions
  # cv
  plot_cv <- function(dat) {
    cc <- dat[, c("CVe", "CVtilde", "CVa", "CVc")]
    plot(x = 1:4, y = NULL, type = "n", xaxt = "n", xlab = NA, ylab = NA,
         xlim = c(1, 4), ylim = c(0.9 * min(cc), 1.1 * max(cc)))
    for (i in seq_len(nrow(cc))) {
      graphics::lines(x = 1:4, y = cc[i, ], type = "b", col = i, pch = 19, cex = 1.5)
    }
    graphics::axis(side = 1, at = 1:4,
                   labels = c(expression(CV[e]), expression(widetilde(CV)),
                              expression(CV[a]), expression(CV[com])))
    graphics::mtext(text = c("Dominance", "Asynchrony", "Averaging"),
                    side = 1, at = 1.5:3.5, line = -1.1, cex = 1)
  }
  # relative effects
  plot_ternary <- function(dat) {
    # set colors (skip if no relative effects to match other plot)
    colors <- seq_len(nrow(dat)) * (dat$delta_rel / dat$delta_rel)
    colors <- colors[!is.na(colors)]
    # empty plot
    isopleuros::ternary_plot(NULL,
                             ann = FALSE,
                             axes = FALSE,
                             panel.first = isopleuros::ternary_grid(10),
                             
    )
    # change axis names and matching colors
    isopleuros::ternary_axis(side = 1, col = "#BB5566", at = seq(0.1, 0.9, 0.1))
    isopleuros::ternary_title(xlab = "Dominance",  col.lab = "#BB5566")
    isopleuros::ternary_axis(side = 2, col = "#004488", at = seq(0.1, 0.9, 0.1))
    isopleuros::ternary_title(ylab = "Asynchrony", col.lab = "#004488")
    isopleuros::ternary_axis(side = 3, col = "#DDAA33", at = seq(0.1, 0.9, 0.1))
    isopleuros::ternary_title(zlab = "Averaging",  col.lab = "#DDAA33")
    # draw points
    isopleuros::ternary_points(
      x = dat$delta_rel,
      y = dat$psi_rel,
      z = dat$omega_rel,
      pch = 19, 
      col = colors
    )
  }
  
  # warn if relative effects not available and avoid plotting them
  if ( nrow(dat) == 1 & any(is.na(dat$delta_rel)) ){
    relative <- FALSE
    warning("Relative effects for some communities could not computed and were not plotted.")
  }
  if ( nrow(dat) > 1 & all(is.na(dat$delta_rel)) ){
    relative <- FALSE
    warning("Relative effects could not computed and cannot be plotted.")
  }
  
  # set layout and plot
  graphics::par(xpd = NA) # allow plotting outside area for axis names
  if ( isTRUE(change) && isTRUE(relative) ) {
    # two columns
    graphics::layout(matrix(c(1, 2), nrow = 1, ncol = 2))
    
    plot_cv(dat)
    plot_ternary(dat)
  } else {
    # one column
    graphics::layout(matrix(1, nrow = 1, ncol = 1))
    
    if ( isTRUE(change) ) {
      plot_cv(dat) 
    } else {
      plot_ternary(dat)
    }
  }
  # reset graphics
  graphics::par(mfrow = c(1, 1), xpd = FALSE, mar = c(5.1, 4.1, 4.1, 2.1))
}