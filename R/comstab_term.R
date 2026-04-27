# Functions to decompose community stability

#' Decompose community stability
#'
#' `comstab_term()` partitions the temporal coefficient of variation of a community into the variability of the average species and three stabilizing effects: the dominance, asynchrony and averaging effects. It allows standard estimates of variance and CV as well as their detrended versions using Hill's two and three term local quadratic variance estimates (see Details).
#' 
#' @param x A data.frame. A community matrix of species abundances with time in rows and taxa in columns. Optionally it can include community and time columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' 
#' @details The analytic framework is described in detail in Segrestin *et al.* (2024).
#' In short, the partitioning relies on the following equation: \deqn{CV_{com} = CV_e \Delta \Psi \omega} 
#' where \eqn{CV_{com}} is the community coefficient of variation (reciprocal of community stability), 
#' \eqn{CV_e} is the expected community CV when controlling for the dominance structure and species temporal synchrony,
#' \eqn{ \Delta} is the dominance effect, \eqn{ \Psi} is the asynchrony effect, and \eqn{ \omega} is the averaging effect.
#' 
#' @returns An object of class `comstab`, a list of named vectors containing the following components:
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
    stop("The community CV is zero. This analysis does not apply to \n perfectly stable communities.")
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
  CVe <- 10^TPL["alpha"] * (mean(x)^TPL["beta"]) # Predict CVe from mean abundance and TPL coefficients (backtransformed from log scale)

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

#' @export
print.comstab <- function(x, ...){
  cat("\nPartitionning of the community temporal variability (CV)")
  cat("\n")
  cat(paste0("Community CV = ", round(x$CVs["CVc"], 2),
             "\nTotal stabilization = ", round(x$Stabilization["tau"], 2),
             "\nDominance effect = ", round(x$Stabilization["Delta"], 2),
             "\nAsynchrony effect = ", round(x$Stabilization["Psi"], 2),
             "\nAveraging effect = ", round(x$Stabilization["omega"], 2)))
  cat("\n")
  cat("\nRelatives contributions:")
  cat(paste0("\n% Dominance = ", round(x$Relative["Delta_cont"], 2)),
      paste0("\n% Asynchrony = ", round(x$Relative["Psi_cont"], 2)),
      paste0("\n% Averaging = ", round(x$Relative["omega_cont"], 2)))
}

#' @export
as.data.frame.comstab <- function(x, ...) {
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
