# Functions to calculate detrended variances and covariances

#' Hill's two term local quadrat variance
#'
#' @param x Numeric. A vector of values to estimate the 2 term local quadrat variance.
#'
#' @details
#' The two term local quadrat variance is calculated using the following expression:
#' \deqn{TTLQV_{2} = \dfrac{\sum_{i=1}^{t-1}{(x_{i+1}-x_{i}})^2}{2(t-1)}}
#' 
#' Where \eqn{x_{i}} is the abundance of species \eqn{i} along time \eqn{t}. 
#'
#' @returns The 2 term local quadrat variance.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @references
#' - Hill, M. O. (1973). The intensity of spatial pattern in plant communities. The Journal of Ecology, 225-235.
#' 
#' @export
var_t2 <- function(x){
  # Compute the 2 term local variance between consecutive observations
  n <- length(x) # sample size
  vtwo_t <- diff(x)^2 # Squared difference between 2 consecutive data points xi and xj 
  vtwo <- sum(vtwo_t)/(2*(n-1)) # Average of all mean squared differences
  return(vtwo)
}

#' Hill's three term local quadrat variance 
#'
#' @param x Numeric. A vector of values to estimate their 3 term local quadrat variance.
#'
#' @details
#' The three term local variance is calculated using the following expression:
#' \deqn{TTLQV_{3} = \dfrac{\sum_{i=1}^{t-2}{(x_{i}-2x_{i+1}+x_{i+2}})^2}{6(t-2)}}
#' 
#' Where \eqn{x_{i}} is the abundance of species \eqn{i} along time \eqn{t}. 
#' 
#' @returns Numeric. The 3 term local quadrat variance
#'
#' @author Lars Götzenberger, \email{}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#'
#' @references
#' - Hill, M. O. (1973). The intensity of spatial pattern in plant communities. The Journal of Ecology, 225-235.
#' - Lepš, J., Götzenberger, L., Valencia, E., & de Bello, F. (2019). Accounting for long‐term directional trends on year‐to‐year synchrony in species fluctuations. Ecography, 42(10), 1728-1741.
#' 
#' @export
var_t3 <- function(x) {
  n <- length(x) # sample size
  multip  <- c(1, -2, 1) # signs for the sum in the numerator
  xsq <- c() # vector to populate
  for (i in 1:(n - 2)) {
    # calculate the squared sum at each step of the moving window
    xisq <- x[i:(i + 2)] * multip
    xsq[i] <- sum(xisq)^2
  }
  v <- sum(xsq) / (6 * (n - 2)) # three term local quadratic variance of x
  return(v)
}

#' Detrended variance using linear regression 
#'
#' @param x Numeric. A vector of values to estimate the detrended variance using a linear regression approach.
#'
#' @returns Numeric. The detrened variance using linear regression.
#' @export
var_linear <- function(x) {
  y <- 1:length(x)
  mod <- stats::lm(x ~ y)
  det_x <- stats::resid(mod)
  return( stats::var(det_x) )
}

#' Covariance between two variables using variance or Hill's 2 and 3 term variance
#'
#' @param x,y Numeric. A pair of vectors.
#' @param term One of "var", "two" or "three" for variance (default), Hill's quadratic variance of terms 2 or 3.
#'
#' @returns Numeric. The covariance between variables x and y.
#' 
#' @export
cov_term <- function(x, y, term = "var") {
  # Match argument for variance function to use
  options <- data.frame(term = c("var", "two", "three", "linear"),
                        var = c("var", "var_t2", "var_t3", "var_linear"))
  # Get choice
  var_func <- match.fun(options[options$term == term,]$var)
  
  cov <- (var_func(x + y) - var_func(x) - var_func(y)) / 2
  return(cov)
}

#' Standard, 2 and 3 term variance covariance matrix
#'
#' @param x A matrix with columns as variables (species) and rows as years
#' @param term Term to estimate the variance. One of "var" (for standard variance and covariane), "two" or "three" for Hills' two or three term local quadrat variance and covariance.
#'
#' @returns A matrix variance/covariance
#' 
#' @export
vcov_term <- function(x, term = "var") {
  
  # Match argument for variance function to use
  options <- data.frame(term = c("var", "two", "three", "linear"),
                        var = c("var", "var_t2", "var_t3", "var_linear"))
  var_func <- match.fun(options[options$term == term,]$var)
  
  # Calculate variances for all columns
  n <- ncol(x)
  vars <- sapply(1:n, function(i) var_func(x[, i]))
  
  # Covariance
  S <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      cov_ij <- (var_func(x[, i] + x[, j]) - vars[i] - vars[j]) / 2
      S[i, j] <- cov_ij
      S[j, i] <- cov_ij  # Symmetric matrix
    }
  }
  return(S)
}
