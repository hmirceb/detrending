# Functions to calculate community (a)synchrony

#' Psi synchrony index
#' 
#' Estimate \eqn{\psi}, 
#' 
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional, by default assumes that rows are in chronological order.
#' 
#' @details
#' 
#' \eqn{\psi} is estimated as:
#' 
#' \deqn{\psi = \sqrt{ \phi } ^ \alpha = \left( \dfrac{ \sigma_{x_{T}} }{ \sum_{i=1}^{S}{\sigma_{x_{i}}} } \right) ^ \alpha } 
#' 
#' Where \eqn{\alpha} is a scaling factor estimated as:
#' 
#' \deqn{\alpha = \dfrac{ \log_{10}{1 / 2} }{ \log_{10}{ \sum_{i=1}^{S}{\sigma_{x_{i}}^2} / (\sum_{i=1}^{S}{\sigma_{x_{i}}})^2 } } } 
#' 
#' Where \eqn{\sigma_{x}} is the standard deviation of a vector of abundances \eqn{x}, \eqn{S} is the number of species in the community, \eqn{x_{i}} is the abundance of species \eqn{i} across time steps, \eqn{x_{T}} is the sum of species abundances for each time step.
#' 
#' @returns A numeric value.
#' 
#' @references
#' - Segrestin, J., Götzenberger, L., Valencia, E., de Bello, F., & Lepš, J. (2024). A unified framework for partitioning the drivers of stability of ecological communities. Global Ecology and Biogeography, 33(5), e13828.
#' 
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' psi_segrestin(x = comm_df$sim_data)
#' @export
psi_segrestin <- function(x, term = "var", time_col = "time"){
  # Match variance function
  var_func <- switch(
    match.arg(term, choices = c("var", "two", "three")),
    var = stats::var,
    two = var_t2,
    three = var_t3
  )
  
  # Check if a time column was specified for detrending methods
  x <- check_time(x = x, time_col = time_col, term = term, rm = TRUE)
  
  # Set NAs as 0 and remove species with 0 abundance
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  varsum <- var_func(rowSums(x)) # variance of sum of abundances
  sdsum <- sqrt(varsum) # Square root of sum of variances (SD). Equivalent to SD of the sum of yearly abundances (whole community)
  
  # Variance of each species
  vari <- apply(X = x, MARGIN = 2, FUN = var_func) 
  sumsd <- sum(sqrt(vari))
  
  rootPhi <- sdsum / sumsd # Ratio between SD of whole community vs sum of individual SDs
  sumvar <- sum(vari) # Sum of individual variances
  alpha <- log10(1/2) / (log10(sumvar/(sumsd^2))) # Scaling coefficient (eq. 7)
  Psi <- rootPhi^alpha # Asynchrony effect
  
  return(Psi)
}

#' Phi synchrony index
#' 
#' This function estimates Loreau & Mazancourts 2008 Phi synchrony index using standard and detrended versions of variances based on Hill's 2 and 3 terms local quadratic variance. 
#' 
#' @param x A community matrix of species abundance with years as rows and species as columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional, by default assumes that rows are in chronological order.
#' 
#' @returns A numeric value.
#'
#' @references
#' - Loreau, M., & de Mazancourt, C. (2008). Species synchrony and its drivers: neutral and nonneutral community dynamics in fluctuating environments. The American Naturalist, 172(2), E48-E66.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' phi_loreau(x = comm_df$sim_data)
#' @export
phi_loreau <- function(x, term = "var", time_col = "time") {
  # Match variance function
  var_func <- switch(
    match.arg(term, choices = c("var", "two", "three")),
    var = stats::var,
    two = var_t2,
    three = var_t3
  )
  
  # Check if a time column was specified for detrending methods
  x <- check_time(x = x, time_col = time_col, term = term, rm = TRUE)
  
  # Set NAs as 0 and remove species with 0 abundance
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  var_com <- var_func(rowSums(x)) # Variance of sum of abundances
  var_sps <- sum(apply(x, 2, function(y) sqrt(var_func(y))))^2 # Squared sum of individual SDs 
  sync <- var_com/var_sps # Ratio between variance of whole community and sum of individual vars
  return(sync)
}

#' Eta synchrony index
#'
#' This function estimates Gross et al. 2014 synchrony index Eta and its weighted version by Blüthgen et al. 2016 using a detrended version of variances based on Hill's 2 and 3 terms local quadratic variance.
#'
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional, by default assumes that rows are in chronological order.
#' @param weighted Boolean. Weight index by average species relative abundances. Default FALSE.
#' 
#' @returns A numeric value.
#'
#' @references
#' - Gross, K., Cardinale, B. J., Fox, J. W., Gonzalez, A., Loreau, M., Wayne Polley, H., ... & van Ruijven, J. (2014). Species richness and the temporal stability of biomass production: a new analysis of recent biodiversity experiments. The American Naturalist, 183(1), 1-12.
#' - Blüthgen, N., Simons, N. K., Jung, K., Prati, D., Renner, S. C., Boch, S., ... & Gossner, M. M. (2016). Land use imperils plant and animal community stability through changes in asynchrony rather than diversity. Nature Communications, 7(1), 10697.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' eta_gross(x = comm_df$sim_data)
#' @export
eta_gross <- function(x, term = "var", time_col = "time", weighted = FALSE) {
  # Match variance function
  var_func <- switch(
    match.arg(term, choices = c("var", "two", "three")),
    var = stats::var,
    two = var_t2,
    three = var_t3
  )
  
  # Check if a time column was specified for detrending methods
  x <- check_time(x = x, time_col = time_col, term = term, rm = TRUE)
  
  # Set NAs as 0 and remove species with 0 abundance
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # Correlations
  corrs <- c() # Vector to store correlations
  for (i in seq_len(ncol(x))) { # Loop over columns
    sp_ab <- x[,i] # Abundance of focal species
    
    # Sum of all other species
    # use appropriate sum if there are only 2 species in the community 
    if ( ncol(x) == 2) {
      com_ab <- x[,-i]
    } else {
      com_ab <- rowSums(x[,-i])
    }
    # Check weighted argument
    if( isTRUE(weighted) ){
      # Relative abundance of species i across all years
      w <- colSums(x)[i] / sum(colSums(x))
    } else {
      # Equal weights otherwise
      w <- 1 / ncol(x)
    }
    # Multiply by weights
    corrs[i] <- w * cor_term(sp_ab, com_ab, term = term) # Correlation between species i and rest of community
  }
  # Mean correlation
  sync <- sum(corrs)
  return(sync)
}

#' Log variance ratio synchrony index
#' 
#' This function lets you estimate the synchrony index using a detrended version of variances based on Hill's 2 and 3 terms local quadratic variance. 
#' 
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns. 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional, by default assumes that rows are in chronological order.
#' @param log Boolean. Apply the natural logarithm to the variance ratio. Default TRUE.
#' @returns A numeric value.
#' 
#' @references
#' - Lepš, J., Májeková, M., Vítová, A., Doležal, J., & de Bello, F. (2018). Stabilizing effects in temporal fluctuations: Management, traits, and species richness in high‐diversity communities. Ecology, 99(2), 360-371.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' logvar_ratio(x = comm_df$sim_data)
#' @export
logvar_ratio <- function(x, term = "var", time_col = "time", log = TRUE) {
  # Match variance function
  var_func <- switch(
    match.arg(term, choices = c("var", "two", "three")),
    var = stats::var,
    two = var_t2,
    three = var_t3
  )
  
  # Check if a time column was specified for detrending methods
  x <- check_time(x = x, time_col = time_col, term = term, rm = TRUE)
  
  # Set NAs as 0 and remove species with 0 abundance
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # Variance of sum of abundances
  var_com <- var_func(rowSums(x))
  # Sum of individual variances 
  var_sps <- sum(apply(x, 2, 
                       function(y) var_func(y))) 
  
  # Compute ratio
  v_ratio <- var_com/var_sps
  
  # Decadic logarithm of ratio variance of whole community and sum of individual vars
  if ( isTRUE(log) ) {
    sync <- log(v_ratio)
  } else {
    sync <- v_ratio
  }

  return(sync)
}


#' Detrended variance ratio synchrony index based on residuals from linear regression 
#' 
#' This function fits linear regression models to each species in the community. Then it estimaes the synchrony of the community using the variance ratio of the observed values, the fitted values from the regression and their residuals.
#' 
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns. 
#' @param time_col Character. Name of the column with time variable. Optional, by default assumes that rows are in chronological order.
#' @returns A named vector of numeric values with the synchrony of the community whole community, the fitted values and their residuals.
#' 
#' @references
#' - Lepš, J., Götzenberger, L., Valencia, E. & de Bello, F. (2019). Accounting for long-term directional trends on year-to-year synchrony in species fluctuations. Ecography, 42, 1728-1741.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' sumsq_ratio(x = comm_df$sim_data)
#' @export
sumsq_ratio <- function(x, time_col = "time") {
  
  # Check if a time column was specified for detrending methods
  x <- check_time(x = x, time_col = time_col, term = "two", rm = TRUE)
  
  # Set NAs as 0 and remove species with 0 abundance
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # Variance of sum of abundances
  var_com <- var_linear(rowSums(x))
  # Sum of individual variances of fitted values and residuals
  var_fitres <- rowSums(apply(dat, 2, var_linear))
  # Sum of individual variances 
  var_sps <- sum(apply(x, 2, var)) 
  
  # Compute ratio
  v_ratio <- (var_com-var_fitres)/var_sps
  
  # Include total synchrony
  sync <- c(sum(v_ratio), v_ratio)
  names(sync) <- c("sync_total", "sync_fitted", "sync_det")
  return(sync)
}

#' Compute community synchrony
#' 
#' `sync_term()` estimates one or several community synchrony indices (see Details) using standard estimates of variance as well as its detrended versions using Hill's two and three term local quadratic variance estimates.
#'
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns. 
#' @param index Character. Synchrony index to calculate. One of "psi", (Segrestin *et al.* 2024), "phi" (Loreau & Mazancourt 2008), "eta" (Gross *et al.* 2014) or "logvar" (Leps *et al.* 2018). 
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param linear Boolean. Decompose synchrony into its trend and non-trend components using linear regression. Default FALSE. 
#' @param time_col Character. Name of the column with time variable. Optional, if not provided the function assumes that rows are in order  by default.
#' @param weighted Boolean. Weight the contribution of each species to Gross *et al.*'s \eqn{\eta} index by its average abundance in the community following Blüthgen *et al.* (2016). Default FALSE.
#'
#' @details
#' There are five synchrony indices available:
#'
#' - Segrestin *et al.*'s (2024) \eqn{\psi}: 
#' \deqn{\psi = \sqrt{ \phi } ^ \alpha = \left( \dfrac{ \sigma_{x_{T}} }{ \sum_{i=1}^{S}{\sigma_{x_{i}}} } \right) ^ \alpha } 
#' \deqn{\alpha = \dfrac{ \log_{10}{1 / 2} }{ \log_{10} \left( \sum_{i=1}^{S}{\sigma_{x_{i}}^2} / (\sum_{i=1}^{S}{\sigma_{x_{i}}})^2 \right) } } 
#' 
#' - Loreau & Mazancourt's (2008) \eqn{\phi}:
#' \deqn{\phi = \dfrac{ \sigma_{x_{T}}^2 }{ (\sum_{i=1}^{S}{\sigma_{x_{i}}})^2 }}
#' 
#' - Gross *et al.*'s (2014) \eqn{\eta}:
#' \deqn{\eta = \dfrac{ 1 }{ S }\sum_{i=1}^{S}{corr(x_{i},\sum_{j \neq i}^{S}{x_{j})}}}
#' 
#' - Blüthgen *et al.*'s (2016) weighted version of \eqn{\eta}, \eqn{\eta_{w}}:
#' \deqn{\eta_{w} = \sum_{i=1}^{S}{p_{i} corr(x_{i},\sum_{j \neq i}^{S}{x_{j})}}}
#' 
#' - Lepš *et al.* (2018) *logvar* ratio:
#' \deqn{logvar = \log_{10} \left( \dfrac{ var(\sum_{i=1}^{S}{x_{i}}) }{ \sum_{i=1}^{S}{var(x_{i})} } \right)}
#' 
#' - Lepš *et al.* (2019) variance ratio of residuals from linear regression:
#' \deqn{Stotal = Strend + Sdetrended = Sfitted + Sresiduals}
#' \deqn{Stotal = \dfrac{ var\left( \sum_{i=1}^{S}{x_{i}} \right) - \sum_{i=1}^{S}{var(x_{i})}}{ \sum_{i=1}^{S}{var(x_{i})} }} 
#' \deqn{Stotal = \dfrac{ var\left( \sum_{i=1}^{S}{fitted_{i}} \right) - \sum_{i=1}^{S}{var(fitted_{i})}}{ \sum_{i=1}^{S}{var(fitted_{i})} } + \dfrac{ var\left( \sum_{i=1}^{S}{residual_{i}} \right) - \sum_{i=1}^{S}{var(residual_{i})}}{ \sum_{i=1}^{S}{var(residual_{i})} }}
#' 
#' And:
#' \deqn{ logvar = \ln (1+Stotal) }
#' 
#' Where \eqn{\sigma_{x}} is the standard deviation of a vector of abundances \eqn{x}, \eqn{S} is the number of species in the community, \eqn{x_{i}} is the abundance of species \eqn{i} across time steps, \eqn{x_{T}} is the sum of species abundances for each time step and \eqn{p_{i}} the average relative abundance of species \eqn{i}. For `Stotal`, \eqn{fitted_{i}} and \eqn{residual_{i}} correspond, respectively, to the fitted and residual values of a linear regression between abundances and time. 
#' 
#' @returns A named vector of length equal to the number of indices calculated.
#'
#' @references
#' - Segrestin, J., Götzenberger, L., Valencia, E., de Bello, F., & Lepš, J. (2024). A unified framework for partitioning the drivers of stability of ecological communities. Global Ecology and Biogeography, 33(5), e13828.
#' - Loreau, M., & de Mazancourt, C. (2008). Species synchrony and its drivers: neutral and nonneutral community dynamics in fluctuating environments. The American Naturalist, 172(2), E48-E66.
#' - Gross, K., Cardinale, B. J., Fox, J. W., Gonzalez, A., Loreau, M., Wayne Polley, H., ... & van Ruijven, J. (2014). Species richness and the temporal stability of biomass production: a new analysis of recent biodiversity experiments. The American Naturalist, 183(1), 1-12.
#' - Blüthgen, N., Simons, N. K., Jung, K., Prati, D., Renner, S. C., Boch, S., ... & Gossner, M. M. (2016). Land use imperils plant and animal community stability through changes in asynchrony rather than diversity. Nature Communications, 7(1), 10697.
#' - Lepš, J., Májeková, M., Vítová, A., Doležal, J., & de Bello, F. (2018). Stabilizing effects in temporal fluctuations: Management, traits, and species richness in high‐diversity communities. Ecology, 99(2), 360-371.
#' - Lepš, J., Götzenberger, L., Valencia, E. & de Bello, F. (2019). Accounting for long-term directional trends on year-to-year synchrony in species fluctuations. Ecography, 42, 1728-1741.
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' sync_term(x = comm_df$sim_data)
#' @export
sync_term <- function(x,
                      index = c("psi", "phi", "eta", "logvar"),  
                      term = "var",
                      linear = FALSE,
                      time_col = "time", 
                      weighted = FALSE) {
  
  # Check if proper indices were selected
  if(sum(index %in% c("psi", "phi", "eta", "logvar")) == 0) {
    stop("Please choose an appropriate synchrony index")
  }
  # Match arguments with multiple synchrony functions
  options <- data.frame(index = c("psi", "phi", "eta", "logvar"),
                        fun = c("psi_segrestin", "phi_loreau", "eta_gross", "logvar_ratio"))
  index_func <- options[options$index %in% index,]$fun
  
  # Check if a time column was specified for detrending methods
  x <- check_time(x = x, time_col = time_col, term = term, rm = FALSE)
  
  sync <- sapply(index_func, function(f){
    syn_func <- match.fun(f) # match function
    # Run functions
    # Correctly use the weighted argument
    if( f == "eta_gross" & isTRUE(weighted) ){
      sync <- syn_func(x = x, 
                       term = term,
                       time_col = time_col,
                       weighted = weighted)
    } else{
      sync <- syn_func(x = x,
                       time_col = time_col,
                       term = term)
    }
    return(sync)
  })
  
  # Set names with indices and term
  if( isTRUE(weighted) ){
    names(sync) <- gsub("eta", "etaw", names(sync))
  }
  names(sync) <- paste(sep = "_", names(sync), term)
  
  # estimate linear synchrony decomposition
  if( isTRUE(linear) ){
    linear_sync <- sumsq_ratio(x, time_col = time_col)
    sync <- c(sync, linear_sync)
  }
  
  return(sync)
}

