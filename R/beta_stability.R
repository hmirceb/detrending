#' Chord transformation for a matrix of abundances
#'
#' Applies chord transformation to a matrix (see Details).
#'
#' @param x A community abundance matrix.
#'
#' @details
#' Chord transformation standardizes values by dividing values in each row by the sample norm (*SN*) of the corresponding row following: 
#' \deqn{SN = \sqrt{\sum_{i=1}^{S} x_i^2}}
#'
#' @returns A matrix of community abundance after chord transformation.
#'
#' @references 
#' - Orlóci, L. (1967). An agglomerative method for classification of plant communities. The Journal of Ecology, 193-206.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Load and clean data
#' data(example_data_wide)
#' metacomm_df <- clean_community(x = example_data_wide,
#'                                input_format = "wide",
#'                                community_col = "comm",
#'                                time_col = "time",
#'                                filter_transient = FALSE)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)] # Select only community 1
#' # Transform
#' chord_transform(x = comm_df)
#' @export
chord_transform <- function(x) {
  c_t <- x / sqrt(rowSums(x^2))
  return(c_t)
}

#' Multivariate variance of community composition
#'
#' @param x A community abundance matrix.
#' @param method Character. Community dissimilarity metric to use. One of "euclidean" or "chord". Default "euclidean".
#'
#' @details
#' The multivariate variance of community composition is defined as the average square Euclidean distance between annual observations (\eqn{X_{i}}) and the average community composition (\eqn{\overline{X}}) following: 
#' \deqn{var_{mv} = \dfrac{ \sum_{i=1}^{t}{ED(X_{i}, \overline{X})^2} }{t-1}}
#'
#' @returns A numeric value. The multivariate variance of the community.
#' 
#' @author Jan Lepš, \email{suspa@@prf.jcu.cz}
#' @author Aleš Lisner, \email{lisnea00@@jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Load and clean data
#' data(example_data_wide)
#' metacomm_df <- clean_community(x = example_data_wide,
#'                                input_format = "wide",
#'                                community_col = "comm",
#'                                time_col = "time",
#'                                filter_transient = FALSE)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate multivariate variance
#' var_mv(x = comm_df, method = "euclidean") # 183.752
#' var_mv(x = comm_df, method = "chord") # 0.399
#' @export
var_mv <- function(x, method = c("euclidean", "chord")){
  # Create DF with average abundance values per species in first row
  mean_vec <- colMeans(x)
  data_merged <- rbind(mean_vec, x)
  
  # Compute distances
  dist_mat <- as.matrix( vegan::vegdist(data_merged, method = match.arg(method)) )
  dist_eu <- dist_mat[-1,1]
  # First column is distances between average community and each year (remove first obs)
  sum_sq <- sum(dist_eu^2)
  mv_var <- sum_sq / (nrow(x)-1)
  
  return(mv_var)
}

#' Multivariate two term local quadratic variance of community composition
#'
#' @param x A community abundance matrix.
#' @param method Character. Community dissimilarity metric to use. One of "euclidean" or "chord". Default "euclidean".
#'
#' @details
#' The multivariate two term local quadratic variance (\eqn{TTLQV_{mv2}}) of community composition is the detrended version of multivariate variance (\eqn{var_{mv}}): 
#' \deqn{TTLQV_{mv2} = \dfrac{ \sum_{i=1}^{t-1}{ED(X_{i}, X_{i+1})^2} }{2(t-1)}}
#' Where \eqn{X_{i}} is the composition of the community at time \eqn{i} and \eqn{ED(X_{i}, X_{i+1})} the Euclidean distance between succesive time points.
#'
#' @returns A numeric value. The multivariate two term quadratic variance of the community.
#' 
#' @author Jan Lepš, \email{suspa@@prf.jcu.cz}
#' @author Aleš Lisner, \email{lisnea00@@jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Load and clean data
#' data(example_data_wide)
#' metacomm_df <- clean_community(x = example_data_wide,
#'                                input_format = "wide",
#'                                community_col = "comm",
#'                                time_col = "time",
#'                                filter_transient = FALSE)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate multivariate variance
#' var_t2mv(x = comm_df, method = "euclidean") # 193.893
#' var_t2mv(x = comm_df, method = "chord") # 0.370
#' @export
var_t2mv <- function(x, method = c("euclidean", "chord")){
  # Compute distances
  dis <- as.matrix( vegan::vegdist(x, method = match.arg(method)) )
  # Get superdiagonal (distance between consecutive years)
  dis <- dis[row(dis) == col(dis) + 1]
  # TTQVmv
  mv_var <- sum(dis^2) / (2*(nrow(x)-1))
  
  return(mv_var)
}

#' Multivariate coefficient of variation
#'
#' @param x A community abundance matrix.
#' @param method Character. Community dissimilarity metric to use. One of "euclidean" or "chord". Default "euclidean".
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance) "two" for Hill's two term local quadrat variance. Default "var".
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#'
#' @details The multivariate coefficient of variation is estimated by dividing the multivariate variance of community composition by the sample norm of the average composition.
#' 
#' @returns A numeric value. The multivariate coefficient of variation of the community estimated using the variance function of choice.
#' 
#' @author Jan Lepš, \email{suspa@@prf.jcu.cz} 
#' @author Aleš Lisner, \email{lisnea00@@jcu.cz}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Load and clean data
#' data(example_data_wide)
#' metacomm_df <- clean_community(x = example_data_wide,
#'                                input_format = "wide",
#'                                community_col = "comm",
#'                                time_col = "time",
#'                                filter_transient = FALSE)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate multivariate variance
#' cv_mv(x = comm_df, method = "euclidean", term = "var") # 0.958
#' cv_mv(x = comm_df, method = "euclidean", term = "two") # 0.961
#' @export
cv_mv <- function(x, time_col = "time", method = c("euclidean", "chord"), term = c("var", "two")) {
  method <- match.arg(method)
  term   <- match.arg(term)
  
  # Check if a time column was specified for detrending methods and order rows
  x <- check_time(x, time_col = time_col, term = term, rm = TRUE)
  
  # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
  x <- remove_empty_sps(x = x, time_col = time_col)
  
  # Match variance function
  var_func <- switch(
    term,
    var = var_mv,
    two = var_t2mv
  )
  
  # Estimate MV variance
  vv <- var_func(x = x, 
                 method = method)
  # Average composition
  mean_vec <- colMeans(x)
  
  # Sample norm of average composition
  SN <- sqrt(sum(mean_vec^2))
  # CVmv
  CVmv <- switch(method,
                 chord     = sqrt(vv),
                 euclidean = sqrt(vv) / SN
  )
  names(CVmv) <- paste0("CV_", method)
  
  return(CVmv)
}