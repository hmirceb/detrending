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
#' metacomm_df <- clean_community_wide(x = example_data_wide)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)] # Select only community 1
#' # Transform
#' chord_transform(x = comm_df)
#' 
#' @export
chord_transform <- function(x) {
  c_t <- x / sqrt(rowSums(x^2))
  return(c_t)
}

#' Multivariate variance of community composition
#'
#' @param x A community abundance matrix.
#' @param d Character. Community dissimilarity metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
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
#' metacomm_df <- clean_community_wide(x = example_data_wide)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate multivariate variance
#' var_mv(x = comm_df, d = "euclidean") # 183.752
#' var_mv(x = comm_df, d = "chord") # 0.399
#' 
#' @export
var_mv <- function(x, d = "euclidean"){
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  
  # Create DF with average abundance values per species in first row
  mean_vec <- colMeans(x)
  data_merged <- rbind(mean_vec, x)
  
  # Apply chord transformation if necessary
  if( d == "chord" ) {
    data_merged <- chord_transform(data_merged)
  }
  # Compute distances
  dist_mat <- as.matrix( stats::dist(data_merged) )
  dist_eu <- dist_mat[-1,1]
  # First column is distances between average community and each year (remove first obs)
  sum_sq <- sum(dist_eu^2)
  mv_var <- sum_sq / (nrow(x)-1)
  
  return(mv_var)
}

#' Multivariate two term local quadratic variance of community composition
#'
#' @param x A community abundance matrix.
#' @param d Character. Community dissimilarity metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
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
#' metacomm_df <- clean_community_wide(x = example_data_wide)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate multivariate variance
#' var_t2mv(x = comm_df, d = "euclidean") # 193.893
#' var_t2mv(x = comm_df, d = "chord") # 0.370
#' 
#' @export
var_t2mv <- function(x, d = "euclidean"){
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  
  # Apply chord transformation if necessary
  if( d == "chord" ) {
    x <- chord_transform(x)
  }
  # Compute distances
  dis <- as.matrix( stats::dist(x) )
  # Get superdiagonal (distance between consecutive years)
  dis <- dis[row(dis) == col(dis) + 1]
  # TTQVmv
  mv_var <- sum(dis^2) / (2*(nrow(x)-1))
  
  return(mv_var)
}

# #' Multivariate three term local quadratic variance of community composition
# #'
# #' @param x A community abundance matrix.
# #' @param d Character. Community dissimilarity metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
# #'
# #' @returns
# #'
# #' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
# #'
# #' @export
# #'
# var_t3mv <- function(x, method = c("euclidean", "chord")){
#   if( !d %in% c("euclidean", "chord") ){
#     stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
#   }
# 
#   # Apply chord transformation if necessary
#   if( method == "chord" ) {
#     x <- chord_transform(x)
#   }
#   # Compute distances
#   dis <- as.matrix(dist(x))
#   # Get superdiagonal (distance between consecutive years)
#   dis <- dis[row(dis) == col(dis) + 1]
#   mv_var <- var_t3(dis)
# 
#   return(mv_var)
# }

#' Multivariate coefficient of variation
#'
#' @param x A community abundance matrix.
#' @param d Character. Community dissimilarity metric to use. One of 'euclidean' or 'chord'. Default 'euclidean'.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance) "two" for Hill's two term local quadrat variance. Default "var".
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
#' metacomm_df <- clean_community_wide(x = example_data_wide)
#' comm_df <- metacomm_df[metacomm_df$comm == 1,][,-c(1:2)]
#' 
#' # Calculate multivariate variance
#' cv_mv(x = comm_df, d = "euclidean", term = "var") # 0.785
#' cv_mv(x = comm_df, d = "euclidean", term = "two") # 0.807
#' 
#' @export
cv_mv <- function(x, d = "euclidean", term = "var") {
  if( !d %in% c("euclidean", "chord") ){
    stop("Unsuitable distance metric. Please choose one of 'euclidean' or 'chord'")
  }
  if( !term %in% c("var", "two") ){
    stop("Unsuitable variance term. Please choose one of 'var' or 'two'")
  }
  
  # Match variance function
  var_func <- switch(
    term,
    var = var_mv,
    two = var_t2mv
  )
  
  # Estimate MV variance
  vv <- var_func(x, d = d)
  # Average composition
  mean_vec <- colMeans(x)
  
  # Sample norm of average composition
  SN <- sqrt(sum(mean_vec^2))
  # CVmv
  if (d == "chord") {
    CVmv <- sqrt(vv)
  } else {
    CVmv <- sqrt(vv) / SN
  }
  return(CVmv)
}