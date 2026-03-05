# Auxiliary functions

#' Check order of observations
#' 
#' This is an auxiliary function that checks if the observations (rows) in a community matrix are in chronological order.
#'
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns.
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' @param term Character. Term used to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param rm Boolean. If TRUE, removes the time column from the returned data.frame.
#' 
#' @returns A data.frame of community data.
#'
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
check_time <- function(x, time_col = "time", term = NULL, rm = TRUE) {
  
  # Check if a time column was specified for detrending methods
  if ( !time_col %in% colnames(x) & 
       term %in% c("two", "three") ) {
    warning("Missing time column. Rows are assumed to be in order for detrending.",
            call. = FALSE)
  } 
  
  # Add time column if missing
  if ( !time_col %in% colnames(x) ) {
    x <- as.data.frame( cbind(time = 1:nrow(x), x) )
    colnames(x)[1] <- time_col
  }
  
  # Reorder according to time
  #x <- x[with(x, order(x[, time_col])),]
  x <- x[order(x[[time_col]]), ]
  
  if ( isTRUE(rm) ) {
    # Remove time column once df is ordered
    id_cols <- colnames(x) %in% c(time_col)
    x <- x[,!id_cols]
  }
  
  return(x)
}

#' Remove empty columns
#'
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns.
#' @param time_col Character. Name of the column with time variable. Optional with default "time".
#' @param community_col Character. Name of the column with community variable. Optional with default "comm".
#'
#' @returns A data.frame with species with 0 abundance removed.
#'
remove_empty_sps <- function(x, time_col = "time", community_col = "comm") {
  # Set NAs as 0
  x[is.na(x)] <- 0
  # Get index of time column
  id_cols <- colnames(x) %in% c(community_col, time_col)
  # Sum abundances of species and check which ones are 0
  sps_to_remove <- colSums(x[, !id_cols]) == 0
  # Get names of species to remove
  sps_to_remove <- names(sps_to_remove)[sps_to_remove]
  # Remove species and time column from table
  x <- x[, !colnames(x) %in% sps_to_remove]
  return(x)
}

#' Jennings–Fischer formula estimates a combined value of plant cover assuming overlap between plants.
#'
#' @param x Numeric. A vector of cover values.
#' @param perc Boolean. If the cover values are expressed as percentages (0-100) or proportions (0-1). Default FALSE.
#'
#' @returns Numeric. 
#'
#' @references
#' - Jennings, M. D., Faber-Langendoen, D., Loucks, O. L., Peet, R. K., & Roberts, D. (2009). Standards for associations and alliances of the US National Vegetation Classification. Ecological Monographs, 79(2), 173-199.
#' - Fischer, H. S. (2015). On the combination of species cover values from different vegetation layers. Applied Vegetation Science, 18(1), 169-170.
#'
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#'
jenfish <- function(x, 
                    perc = FALSE) {
  # Check if cover is percentage
  if( isTRUE(perc) ) {
    x <- x/100
    c = (1 - prod(1 - x) )*100
  }
  c = (1 - prod(1 - x) )
  return(c)
}

#' Get dominant species
#'
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns.
#' @param q Numeric. Threshold of relative abundance to consider a species dominant.
#' @param plot Boolean. Plot the species-abundance curve of the community.
#'
#' @returns A data.frame with each species in the community, its mean abundance and if it is dominant or not.
#' @export
#'
get_dominants <- function(x, q = 0.9, plot = F) {
  # Sort species by their mean abundance across years
  sps_sorted <- sort(apply(x, 2,
                           function(y) mean(y[y > 0])), decreasing = T)
  
  # Get specified quantile
  qu = stats::quantile(sps_sorted, probs = q)[1]
  # Make DF and specify dominant species
  df = data.frame(taxon = names(sps_sorted),
                  abund = sps_sorted, 
                  dominant = as.factor(ifelse(sps_sorted >= qu, "yes", "no")))
  # Plot
  if( isTRUE(plot) ){
    # Points
    plot(df$abund, cex = 1.3, pch = 21, bg = df$dominant)
    # Labels
    graphics::text(sps_sorted, names(sps_sorted), cex=0.6, pos=1, col="red")
    # Threshold line 
    graphics::abline(h = qu, lty = "dashed")
  }
  return(df)
}

#' Check if dominant species have missing data
#'
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns.
#' @param q Numeric. Threshold of relative abundance to consider a species dominant.
#'
#' @returns A character vector. Names of the dominant species with missing values.
#' @export
#'
check_dominants <- function(x, q = 0.9) {
  # Get dominant species
  doms <- get_dominants(x = x, q = q, plot = FALSE)
  doms <- doms[doms$dominant == "yes",]$taxon
  # Check if they have missing data
  if( length(doms) == 1) {
    with_missing <- any( x[,colnames(x) %in% doms] == 0 )
    names(with_missing) <- doms
  } else {
    with_missing <- apply(x[,colnames(x) %in% doms], 2, 
                          function(y) any(y == 0))
  }
  
  # Get names of species
  missing_names <- names(with_missing[with_missing])
  if( length(missing_names) == 0) {
    message("No dominant species with missing values.")
  } else {
    message("The following dominant species have missing values.")
    return(missing_names)
  }
}

is_even <- function(x) {
  return( x %% 2 == 0 )
}

#' Plot a community time series
#'
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns.
#'
#' @returns A plot.
#' @export
#'
plot_com <- function(x) {
  # Plot first species
  plot(y = x[,1],
       x = 1:nrow(x),
       pch = 19,
       col = 1,
       ylim = c(min(x), max(x)),
       xlab = "time",
       ylab = "abundance")
  # Add additional points
  for (i in 2:ncol(x)) {
    graphics::points(y = x[,i],
           x = 1:nrow(x),
           pch = 19,
           col = i)
  }
}

#' Estimate of Taylor's Power Law
#'
#' @param vari Numeric. A vector of variances.
#' @param meani Numeric. A vector of means.
#'
#' @returns A named vector with coefficients a and b for the relation between variance and mean according to Taylor's Power Law.
#' @export
#'
tpl <- function(vari, meani) {
  # Check species with variance 0 because their log cannot be computed
  inds <- vari != 0
  # Remove them
  y <- vari[inds]
  x <- meani[inds]
  # Estimate TPL
  coefs <- stats::coef(stats::lm(log10(y) ~ log10(x)))
  names(coefs) <- c("alpha", "beta")
  return(coefs)
}
