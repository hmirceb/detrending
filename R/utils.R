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
#' @noRd
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
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#'
#' @noRd
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

#' Calculate total plant cover using Jennings–Fischer formula
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
#' @noRd
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
#' @param plot Boolean. Plot the species-abundance curve of the community. Default FALSE.
#'
#' @returns A data.frame with each species in the community, its mean abundance and if it is dominant or not.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @noRd
get_dominants <- function(x, q = 0.9, plot = FALSE) {
  # Sort species by their mean abundance across years
  sps_sorted <- sort(apply(x, 2,
                           function(y) mean(y[y > 0], na.rm = TRUE)), 
                     decreasing = TRUE)
  
  # Get specified quantile
  qu <- stats::quantile(sps_sorted, probs = q, na.rm = TRUE)[1]
  # Make DF and specify dominant species
  df <- data.frame(taxon = names(sps_sorted),
                  abund = sps_sorted, 
                  dominant = as.factor(ifelse(sps_sorted >= qu, "yes", "no")))
  # Plot
  if( isTRUE(plot) ){
    # Points
    plot(df$abund, cex = 1.3, pch = 21, bg = df$dominant, 
         ylim = c(min(df$abund, na.rm = TRUE), 1.1 * max(df$abund, na.rm = TRUE)))
    # Labels
    graphics::text(sps_sorted, labels = short_names(names(sps_sorted)), 
                   cex=0.6, 
                   pos=3,
                   col="black")
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
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @noRd
check_dominants <-  function(x, q = 0.9) {
  # Get dominant species
  doms <- get_dominants(x = x, q = q, plot = FALSE)
  doms <- doms[doms$dominant == "yes",]$taxon
  ab_doms <- x[,colnames(x) %in% doms]
  # Check if they have missing data
  if( length(doms) == 1 ) {
    miss <- sum(ab_doms == 0 | is.na(ab_doms))
    with_missing <- data.frame(taxon = doms, 
                               n_missing = miss)
  } else {
    miss <- apply(ab_doms, 2, 
                  function(y) sum(y == 0 | is.na(y)),
                  simplify = TRUE)
    with_missing <- data.frame(taxon = names(miss), 
                               n_missing = miss)
    rownames(with_missing) <- NULL # remove row names
  }
  
  # Get names of species
  missing_names <- unique(with_missing[with_missing$n_missing > 0,]$taxon)
  if( length(miss) == 0) {
    message("No dominant species with missing values.")
    return(with_missing)
  } else {
    message(paste0("The following dominant species have missing values: ",
                   paste(missing_names, collapse = ", ")))
    return(rbind(with_missing))
  }
}

#' Even number
#'
#' @param x Numeric. A number
#'
#' @returns Boolean.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @noRd
is_even <- function(x) {
  return( x %% 2 == 0 )
}

#' Plot a community time series
#'
#' @param x A data.frame. A community matrix of abundances with time in rows and taxa in columns.
#' @param total Boolean. Also plot the total abundance of the community by timestep. Default FALSE.
#' @param title Character. Title for the plot. Default NULL.
#' @returns A plot.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Simulate data
#' comm_data <- sim_mvcomm()
#' plot_com(comm_data$sim_data)
#' 
#' @export
plot_com <- function(x, total = FALSE, title = NULL) {
  if ( isTRUE(total) ) {
    # Plot total abundance
    plot(y = rowSums(x),
         x = 1:nrow(x),
         pch = 19,
         type = "l",
         lwd = 2,
         ylim = c(min(x, na.rm = T), max(rowSums(x), na.rm = T)),
         xlab = "Time",
         ylab = "Abundance")
    # Add additional points
    for (i in 1:ncol(x)) {
      graphics::lines(y = x[,i],
                       x = 1:nrow(x),
                       pch = 19,
                       col = i)
    }
  } else {
    # Plot first species
    plot(y = x[,1],
         x = 1:nrow(x),
         type = "l",
         pch = 19,
         col = 1,
         ylim = c(min(x, na.rm = T), max(x, na.rm = T)),
         xlab = "Time",
         ylab = "Abundance")
    title(title, adj = 0, line = 0.5)
    # Add additional points
    for (i in 2:ncol(x)) {
      graphics::lines(y = x[,i],
                       x = 1:nrow(x),
                       pch = 19,
                       col = i)
    }
  }
}

#' Estimate of Taylor's Power Law
#'
#' @param vari Numeric. A vector of variances.
#' @param meani Numeric. A vector of means.
#'
#' @returns A named vector with coefficients a and b for the relation between variance and mean according to Taylor's Power Law.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @noRd
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

#' Report a single warning if there are several
#'
#' @noRd
warn_once <- function(expr) {
  warned <- FALSE # initial state (no warnings yet)
  withCallingHandlers(
    expr, # run whatever function you want to suppress
    warning = function(w) {
      if ( isFALSE(warned)  ) {
        warned <<- TRUE # <<- changes the variable outside warning()
        warning(conditionMessage(w), call. = FALSE) # show warning from function
      }
      # suppress next warnings
      invokeRestart("muffleWarning")
    }
  )
}

#' Short species names for plots 
#'
#' @noRd
short_names <- function(x) {
  # split by " " or "_"
  split_names <- strsplit(x, " |_")
  # put genus in upper case and only first letter
  genera <- sapply(split_names, function(y) {
    # In only genus is available do not shorten
    if( length(y) == 1 ) {
      paste0(toupper(substring(y[1], 1, 1)),
             substring(y[1], 2, nchar(y[1])))
    } else {
      toupper(substring(y[1], 1, 1))
    }
  })
  # paste rest
  others <- sapply(split_names, function(y) paste(y[-1], collapse = " "))
  # paste together
  split_pasted <- paste(genera, others, sep = " ")
  return(split_pasted)
}

#' Identify transient species
#'
#' @noRd
get_transient <- function(x, threshold = 0.3) {
  # Get number of zeros (0) and missing years by species
  missing <- colSums(x == 0 | is.na(x))
  
  missing_n <- data.frame(taxon = names(missing),
                          n_missing = missing,
                          p_missing = missing / nrow(x))
  missing_n$transient <- ifelse(missing_n$p_missing > threshold, "x", "")
  
  rownames(missing_n) <- NULL
  return(missing_n)
}