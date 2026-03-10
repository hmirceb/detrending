#' Species richness
#'
#' @param x Numeric. A vector of abundances.
#'
#' @returns A numeric value.
#' 
#' @noRd
richness <- function(x){
  sum(x > 0)
}

#' Shannon-Wiener index
#'
#' @param x Numeric. A vector of abundances.
#' @param relative Boolean. If abundance is expressed relative to a total instead of counts or biomass. Default FALSE.
#' 
#' @returns A numeric value.
#' 
#' @noRd
shannon <- function(x, relative = FALSE) {
  a <- x[x != 0]
  if ( isFALSE(relative)) {
    b <- a/sum(a)
  } else {
    b <- a 
  }
  H <- -1*sum(b*log(b))
  return(H)
}

#' Shannon-Wiener index
#'
#' @param x Numeric. A vector of abundances.
#' @param relative Boolean. If abundance is expressed relative to a total instead of counts or biomass. Default FALSE.
#' 
#' @returns A numeric value.
#' 
#' @noRd
simpson <- function(x, relative = FALSE) {
  a <- x[x != 0]
  if ( isFALSE(relative) ) {
    b <- a/sum(a)
  } else {
    b <- a 
  }
  D <- sum(b*b)
  return(D)
}

#' Pielou's evenness index
#'
#' @param x Numeric. A vector of abundances.
#'
#' @returns A numeric value.
#' 
#' @noRd
pielou <- function(x) {
  J <- shannon(x)/log(richness(x))
  return(J)
}

#' Explore community data
#'
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns.
#' @param by_timestep Boolean. Get community information at each time step. Default FALSE.
#' @param total Character. Wether to compute diversity indices from the average relative abundance of species across years (overall) or the average of annula diversity indices.  
#' @param community_col Character. Name of column with the community identifier.
#' @param time_col Character. Name of column with time variable.
#' @param trend Character. Method to check for trends in species abundance
#'
#' @returns A data.frame wiht.
#' 
#' @export
comm_expl <- function(x,
                     by_timestep = FALSE,
                     total = c("average", "overall"),
                     community_col = "comm",
                     time_col = "time",
                     trend = c("none", "dennis", "loglinear")){
  
  # Check community column, if not present create one and assume a single community
  if( !community_col %in% colnames(x) ) {
    warning("Missing 'community' column. Data are assumed to belong to a single community.",
            call. = FALSE)
    community_col <- "comm"
    x <- cbind(comm = as.character(rep(1, times = nrow(x))), x)
  }
  
  # Split data by timestep
  if( isFALSE(by_timestep) ){
    splitting_factor <- as.character(x[, community_col])
  } else {
    splitting_factor <- paste(sep = "_", x[, community_col], x[, time_col])
  }
  # split data by community
  c_list <- split(x, f = splitting_factor)
  
  # Get info per community
  info_by_comm <- lapply(c_list, function(c_com){
    # Get columns with community ids
    sps_index <- !colnames(c_com) %in% c(community_col, time_col)
    # Number of years
    ny <- nrow(c_com[,sps_index])
    
    # Calculate diversity from average annual diversity values 
    if( total == "average" ) {
      # Average species richness
      S <- mean(apply(c_com[,sps_index], MARGIN = 1, FUN = richness))
      # Shannon's index
      H <- mean(apply(c_com[,sps_index], MARGIN = 1, FUN = shannon))
      # Pielou`s evenness
      J <- mean(apply(c_com[,sps_index], MARGIN = 1, FUN = pielou))
    }
    
    # Calculate diversity from average relative abundance of each species over time
    if( total == "overall" ) {
      # Average relative abundance per species
      b <- colMeans(
        c_com[,sps_index] / rowSums(c_com[,sps_index])
      )
      b <- b[b != 0] # remove 0s
      # Average species richness
      S <- richness(b)
      # Shannon's index from average relative abundances
      H <- -1*sum(b*log(b))
      # Pielou`s evenness
      J <- pielou(colMeans(c_com[,sps_index]))
    }
    
    # Format
    if( isFALSE(by_timestep) ){
      res <- data.frame(comm = unique(c_com[, community_col]),
                       ny = ny, 
                       S = S, 
                       H = H, 
                       J = J)
      # Reorder result
      res <- res[with(res, order(res[, community_col])),]
    } else { # Information by timestep
      res <- data.frame(comm = unique(c_com[, community_col]),
                       time = as.numeric(unique(c_com[, time_col])),
                       S = S,
                       H = H, 
                       J = J)
      # Reorder result
      res <- res[with(res, order(res[, community_col], 
                                 res[, time_col])),]
      
    }
    return(res)
    }
    )
  # Join community-wise results into single df
  info_df <- do.call("rbind", info_by_comm)
  # Remove row names
  rownames(info_df) <- NULL
  
  return(info_df)
  }
