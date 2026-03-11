#' Species richness
#'
#' @param x Numeric. A vector of abundances.
#'
#' @returns A numeric value.
#' 
#' @noRd
richness <- function(x){
  sum(x > 0, na.rm = T)
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
  a <- x[x != 0 & !is.na(x)]
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
  a <- x[x != 0 & !is.na(x)]
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
#' @param trends Boolean. Check for trends in species using linear regression on log-transformed abundances. Default FALSE.
#'
#' @returns A named list:
#'  - `diversity`: A data.frame with several diversity metrics for each community.
#'  - `trends`: A data.frame with the estimated mean abundance trends of the species in each community.
#' @export
comm_expl <- function(x,
                      by_timestep = FALSE,
                      total = "average",
                      community_col = "comm",
                      time_col = "time",
                      trends = FALSE){
  # Check arguments
  if( !total %in% c("average", "overall") ) {
    stop("Argument 'total' must be one of 'average' or 'overall'")
  }
  
  # Check community column, if not present create one and assume a single community
  if( !community_col %in% colnames(x) ) {
    warning("Missing community column. Data are assumed to belong to a single community.",
            call. = FALSE)
    community_col <- "comm"
    x <- cbind(comm = as.character(rep(1, times = nrow(x))), x)
  }
  
  # split data by community (or community and time step)
  c_list <- split(x, f = as.character(x[, community_col]))
  
  if ( isTRUE(trends) ) {
    # Estimate trends (comm_trend already checks the time column)
    trends_df <- lapply(c_list, function(t_com){
      sps_index <- !colnames(t_com) %in% c(community_col, time_col)
      cbind(comm = unique(t_com[, community_col]),
            comm_trend(x = t_com[sps_index], 
                       method = "loglinear", 
                       plot = F))
    }
    )
    trends_df <- do.call("rbind", trends_df)
    rownames(trends_df) <- NULL
  }
  
  # Get info per community
  info_by_comm <- warn_once( # avoid a warnring for each comm if time_col is missing
    lapply(c_list, function(c_com){
      if ( isTRUE(by_timestep) ) {
        # Check time column or add one if missing
        c_com <- check_time(c_com, time_col = time_col, term = "two", rm = FALSE)
        c_com <- remove_empty_sps(c_com, time_col = time_col, community_col = community_col)
        # Get columns with community ids
        sps_index <- !colnames(c_com) %in% c(community_col, time_col)
        # Number of years
        ny <- nrow(c_com[,sps_index])
        # Average species richness
        S <- apply(c_com[,sps_index], MARGIN = 1, FUN = richness)
        # Shannon's index
        H <- apply(c_com[,sps_index], MARGIN = 1, FUN = shannon)
        # Pielou's evenness
        J <- apply(c_com[,sps_index], MARGIN = 1, FUN = pielou)
        # Format results into df
        res <- data.frame(comm = unique(c_com[, community_col]),
                          time = unique(c_com[, time_col]),
                          S = S,
                          H = H, 
                          J = J)
      } else {
        # remove columns with empty species (just in case)
        c_com <- remove_empty_sps(c_com, time_col = time_col, community_col = community_col)
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
            c_com[,sps_index] / rowSums(c_com[,sps_index], na.rm = TRUE),
            na.rm = TRUE
          )
          b <- b[b != 0] # remove 0s
          # Species richness
          S <- richness(b)
          # Shannon's index from average relative abundances
          H <- -1*sum(b*log(b))
          # Pielou's evenness
          J <- H / log(S)
        }
        # Format results into df
        res <- data.frame(comm = unique(c_com[, community_col]),
                          ny = ny, 
                          S = S, 
                          H = H, 
                          J = J)
      }
      return(res)
    }
    )
  )
  # Join community-wise results into single df
  info_df <- do.call("rbind", info_by_comm)
  # Remove row names
  rownames(info_df) <- NULL
  
  # Reorder result
  if ( isTRUE(by_timestep) ) {
    info_df <- info_df[with(info_df, order(info_df[, community_col], 
                                           info_df[, time_col])),]
  } else {
    info_df <- info_df[with(info_df, order(info_df[, community_col])),]
  }
  
  # Create resulting list 
  if ( isTRUE(trends) ) {
    res_list <- list(diversity = info_df, trends = trends_df)
  } else {
    res_list <- list(diversity = info_df)
  }
  
  return(res_list)
}
