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
#' @param total Character. Whether to compute diversity indices from the average relative abundance of species across years (overall) or the average of annual diversity indices.  
#' @param community_col Character. Name of column with the community identifier.
#' @param time_col Character. Name of column with time variable.
#' @param trends Boolean. Check for trends in species using linear regression on log-transformed abundances. Default FALSE.
#' @param check_dominants Boolean. Check if dominant species according to threshold *q* have missing data. Default FALSE.
#' @param q Numeric. A number between 0 and 1 indicating the abundance threshold to consider a species as dominant.
#'  
#' @returns A named list:
#'  - `diversity`: A data.frame indicaitng the number of timesteps in each community, along with their species richness (S), Shannon's (H) and Pielou's indices (P).
#'  - `trends`: A data.frame with the estimated mean abundance trends of the species in each community as returned by `comm_trend()`.
#'  - `dominant_taxa`: A data.frame with the species considered dominant in each species and the number of missing data points in the time series.
#'  
#' @export
comm_expl <- function(x,
                       by_timestep = FALSE,
                       total = "average",
                       community_col = "comm",
                       time_col = "time",
                       trends = FALSE,
                       check_dominants = FALSE,
                       q = 0.7){
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
  
  # split data by community
  c_list <- split(x, f = as.character(x[, community_col]))
  
  # Estimate trends (comm_trend already checks the time column)
  if ( isTRUE(trends) ) {
    trends_df <- lapply(c_list, function(t_com){
      sps_index <- !colnames(t_com) %in% c(community_col, time_col)
      cbind(comm = unique(t_com[, community_col]),
            comm_trend(x = t_com[sps_index], 
                       method = "loglinear", 
                       plot = F))
    }
    )
    trends_df <- do.call("rbind", trends_df)
    colnames(trends_df)[1] <- community_col
    rownames(trends_df) <- NULL
  }
  
  # Check dominants with missing values
  if (isTRUE(check_dominants) ) {
    dom_check <-lapply(c_list, function(c_com) {
      # Check time column or add one if missing
      suppressWarnings(
        c_com <- check_time(c_com, time_col = time_col, term = "two", rm = FALSE)
      )
      # Remove species with no abundance
      c_com <- remove_empty_sps(c_com, time_col = time_col, community_col = community_col)
      # Get columns with species ids
      sps_index <- !colnames(c_com) %in% c(community_col, time_col)
      
      # Check if dominant species have missing years
      suppressMessages(
        dom_check <- check_dominants(c_com[sps_index], q = q)
      )
      dom_check <- cbind(comm = unique(c_com[, community_col]),
                           dom_check)
      return(dom_check)
    }
    )
    dom_check <- do.call("rbind", dom_check)
    rownames(dom_check) <- NULL
    colnames(dom_check)[1] <- community_col
    
    if( nrow(dom_check) == 0){
      dom_check <- "No dominant species with missing values."
    }
  }
  
  # Get info per community
  info_by_comm <- warn_once( # avoid a warning for each comm if time_col is missing
    lapply(c_list, function(c_com){
      # Check time column or add one if missing
      c_com <- check_time(c_com, time_col = time_col, term = "two", rm = FALSE)
      # Remove species with no abundance
      c_com <- remove_empty_sps(c_com, time_col = time_col, community_col = community_col)
      
      # Get columns with species ids
      sps_index <- !colnames(c_com) %in% c(community_col, time_col)
      # Number of years
      nt <- nrow(c_com[,sps_index])
      
      if ( isTRUE(by_timestep) ) {
        # Species richness
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
        colnames(res)[1:2] <- c(community_col, time_col)
        
      } else {
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
                          nt = nt, 
                          S = S, 
                          H = H, 
                          J = J)
        colnames(res)[1] <- community_col
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
    if ( isTRUE(check_dominants) ) {
      res_list <- list(diversity = info_df, trends = trends_df, dominant_taxa = dom_check)
    } else {
      res_list <- list(diversity = info_df, trends = trends_df) 
    }
  } else {
    if ( isTRUE(check_dominants) ) {
      res_list <- list(diversity = info_df, dominant_taxa = dom_check)
    } else {
      res_list <- list(diversity = info_df) 
    }
  }
  
  return(res_list)
}
