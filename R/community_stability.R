#' Estimate stability and sychrony metrics for one or more communities
#'
#' @param x A data.frame. A community matrix of species abundance with years as rows and species as columns.
#' @param time_col Character. Name of column with time variable.
#' @param community_col Character. Name of column with the community identifier.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param cv Character. Coefficient of variation, one or both of "com" (for the CV of the sum of species' abundances) or "pop" (average of species' CVs).
#' @param weighted Boolean. Weight population CV and stability metrics by species mean abundance. Default TRUE.
#' @param synchrony Boolean. Estimate synchrony metrics. Default TRUE.
#' @param index Character. Synchrony index to calculate. One of "psi", (Segrestin *et al.* 2024), "phi" (Loreau & Mazancourt 2008), "eta" (Gross *et al.* 2014) or "logvar" (Leps *et al.* 2018). 
#'
#' @returns A data frame with the selected metrics for each community, their value, the variance term used to estimate them and whether weighting was used.
#'
#' @examples
#' #' require(detrending)
#' 
#' # Simulate community data
#' comm_df <- sim_mvcomm()
#' community_stability(x = comm_df$sim_data)
#' 
#' @export
community_stability <- function(
    x, 
    time_col = "time",
    community_col = "comm",
    term = "var",
    cv = c("com", "pop"),
    weighted = TRUE,
    synchrony = TRUE,
    index = c("psi", "phi", "eta", "logvar")
    ) {
  
  # Check community column, if not present create one and assume a single community
  if( !community_col %in% colnames(x) ) {
    warning("Missing community column. Data are assumed to belong to a single community.",
            call. = FALSE)
    community_col <- "comm"
    x <- cbind(comm = as.character(rep(1, times = nrow(x))), x)
  }
  
  # split data by community
  c_list <- split(x, f = as.character(x[, community_col]))
  
  res_list <- lapply(c_list, function(y) {
    # Check if a time column was specified for detrending methods and order rows
    y <- check_time(y, time_col = time_col, term = term, rm = FALSE)
    # Replace NAs with 0 and remove columns (species) with 0 abundance across all years 
    y <- remove_empty_sps(x = y, time_col = time_col, community_col = community_col)
    # get community id
    comm_id <- unique(y[,colnames(y) %in% c(community_col)])
    # remove community column
    y <- y[,!colnames(y) %in% c(community_col)]
    
    # synchrony metrics
    if ( isTRUE(synchrony) ){
      synchrony <- sync_term(y, 
                             index = index, 
                             term = term, 
                             time_col = time_col, 
                             weighted = weighted)
    }
    
    # cv metrics
    if ( any(cv %in% "com") ){
      cv_com <- cv_com_term(y, total = TRUE, 
                        time_col = time_col, 
                        term = term)
    } else { cv_com <- NA }
    
    if ( any(cv %in% "pop") ){
      cv_pop <- cv_com_term(y, total = F, 
                        time_col = time_col, 
                        term = term, 
                        weighted = weighted)$CV
    } else { cv_pop <- NA }
    
    res <- c(synchrony, cv_com, cv_pop)
    res_df <- data.frame(metric = names(res), 
                         value = res,
                         term = term)
    res_df$metric <- gsub(paste0("_", term), "", res_df$metric)
    res_df <- cbind(comm = comm_id,
                    res_df)
    return(res_df)
  }
  )
  comm_stability <- do.call("rbind", res_list)
  comm_stability <- comm_stability[comm_stability$metric != "", ]
  rownames(comm_stability) <- NULL
  
  # tag weighted metrics
  comm_stability$weighted <- isTRUE(weighted) & grepl("w", comm_stability$metric)
  # fix some metric names 
  comm_stability$metric <- gsub("w", "", comm_stability$metric)
  comm_stability$metric <- ifelse(comm_stability$metric == "CVt", "CV_com", comm_stability$metric)
  
  return(comm_stability)
}
