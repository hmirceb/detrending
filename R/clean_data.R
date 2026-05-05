#' Clean community data in wide format
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param community_col Character. Name of column with the community identifier. Default "comm".
#' @param time_col Character. Name of column with time variable. Default "time".
#' @param na_zero Boolean. Replace missing values (NAs) with zeros (0). Default TRUE.
#' @param filter_transient Boolean. Filter transient species Default TRUE.
#' @param empty_years Boolean. Remove empty years. Default FALSE.
#' @param threshold Numeric. Minimum proportion (between 0 and 1) of time points with valid data to consider a species as transient. Default 0.3.
#'
#' @returns A data.frame with the community data in wide format.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' data(example_data_wide) # Load sample data
#' clean_community_wide(x = example_data_wide,
#'                      community_col = "comm",
#'                      time_col = "time")
#' @export
clean_community_wide <- function(x, 
         community_col = "comm",
         time_col = "time",
         na_zero = TRUE,
         filter_transient = TRUE,
         empty_years = FALSE,
         threshold = 0.3) {
  # Set data as DF just in case its a tibble
  x <- as.data.frame(x)
  
  # Check community column, if not present create one and assume a single community
  if( !community_col %in% colnames(x) ) {
    warning("Missing 'community' column. Data are assumed to belong to a single community.",
            call. = FALSE)
    community_col <- "comm"
    x <- cbind(comm = as.character(rep(1, times = nrow(x))), x)
  }
  # Check time column
  if( !time_col %in% colnames(x) ) {
    stop("Missing 'time' column.", 
         call. = FALSE)
  }
  # Check if any time point has repeated data
  if( sum(duplicated(paste(sep = "_", x[,community_col], x[,time_col]))) != 0 ) {
    stop("Multiple rows for a single community and time point are not allowed.")
  }

  # Set NAs to 0
  if( isTRUE(na_zero) ){
    x[is.na(x)] <- 0
  }
  
  # Split data into communities
  splitting_factor <- as.character(x[, community_col])
  comm_list <- split(x = x, f = splitting_factor)
  
  # Filter data if necessary
  if( isTRUE(filter_transient) ) {
    
    filtered_comms_list <- lapply(comm_list, FUN = function(y) {
      # Get indices of species columns
      sps_index <- !colnames(y) %in% c(community_col, time_col)
      # Get number of zeros (0) and missing years by species
      missing <- colSums(y[, sps_index] == 0 | is.na(y[, sps_index]))/nrow(y[, sps_index])
      
      # Get transient species (missing propoportion higher than threshold)
      transient_sps <- which(missing > threshold)+
        length(c(community_col, time_col)) # this offsets the missing comm and time columns
      
      # If no transient species are detected keep everything
      if( length(transient_sps) == 0 ){ 
        filtered_comm <- y
      } else {
        # Else set their abundance to 0
        y[,transient_sps] <- 0
        filtered_comm <- y
      }
      
      # Reshape to long format to facilitate joining later
      filtered_comm_long <- stats::reshape(data = filtered_comm,
                     direction = "long",
                     varying = which(!names(filtered_comm) %in% c(community_col, time_col)),
                     v.names = "abundance",
                     times = colnames(filtered_comm)[!colnames(filtered_comm) %in% c(community_col, time_col)],
                     timevar = "species")
      return(filtered_comm_long)
    }
    )
    
    # Join community data into single matrix 
    data <- do.call("rbind", filtered_comms_list)
    # Add ID (comm by year) col
    data$id_comm <- paste(sep = "_",
                    unlist(lapply(strsplit(rownames(data), "\\."), function(x) x[1])),
                    unlist(lapply(strsplit(rownames(data), "\\."), function(x) x[2])))
    # Add species col
    data$species <- unlist(lapply(strsplit(rownames(data), "\\."), function(x) x[length(x)]))
    
    # Pivot data to years by species matrix
    d_wide <- stats::reshape(data,
                      direction = "wide",
                      timevar = "species",
                      idvar = "id_comm",
                      v.names = "abundance")
    colnames(d_wide) <- gsub("abundance.", "", colnames(d_wide))
    rownames(d_wide) <- NULL
  } else {
    d_wide <- x
  }
  
  # Remove species with 0 abundance after filtering by threshold
  id_cols <- colnames(d_wide) %in% c(community_col, time_col, "id", "id_comm")
  sps_to_remove <- colSums(d_wide[, !id_cols]) == 0
  sps_to_remove <- names(sps_to_remove)[sps_to_remove]
  d_wide <- d_wide[, !colnames(d_wide) %in% sps_to_remove]
  
  # Remove empty years if necessary
  if( isFALSE(empty_years) ) {
    # Get years with total abundance = 0
    year_with_data <- rowSums(d_wide[which(!colnames(d_wide) %in% c(community_col, time_col, "id_comm"))], na.rm = TRUE) > 0
    # Remove them
    d_wide <- d_wide[year_with_data,]
  }
  
  # Sort by community and time
  d_wide <- d_wide[with(d_wide, order(d_wide[, community_col], d_wide[, time_col])),]
  # Remove id columns generated in the process
  d_wide <- d_wide[, !colnames(d_wide) %in% c("id", "id_comm")]
  
  return(d_wide)
}

#' Clean community data in long format
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param community_col Character. Name of column with the community identifier. Default "comm".
#' @param time_col Character. Name of column with time variable. Default "time".
#' @param taxa_col Character. Name of column with taxa names. Default "species".
#' @param abundance_col Character. Name of column with abundance values. Default "abundance".
#' 
#' @returns A data.frame with the community data in long format.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' data(example_data_long) # Load sample data
#' dim(example_data_long) # 300 4
#' data_wide <- clean_community_long(x = example_data_long,
#'                      community_col = "comm",
#'                      time_col = "time",
#'                      taxa_col = "species",
#'                      abundance_col = "abundance")
#' dim(data_wide) # 30 12
#' @export
clean_community_long <- function(x, 
                                 community_col = "comm",
                                 time_col = "time",
                                 taxa_col = "species",
                                 abundance_col = "abundance") {
  # Set data as DF just in case its a tibble
  x <- as.data.frame(x)
  
  # Check community column, if not create one and assume a single community
  if( !community_col %in% colnames(x) ) {
    warning("Missing 'community' column. Data are assumed to belong to a single community.")
    community_col <- "comm"
    x <- cbind(comm = as.character(rep(1, times = nrow(x))), x)
  }
  # Check time column
  if( !time_col %in% colnames(x) ) {
    stop("Missing 'time' column.")
  }
  # Check species column
  if( !taxa_col %in% colnames(x) ) {
    stop("Missing 'species' column.")
  }
  # Check abundance column
  if( !abundance_col %in% colnames(x) ) {
    stop("Missing 'abundance' column.")
  }
  # Check if any time point has repeated data
  if( any(duplicated(paste(sep = "_", x[,community_col], x[,time_col], x[,taxa_col]))) ) {
    stop("Multiple rows for a single time point are not allowed.")
  }
  
  # Select only columns of interest
  x <- x[, colnames(x) %in% c(community_col, time_col, taxa_col, abundance_col)]
  
  # Set NAs to 0 to avoid pivoting problems
  x[, abundance_col][is.na(x[, abundance_col])] <- 0
  
  # Make ID col for pivoting
  x <- cbind(
    id_comm = paste(sep = "_", 
                    x[, community_col], 
                    x[, time_col]), 
    x)
  x$id_comm <- as.factor(x$id_comm)
  
  # Remove rows with NAs in taxa_col or time_col 
  if( any( is.na(x[, taxa_col])) ){
    mis_tax <- is.na(x[, taxa_col])
    message(paste0(sum(mis_tax), " row(s) with missing values in column '", taxa_col, "' removed."))
    x <- x[!mis_tax,]
  }
  if( any( is.na(x[, time_col])) ){
    mis_time <- is.na(x[, time_col])
    message(paste0(sum(mis_time), " row(s) with missing values in column '", time_col, "' and were removed."))
    x <- x[!mis_time,]
  }
  
  # Pivot data to years by species matrix
  data_wide <- stats::reshape(x,
                       direction = "wide",
                       idvar = "id_comm",
                       timevar = taxa_col,
                       v.names = abundance_col)
  colnames(data_wide) <- gsub(paste0(abundance_col, "."), "", colnames(data_wide))
  rownames(data_wide) <- NULL
  
  # Sort by community and time
  data_wide <- data_wide[with(data_wide, order(data_wide[, community_col], data_wide[, time_col])),]
  # Remove id_comm column
  data_wide <- data_wide[, !colnames(data_wide) == "id_comm"]

  return(data_wide)
}

#' Clean community data to use in other functions
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param input_format Character. Format of the data to clean. One of "wide" or "long".
#' @param community_col Character. Name of column with the community identifier. Default "comm".
#' @param time_col Character. Name of column with time variable. Default "time".
#' @param taxa_col Character. Name of column with taxa names. Default "species".
#' @param abundance_col Character. Name of column with abundance values. Default "abundance".
#' @param na_zero Boolean. Replace missing values (NAs) with zeros (0). Default TRUE.
#' @param filter_transient Boolean. Filter transient species. Default FALSE.
#' @param empty_years Boolean. Remove empty years. Default FALSE.
#' @param threshold Numeric. Minimum proportion (between 0 and 1) of time points with valid data to consider a species as transient. Default 0.3.
#'
#' @returns A data.frame with community data ready to use in other functions.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @examples
#' require(detrending)
#' 
#' # Clean data in long format
#' data(example_data_long)
#' clean_community(x = example_data_long, 
#'                 input_format = "long",
#'                 community_col = "comm",
#'                 time_col = "time",
#'                 taxa_col = "species",
#'                 abundance_col = "abundance")
#' 
#' # Clean data in wide format
#' data(example_data_wide)
#' clean_community(x = example_data_wide, 
#'                 input_format = "wide",
#'                 community_col = "comm",
#'                 time_col = "time")
#' @export
clean_community <- function(x, 
                            input_format = "wide",
                            community_col = "comm",
                            time_col = "time",
                            taxa_col = "species",
                            abundance_col = "abundance",
                            na_zero = TRUE,
                            filter_transient = FALSE,
                            empty_years = FALSE,
                            threshold = 0.3) {
  # Check input format
  if( !input_format %in% c("long", "wide") ){
    stop("Please provide a valid input format.")
  }
  
  # Choose method for "wide" format
  if( input_format == "wide" ) {
    clean_data <- clean_community_wide(x = x,
                         community_col = community_col,
                         time_col = time_col,
                         na_zero = na_zero,
                         filter_transient = filter_transient,
                         empty_years = empty_years,
                         threshold = threshold)
  }
  # Method for "long" format
  if ( input_format == "long" ) {
    # Process "long" data
    ccl <- clean_community_long(x = x,
                         community_col = community_col,
                         time_col = time_col,
                         taxa_col = taxa_col,
                         abundance_col = abundance_col)
    
    # Transform to "wide" format
    clean_data <- clean_community_wide(x = ccl,
                                       community_col = community_col,
                                       time_col = time_col,
                                       na_zero = na_zero,
                                       filter_transient = filter_transient,
                                       empty_years = empty_years,
                                       threshold = threshold)
  }
  
  return(clean_data)
}

#' Prepare metacommunity data from a wide format data.frame
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param community_col Character. Name of column with the community identifier. Default "comm".
#' @param time_col Character. Name of column with time variable. Default "time".
#' @param taxa_col Character. Name of column with taxa names. Default "species".
#'
#' @returns A data.frame with the community data ready to use by the `metacomstab_term()` function.
#' 
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @noRd
metacoms_data <- function(x, 
                         community_col = "comm",
                         time_col = "time",
                         taxa_col = "species"
) {
  # Get community and time columns
  ids <- colnames(x) %in% c(community_col, time_col)
  
  # Pivot to long format
  data_long <- stats::reshape(
    data = x,
    direction = "long",
    varying = colnames(x[,!ids]),
    v.names = "value",
    timevar = taxa_col,
    times = colnames(x[,!ids]),
    idvar = c(community_col, time_col)
  )
  
  # Pivot to wide format
  data_wide <- stats::reshape(
    data = data_long,
    direction = "wide",
    idvar = c(community_col, taxa_col),
    timevar = time_col
  )
  
  # Clean column and row names
  colnames(data_wide) = gsub("value.", "t", colnames(data_wide))
  rownames(data_wide) = NULL
  
  # Set NAs to 0
  data_wide[is.na(data_wide)] <- 0
  
  return(data_wide)
}