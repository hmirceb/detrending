# Functions to decompose metacommunity stability

#' Partitioning of metacommunity stability in its components
#' 
#' This function partitions the stability of a metacommunity following the methods by Segrestin & Leps (2022) and Hammond *et al.* (2020). It allows standard estimates of variance and CV as well as dentrended versions using Hill's two and three term local quadratic variance estimates. Ideally input data should be processed with `clean_community()` before use. 
#'
#' @usage cv2_decomp(x, community_col = "comm", 
#' time_col = "time", taxa_col = "species", term = "var", nrand = NA)
#'
#' @param x A data.frame. Metacommunity matrix with time in rows, taxa in columns and a column identifying each community.
#' @param community_col Character. Name of column with the community identifier.
#' @param time_col Character. Name of column with time variable.
#' @param taxa_col Character. Name of column with taxa names. Default "species".
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param nrand Numeric. Number of randomizations to estimate confidence intervals. Default NA.
#'
#' @returns An object of class `cv.dec`, a list of named vectors.
#' 
#' @references
#' - Segrestin, J., & Lepš, J. (2022). Towards a better ecological understanding of metacommunity stability: A multiscale framework to disentangle population variability and synchrony effects. Journal of Ecology, 110(7), 1632-1645.
#' - Hammond, M., Loreau, M., De Mazancourt, C., & Kolasa, J. (2020). Disentangling local, metapopulation, and cross‐community sources of stabilization and asynchrony in metacommunities. Ecosphere, 11(4), e03078.
#'
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @export
cv2_decomp_term <- function(x, 
                            community_col = "comm",
                            time_col = "time",
                            taxa_col = "species",
                            term = "var",
                            nrand = NA){
  # Set data as DF just in case its a tibble
  x <- as.data.frame(x)
  
  # Match argument for variance function to use
  # Table of options
  options <- data.frame(term = c("var", "two", "three", "linear"),
                        var = c("var", "var_t2", "var_t3", "var_linear"))
  # Get choice
  var_func <- match.fun(options[options$term == term,]$var)
  
  #-----------------#
  #### Data prep #### 
  #-----------------#
  x <- metacoms_data(x = x,
                       community_col = community_col,
                       time_col = time_col,
                       taxa_col = taxa_col)
  
  # Check column format
  if( !taxa_col %in% colnames(x) ) {
    stop(
      paste0("Column ", taxa_col, " not found in data")
    )
  }
  if( !community_col %in% colnames(x) ) {
    stop(
      paste0("Column ", community_col, " not found in data")
    )
  }
  if( length(unique(x[, community_col])) == 1 ) {
    message(
      "Only 1 community provided"
    )
  }
  # Set communty and species as factor
  x[, community_col] <- as.factor(x[, community_col])
  x[, taxa_col] <- as.factor(x[, taxa_col])
  
  # Matrix with only biomass (no species or community column)
  comtaxa_ids <- colnames(x) %in% c(taxa_col, community_col)
  X <- x[, !comtaxa_ids]
  
  # Check species with 0 biomass across years and remove them
  n <- sum(rowSums(X) == 0)
  if(n > 0) message(paste(n, "row(s) with no biomass were ommited"))
  x <- x[rowSums(X) != 0, ]
  X <- x[, !comtaxa_ids]
  
  #---------------------#
  #### Metacommunity #### 
  #---------------------#
  n <- nrow(x) # Number of species and community combinations
  
  var_metacom <- var_func(colSums(X)) # Variance of the meta-community abundance
  mean2_metacom <- mean(colSums(X))^2 # Squared mean of meta-community abundance
  CV2_metacom <- var_metacom / mean2_metacom # CV squared of metacommunity abundance
  
  #-------------------#
  #### Populations ####
  #-------------------#
  mat_all <- vcov_term(t(X), term = term) # VCOV matrix (standard, 2 or 3 three term) 
  var_all <- diag(mat_all) # Variance of each population
  pop_var <- sum(var_all) / mean2_metacom # CV squared of all populations if they are independent
  diag(mat_all) <- 0 # Set diagonal (variances) to 0
  
  # DF with pairwise covariances between populations (species by community) and max CV 
  cov_all <- data.frame(com1 = rep(x[, community_col], n), sp1 = rep(x[, taxa_col], n), # Combination of community and species 1
                        com2 = rep(x[, community_col], each = n), sp2 = rep(x[, taxa_col], each = n), # Combination of community and species 2
                        cov = as.vector(mat_all) / mean2_metacom, # eq.4 in paper, coefficient of covariation between all combinations of populations
                        max = as.vector(outer(sqrt(var_all), # eq. 2 in paper, coefficient of covariation between 
                                              sqrt(var_all))) / mean2_metacom)
  
  # Unique ID for each population (combination of species and community)
  real_sp <- paste0(x[, taxa_col], x[, community_col])
  
  # Occurrence of each species in each communities (can be FALSE and FALSE)
  cov_all$in1 <- paste0(cov_all$sp1, cov_all$com2) %in% real_sp
  cov_all$in2 <- paste0(cov_all$sp2, cov_all$com1) %in% real_sp
  
  # Within species synchrony (Pop.sync_intra in paper): 
  # covariance terms of populations of the same species in different communities
  within <- cov_all[cov_all$sp1 == cov_all$sp2, ] # Filter 
  within <- c(value = sum(within$cov),
              max = sum(within$max),
              Beta_MP = sum(within$max) - sum(within$cov), # Hammond et al. 2020
              n = nrow(within))
  
  # Direct interaction synchrony (Pop.sync_direct): 
  # covariance terms of populations of different species within the same community
  direct <- cov_all[cov_all$com1 == cov_all$com2 & cov_all$sp1 != cov_all$sp2, ]
  direct <- c(value = sum(direct$cov),
              max = sum(direct$max),
              delta = sum(direct$max) - sum(direct$cov), # Hammond et al. 2020
              n = nrow(direct))
  
  # Indirect interaction synchrony (Pop.sync_indirect): 
  # covariance terms of populations of different species in different communities but where both species cooccur in at least one community
  indirect <- cov_all[cov_all$com1 != cov_all$com2 
                      & cov_all$sp1 != cov_all$sp2
                      & (cov_all$in1 | cov_all$in2), ]
  indirect <- c(value = sum(indirect$cov),
                max = sum(indirect$max),
                Beta_CCi = sum(indirect$max) - sum(indirect$cov), # Hammond et al. 2020
                n = nrow(indirect))
  
  # No interaction synchrony (Pop.sync_no): 
  # covariance terms of populations of different species in different communities but where both species occur in their community only
  no <- cov_all[cov_all$com1 != cov_all$com2 
                & cov_all$sp1 != cov_all$sp2
                & !(cov_all$in1 | cov_all$in2),]
  no <- c(value = sum(no$cov),
          max = sum(no$max),
          Beta_CCno = sum(no$max) - sum(no$cov), # Hammond et al. 2020
          n = nrow(no))
  
  # Add results to list
  res <- list(CV2 = CV2_metacom,
              pop_var = pop_var,
              Pop_sync_direct = direct,
              Pop_sync_intra = within,
              Pop_sync_indirect = indirect,
              Pop_sync_no = no)
  
  #--------------------------#
  #### Randomization test ####
  #--------------------------#
  
  if ( is.numeric(nrand) ){
    rand <- matrix(NA, nrand, 4)
    for(i in 1:nrand){
      randomize <- apply(X, MARGIN = 1, sample, size = ncol(X))
      mat_all <- vcov_term(randomize, term = term)
      diag(mat_all) <- 0
      
      cov_rand <- as.vector(mat_all) / mean2_metacom
      within_rand <- sum(cov_rand[cov_all$sp1 == cov_all$sp2])
      direct_rand <- sum(cov_rand[cov_all$com1 == cov_all$com2
                                  & cov_all$sp1 != cov_all$sp2])
      indirect_rand <- sum(cov_rand[cov_all$com1 != cov_all$com2 
                                    & cov_all$sp1 != cov_all$sp2
                                    & (cov_all$in1 | cov_all$in2)])
      no_rand <- sum(cov_rand[cov_all$com1 != cov_all$com2 
                              & cov_all$sp1 != cov_all$sp2
                              & !(cov_all$in1 | cov_all$in2)])
      rand[i, ] <- c(direct_rand, within_rand, indirect_rand, no_rand)
    }
    res <- c(res, rand = list(rand))
  }
  
  class(res) <- "cv.dec"
  return(res)                    
}

#' Print method for cv2_decomp_term()
#'
#' @param x Result from cv2_decomp_term()
#' @param ... 
#' 
#' @export
print.cv.dec <- function (x, ...) {
  cat("\nDecomposition of the metacommunity squared coefficient of variation")
  cat("\nSee Segrestin & Leps (2022)")
  cat("\n")
  pop_sync <- x$Pop_sync_intra[1] + x$Pop_sync_direct[1] + x$Pop_sync_indirect[1] + x$Pop_sync_no[1]
  if ("rand" %in% names(x)) {
    pop_sync_rand <- quantile(rowSums(x$rand), probs = c(0.025, 0.975))
    pop_sync_rand <- paste0("[", round(pop_sync_rand[1], 4), "; ", round(pop_sync_rand[2], 4), "]")
    cat(paste0("\nCV2 = ", round(x$CV2, 4), ", Pop.var = ", round(x$pop_var, 4), ", Pop.sync = ", round(pop_sync, 4), " ", pop_sync_rand))
  } else {
    cat(paste0("\nCV2 = ", round(x$CV2, 4), ", Pop.var = ", round(x$pop_var, 4), ", Pop.sync = ", round(pop_sync, 4)))
  }
  
  cat("\n")
  names_pop_sync <- c("Pop.sync[direct]", "Pop.sync[intra]", "pop.sync[indirect]", "pop.sync[no]")
  pop_sync_val <- round(unlist(lapply(x[3:6], "[", 1)), 4)
  names_hamm <- c("Delta", "Beta[MP]", "Beta[CCi]", "Beta[CCno]")
  pop_sync_hamm <- round(unlist(lapply(x[3:6], "[", 3)), 4)
  
  if ("rand" %in% names(x)){
    pop_sync_ind <- apply(x$rand, MARGIN = 2, quantile, probs = c(0.025, 0.975))
    pop_sync_ind <- apply(pop_sync_ind, MARGIN = 2, 
                          function(x) paste0("[", round(x[1], 4), "; ", round(x[2], 4), "]"))
    df <- data.frame(paste0(names_pop_sync, " = ", pop_sync_val),
                     pop_sync_ind,
                     paste0(names_hamm, " = ", pop_sync_hamm))
    colnames(df) <- c("Segrestin & Leps (2022)", "Rand (95% CI)", "Hammond et al. (2020)")
  } else {
    df <- data.frame(paste0(names_pop_sync, " = ", pop_sync_val),
                     paste0(names_hamm, " = ", pop_sync_hamm))
    colnames(df) <- c("Segrestin & Leps (2022)", "Hammond et al. (2020)")
  }
  
  cat("\n")
  print(df, row.names = F, right = F)
  cat("\n")
}

#' Decompose metacommunity stability in its components
#' 
#' This function partitions the stability of a metacommunity following the methods by Segrestin & Leps (2022) and Hammond *et al.* (2020). It allows standard estimates of variance and CV as well as dentrended versions using Hill's two and three term local quadratic variance estimates. Ideally input data should be processed with `clean_community()` before use. A wrapper for `cv2_decomp_term()` that returns a data.frame instead.
#'
#' @param x A data.frame. Community matrix with time in rows and taxa in columns.
#' @param community_col Character. Name of column with the community identifier.
#' @param time_col Character. Name of column with time variable.
#' @param term Character. Term to estimate the variance. One of "var" (for standard variance and covariance), "two" or "three" for Hills' two or three term local quadrat variance and covariance. Default "var".
#' @param method Character. Method used to decompose metacommunity stability. One or both of "segrestin" or "hammond".
#' @param nrand Numeric. Number of randomizations to estimate confidence intervals. Default NA.
#' @param conf Numeric. Confidence level to estimate confidence intervals. Default 0.95.
#'
#' @returns A data.frame with one row per metacommunity and one column per stability component.
#'
#' @references
#' - Segrestin, J., & Lepš, J. (2022). Towards a better ecological understanding of metacommunity stability: A multiscale framework to disentangle population variability and synchrony effects. Journal of Ecology, 110(7), 1632-1645.
#' - Hammond, M., Loreau, M., De Mazancourt, C., & Kolasa, J. (2020). Disentangling local, metapopulation, and cross‐community sources of stabilization and asynchrony in metacommunities. Ecosphere, 11(4), e03078.
#'
#' @author Jules Segrestin, \email{jsegrestin@@gmail.com}
#' @author Héctor Miranda-Cebrián, \email{hectorm94@@gmail.com}
#' 
#' @export
metacomstab_term <- function(x, 
                            community_col = "comm",
                            time_col = "time",
                            term = "var",
                            method = c("segrestin", "hammond"),
                            nrand = NA,
                            conf = 0.95) {
  # Run metacommunity decomposition function
  res <- cv2_decomp_term(x, 
                        community_col = community_col,
                        time_col = time_col,
                        term = term,
                        nrand = nrand)
  
  # Compute confidence intervals if necessary
  if ( "rand" %in% names(res) ) {
    
    # Check that confidence level is valid
    if ( conf >= 1 | conf <= 0) {
      stop("Confidence level must be between 0 and 1")
    }
    
    # Get lower and upper CIs
    lower <- apply(res$rand, 2, quantile, (1-conf)/2)
    upper <- apply(res$rand, 2, quantile, conf + (1-conf)/2)
    
    # Create DF with results
    res_df <- data.frame(CV2 = res$CV2,
                        Pop.var = res$pop_var,
                        # Segrestin & Leps
                        Pop.sync = res$Pop_sync_direct[1]+res$Pop_sync_intra[1]+res$Pop_sync_indirect[1]+res$Pop_sync_no[1],
                        Pop.sync_l = sum(lower),
                        Pop.sync_u = sum(upper),
                        direct = res$Pop_sync_direct[1],
                        direct_l = lower[1],
                        direct_u = upper[1],
                        intra = res$Pop_sync_intra[1],
                        intra_l = lower[2],
                        intra_u = upper[2],
                        indirect = res$Pop_sync_indirect[1],
                        indirect_l = lower[3],
                        idirect_u = upper[4],
                        no = res$Pop_sync_no[1],
                        no_l = lower[4],
                        no_u = upper[4],
                        # Hammond et al
                        delta = res$Pop_sync_direct[3],
                        MP = res$Pop_sync_intra[3],
                        CCi = res$Pop_sync_indirect[3],
                        CCno = res$Pop_sync_no[3])
    # Add confidence level to colnames
    colnames(res_df)[grepl("_", colnames(res_df))] = paste0(colnames(res_df)[grepl("_", colnames(res_df))],
                                                                100*conf)
    
  } else {
    res_df <- data.frame(CV2 = res$CV2,
                        Pop.var = res$pop_var,
                        # Segrestin & Leps
                        Pop.sync = res$Pop_sync_direct[1]+res$Pop_sync_intra[1]+res$Pop_sync_indirect[1]+res$Pop_sync_no[1],
                        direct = res$Pop_sync_direct[1],
                        intra = res$Pop_sync_intra[1],
                        indirect = res$Pop_sync_indirect[1],
                        no = res$Pop_sync_no[1],
                        # Hammond et al
                        delta = res$Pop_sync_direct[3],
                        MP = res$Pop_sync_intra[3],
                        CCi = res$Pop_sync_indirect[3],
                        CCno = res$Pop_sync_no[3])
    }
  # Remove rownames
  rownames(res_df) <- NULL
  
  # Return results based on desired method
  # All
  if( sum(method %in% c("segrestin", "hammond")) == 2 ) {
    return(res_df)
  }
  # Segrestin & Leps 2022
  if( method == "segrestin") {
    return(res_df[, -c("delta", "MP", "CCi", "CCno")])
  }
  # Hammond et al. 2020
  if( method == "hammond" ) {
    return(res_df[, c("CV2", "Pop.var", "delta", "MP", "CCi", "CCno")])
  }
  
}
