n_com <- 2
n_sps <- 10
n_year <- 15
d <- round(matrix(data = abs(rnorm(n = n_com*n_sps*n_year)),
                  nrow = n_year*n_com, 
                  ncol = n_sps)*10)
dd = as.data.frame(cbind(rep(seq_len(n_com), each = n_year), 
                         rep(seq_len(n_year), times = n_com),
                         d))
names(dd) <- c("comm", "time", paste0("sp_", seq_len(n_sps)))

# Create matrix with missing data
missing_mat <- matrix(data = sample(c(1, NA),
                                    size = n_sps*n_year*n_com, 
                                    replace = T,
                                    prob = c(0.8, 0.2)),
                      nrow = n_year*n_com, 
                      ncol = n_sps)
dd[-c(1:2)] <- dd[-c(1:2)] * missing_mat

# Shuffle rows to test ordering
dd <- dd[sample(nrow(dd)),]

dd_long <- reshape(data = dd,
                   direction = "long",
                   varying = 3:(n_sps+2),
                   v.names = "abundance",
                   times = names(dd)[3:(n_sps+2)],
                   timevar = "species")

dd_wide = reshape(dd_long,
                  direction = "wide",
                  timevar = "species",
                  idvar = "id",
                  v.names = "abundance")
dd_wide = dd_wide[!names(dd_wide) %in% c("id")]
dd_wide[,5] <- 0

dd_long <- dd_long[, !names(dd_long) == "id"]

