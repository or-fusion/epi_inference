require(MASS)
require(stringr)

sample_county_negbin <- function(countyfile, window = 3, n_samples = 100, read_dir, write_dir){
  dat <- read.csv(paste0(read_dir, countyfile), stringsAsFactors = F)
  range <- c((window + 1) : nrow(dat))     # Instead using a symmetric window until the end, then using past data
  
  # If the county has no cases, keep them all at zero
  if (dat$Confirmed[nrow(dat)] == 0){
    samples_negbin <- as.data.frame(matrix(data=0, nrow=length(range), ncol=n_samples))
  } else {
    daily <- c(dat$Confirmed[1], dat$Confirmed[2:nrow(dat)] - dat$Confirmed[1:(nrow(dat) - 1)])
    params <- as.data.frame(matrix(data=NA, ncol=3))
    samples_negbin <- as.data.frame(matrix(data=NA, ncol=n_samples))
    r <- 1
    
    for (i in (window + 1):nrow(dat)){
      if (i > (nrow(dat) - window)){
        window_data <- daily[(length(daily) - (2*window + 1)) : length(daily)]
      } else {
        # Using a symmetric window (window size is number of days on either side of date of interest)
        window_data <- daily[(i - window):(i + window)]
      }
      if (all(window_data == 0)){
        # Need to force the negative binomial parameters to get a fit in some cases
        p2 <- 0
        p3 <- 1
      } else {
        if (min(window_data) < 2){
          low <- 0.1
        } else {
          low <- 1
        }
        fit <- TRUE
        fit <- tryCatch(fitdistr(window_data, 'Negative Binomial', lower = low), 
                        error = function(cond){
                          return(fitted = FALSE)
                        })
        if (length(fit) == 1){
          p2 <- mean(window_data)
          p3 <- sd(window_data)
        } else {
          p2 <- as.numeric(fitdistr(window_data, 'Negative Binomial', lower = low)$estimate[2])
          p3 <- as.numeric(fitdistr(window_data, 'Negative Binomial', lower = low)$estimate[1])
        }
      }
      params[r,] <- c(as.character(dat$Date[i]), p2, p3)
      samples_negbin[r,] <- rnegbin(n_samples, p2, p3)
      r <- r + 1
    }
  }
  for (j in 1:n_samples){
    new_df <- dat[range,]
    new_df$Confirmed <- cumsum(samples_negbin[,j])
    # Insert one day of zeroes before so that Bill's code works for all counties
    new_df <- rbind(dat[2,], new_df)
    new_df$Confirmed[1] <- 0
    new_df$FIPS <- as.character(new_df$FIPS)
    new_df$FIPS <- str_pad(new_df$FIPS, 5, pad='0')
    
    # Where to write these files
    # write_path <- paste0('../covid-data/formatted_data/county_data_resample/', folder_name, '/sample', str_pad(as.character(j), 3, pad='0'))
    write_path <- paste0(write_dir, '/sample', str_pad(as.character(j), 3, pad='0'))
    if (!dir.exists(file.path(write_path))){
      dir.create(file.path(write_path))
    }
    write.csv(new_df, file=paste0(write_path, '/', countyfile))
  }
}

# main_dir <- '../covid-data/formatted_data/county_data/'
# # main_dir <- './formatted_data/county_data/'
# folder_name <- list.files(main_dir)[length(list.files(main_dir))]
# path <- paste0(main_dir, folder_name, '/')
# county_files <- list.files(path)
# 
# write_path <- paste0('../covid-data/formatted_data/county_data_resample/', folder_name)
# # write_path <- paste0('./formatted_data/county_data_resample/', folder_name)
# if (!dir.exists(file.path(write_path))){
#   dir.create(file.path(write_path))
# }
# 
# lapply(county_files, sample_county_negbin)
