#---------
# Commands
#---------

# Not in
'%ni%' <- Negate('%in%')

#---------
# Special
#---------

# Method declare to help get rid of NaNs
# https://stackoverflow.com/questions/52490552/r-convert-nan-to-na
is_nan_data_frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

#---------
# Functions
#---------

## Get number of unique instances
unique_n <- function(x) {
  return (length(unique(x)))
}

## Get standard error
se <- function(x, na.rm = F) {
  if(na.rm) {x <- na.omit(x)}
  return(sd(x)/sqrt(length(x)))
}

## Function to get SE from mean and CI
get_se_from_ci <- function(lower, upper){
  return((upper - lower)/3.92)
}

## Function for non-negative values
# Returns the valule if positive, returns 0 if negative
positive_or_zero <- function(x) {
  if (is.na(x)) {
    return(NA)
  } else if (x < 0) {
    return(0)
  } else {
    return(x)
  }
}

# Calculate confidence interval based on se
calculate_ci_se <- function(mean, se, level = 0.95) {
  # Calculate the Z-score based on the confidence level
  z_score <- qnorm(1 - (1 - level) / 2)
  
  # Calculate the confidence interval
  lower_bound <- mean - z_score * se
  upper_bound <- mean + z_score * se
  
  return(c(lower_bound, upper_bound))
}


# Formula to get the t easily
get_independent_t <- function(n1, n2, m1, m2, sd1, sd2){
  SD_p <- get_SD_p(n1, sd1, n2, sd2)
  SE_diff <- sqrt(((SD_p^2)/n1)+((SD_p^2)/n2))
  t = (m1 - m2) / SE_diff
  return(t)
}