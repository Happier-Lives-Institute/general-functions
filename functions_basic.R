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
  if (x < 0) {
    return(0)
  } else {
    return(x)
  }
}