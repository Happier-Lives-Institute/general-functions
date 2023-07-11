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

# Get number of unique instances
unique_n <- function(x) {
  return (length(unique(x)))
}

# Get standard error
se <- function(x, na.rm = F) {
  if(na.rm) {x <- na.omit(x)}
  return(sd(x)/sqrt(length(x)))
}