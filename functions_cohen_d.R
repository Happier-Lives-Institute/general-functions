# Functions related to Cohen's d calculations

# Get pooled SD
get_SDp <- function(n1, sd1, n2, sd2){
  return(
    abs(
      sqrt(
        (((n1-1)*(sd1^2))+((n2-1)*(sd2^2))) / (n1 + n2 - 2)
      )
    )
  )
}

# Get Cohen's D
get_d <- function(m1, m2, pooledSD){
  return(
    (m1 - m2) / pooledSD
  )
}

# Get standard error of Cohen's D (see Harrer et al., 2021 for example)
get_d_se <- function(d, n1, n2){
  return(
    sqrt(((n1 + n2)/ (n1*n2)) + ((d^2) / (2*(n1 + n2))))
  )
}

# See Julio et al. 2003
# https://www.um.es/metaanalysis/pdf/7078.pdf
get_d_se_dhh <- function(yes1, no1, yes2, no2){
  return(
    sqrt((3/sqrt(pi))*((1/yes1)+(1/no1)+(1/yes2)+(1/no2)))
  )
}
