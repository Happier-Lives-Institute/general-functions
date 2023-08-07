# Functions related to Cohen's d calculations

## MSD
# Get pooled SD
get_SD_p <- function(n1, sd1, n2, sd2){
  return(
    abs(
      sqrt(
        (((n1-1)*(sd1^2))+((n2-1)*(sd2^2))) / (n1 + n2 - 2)
      )
    )
  )
}

# Get Cohen's D
get_d <- function(m1, m2, SD_pooled){
  return(
    (m1 - m2) / SD_pooled
  )
}

# Get standard error of Cohen's D (see Harrer et al., 2021 for example)
get_d_se <- function(d, n1, n2){
  return(
    sqrt(((n1 + n2) / (n1*n2)) + ((d^2) / (2*(n1 + n2))))
  )
}

## convert OR into d
# See Sánchez-Meca et al. 2003
# https://www.um.es/metaanalysis/pdf/7078.pdf

# HH
get_d_hh <- function(OR) {
  return(log(OR)*(sqrt(3)/pi))
}

get_d_hh_se <- function(yes1, no1, yes2, no2){
  return(
    sqrt((3/sqrt(pi))*((1/yes1)+(1/no1)+(1/yes2)+(1/no2)))
  )
}

# Cox
get_d_Cox <- function(OR) {
  return(log(OR)/1.65)
}

get_d_Cox_se <- function(yes1, no1, yes2, no2){
  return(
    sqrt(0.367*((1/yes1)+(1/no1)+(1/yes2)+(1/no2)))
  )
}

## Get Hedge's g
get_g <- function(d, n1, n2) {
 return(d * (1-(3/(4*(n1+n2)-9))))
}

## Get se corrected for SMDs so that it is independent from effect size
get_corrected_se_smd <- function(n1, n2) {
  return(sqrt(((n1 + n2)/(n1*n2))))
}