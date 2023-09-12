# Functions to helps us read, write, and report answers.

# Round number to string, fills in missing digits (decimals) with 0s
round_c <- function(x, digits = 2) {
  formatC(x, format="f", digits=digits)
}

# Round number to string with % in front
round_per <- function(x) {
  percent_string <- sprintf("%.2f%%", x * 100)
  return(percent_string)
}

# Function to easily write a confidence interval string
present_with_CI <- function(estimate, lower, upper, per = "95%") {
  if(is.na(estimate)) {
    return(NA)
  } else if (is.na(lower) | is.na(upper)) {
    return(sprintf("% .2f", estimate))
  } else {
    return(sprintf("% .2f (%s CI: % .2f, % .2f)", estimate, per, lower, upper))
  }
}

# Function to report a distribution #
MMSDCI.vec <- function(x, ci = .95){
  
  # Get the upper and lower quantiles
  CI.low <- (1 - ci)/2
  CI.upp <- ci + CI.low
  
  # Return the summary
  return(
    paste0(
      "M = ", round_c(mean(x, na.rm=T), 2),
      ", SD = ", round_c(sd(x, na.rm=T), 2),
      ", Median = ", round_c(median(x, na.rm=T), 2),
      ", (95% CI: ", round_c(quantile(x, probs = c(CI.low, CI.upp))[[1]], 2), ", ",
      round_c(quantile(x, probs = c(CI.low, CI.upp))[[2]], 2), ")"
    )
  )
}

MMSDCI <- function(x, ci = .95){
  
  # If it is vector, just report one line
  if (is.vector(x)){
    report <- MMSDCI.vec(x, ci)
  } else {
    # Filter out non-numeric columns
    x <- x %>% select_if(is.numeric)
    
    # Otherwise, report a line for every column
    report <- lapply(x, MMSDCI.vec)
  }
  
  return(report)
}

# Combine PE and SIM
combine_PE_SIM <- function(pe, sim, ci = .95){
  
  # Get the upper and lower quantiles
  CI.low <- (1 - ci)/2
  CI.upp <- ci + CI.low
  
  # Prepare the new dataframe where each row is a variable
  # And we will have a pe, ci, and combined, column
  df <- pe %>% pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "pe"
  )
  
  # Get the CI for each variable
  df$ci <- NA
  for (variable in df$variable) {
    lower <- quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[1]]
    upper <- quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[2]]
    new_ci <-   sprintf("(%I CI: % .2f, % .2f)", (ci * 100), lower, upper)
    df$ci[which(df$variable == variable)] <- new_ci
  }
  
  # Copy into one string, for ease fo copy pasting
  df <- df %>% rowwise() %>% 
    mutate(combined = sprintf("% .2f %s", pe, ci)) %>% 
    ungroup()
  
  return(df)
}

# A function to format p values like it is common to report them
format_p_value <- function(p_value, include_p = T) {
  p_prefix <- if (include_p) "p " else ""
  p_equals <- if (include_p) "p = " else ""
  
  if (p_value < 0.001) {
    return(paste0(p_prefix, "< .001"))
  } else {
    formatted_p <- formatC(p_value, format = "f", digits = 3)
    return(paste0(p_equals, sub("^0\\.", ".", formatted_p)))
  }
}
