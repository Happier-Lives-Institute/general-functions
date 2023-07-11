# Functions to helps us read, write, and report answers.

# Round number to string, fills in missing digits (decimals) with 0s
round_c <- function(x, digits = 2) {
  roundedElement <- round(x, digits)
  roundedElement.char <- as.character(roundedElement)
  
  # If digits are more than 0
  if(digits > 0) {
    # If no decimals
    if(roundedElement %% 1 == 0) {
      roundedElement.char <- paste0(roundedElement.char, ".")
      for (i in 1:digits) {
        roundedElement.char <- paste0(roundedElement.char, "0")
      }
    } else {
      ndecimals <- nchar(strsplit(roundedElement.char, ".", fixed = T)[[1]][2])
      diffToDigits <- digits - ndecimals
      
      if(diffToDigits > 0) {
        for (i in 1:diffToDigits) {
          roundedElement.char <- paste0(roundedElement.char, "0")
        }
      }
    }
  }
  
  return(roundedElement.char)
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
    
    newCI <- paste0(
      "(95% CI: ", 
      round_c(quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[1]], 2), ", ",
      round_c(quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[2]], 2), ")"
    )
    
    df$ci[which(df$variable == variable)] <- newCI
  }
  
  # Copy into one string, for ease fo copy pasting
  df <- df %>% rowwise() %>% 
    mutate(combined = paste0(round_c(pe, 2), " ", ci)) %>% 
    ungroup()
  
  return(df)
}

# Report lines of meta-analyses
ma_report <- function(m) {
  sapply(
    1:nrow(m$b),
    function(x){
      sprintf("(%.0f) %.2f (95%% CI: %.2f, %.2f)", x, m$b[x], m$ci.lb[x], m$ci.ub[x])
    })
}
