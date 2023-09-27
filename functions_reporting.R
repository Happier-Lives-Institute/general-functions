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
  df$ci_num <- NA
  df$ci_lower <- NA
  df$ci_upper <- NA
  for (variable in df$variable) {
    lower  <- quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[1]]
    upper  <- quantile(as.data.frame(sim)[,variable], 
                       probs = c(CI.low, CI.upp))[[2]]
    new_ci <- sprintf("(%i%s CI: % .2f, % .2f)", (ci * 100), "%", lower, upper)
    df$ci_lower[which(df$variable == variable)] <- lower
    df$ci_upper[which(df$variable == variable)] <- upper
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

# A function for lineseparations in rmd kableExtra tables
generate_linesep <- function(table, freq, insert = "\\addlinespace") {
  # Initialize an empty vector
  linesep_vec <- rep("", nrow(table))
  
  # Insert "\\addlinespace" (or other specified) at the specified frequency
  linesep_vec[seq(freq, nrow(table), by = freq)] <- insert
  
  return(linesep_vec)
}

# Function to get the total effect with all the flexible methods we've used up to now
get_total_effect <- function(
    initial_centre, # mean
    initial_spread, # se
    initial_ll = -Inf, # lower limit
    initial_ul = +Inf, # upper limit
    trajectory_centre, # mean
    trajectory_spread, # se
    trajectory_ll = -Inf, # lower limit
    trajectory_ul = +Inf, # upper limit
    duration, # either individual_0, overall_0, or a value
    WELLBY_conversion_rate = 1, # rate at which SDs are converted to WELLBYs
    simulations, # number of simulations
    illustration_fraction = 1, # fraction of simulations to take for graph
    seed = NULL # seed for simulations
) {
  
  # Run a potential seed
  if(!is.null(seed)){set.seed(seed)}
  
  ## Prepare the first elements
  total_effect_pe <- tibble(
    initial = initial_centre,
    trajectory = trajectory_centre
  )
  
  total_effect_sim <- tibble(
    initial = msm::rtnorm(
      n = simulations,
      mean = initial_centre,
      sd = initial_spread,
      lower = initial_ll,
      upper = initial_ul
    ),
    trajectory = msm::rtnorm(
      n = simulations,
      mean = trajectory_centre,
      sd = trajectory_spread,
      lower = trajectory_ll,
      upper = trajectory_ul
    )
  )
  
  ## Set the duration
  # If each simulation gets its own individual duration
  if(duration == "individual_0") {
    total_effect_pe <- total_effect_pe %>% mutate(
      duration = abs(initial/trajectory)
    )
    
    total_effect_sim <- total_effect_sim %>% mutate(
      duration = abs(initial/trajectory)
    )
  } else {
    # If the duration is set to a certain value
    # If it is set to when the point estimates reach 0
    if(duration == "overall_0") {
      duration_value = abs(initial_centre/trajectory_centre)
    } else {
      # otherwise set to value of the user
      duration_value = duration
    }
    
    total_effect_pe <- total_effect_pe %>% mutate(
      duration = duration_value
    )
    total_effect_sim <- total_effect_sim %>% mutate(
      duration = duration_value
    )
  }
  
  ## Calculate the total effect
  total_effect_pe <- total_effect_pe %>% rowwise() %>% mutate(
    total_effect = pracma::integral(function(t) {initial + trajectory * t}, 
                                    0, duration),
    total_effect_wellbys = total_effect*WELLBY_conversion_rate
  ) %>% ungroup()
  
  total_effect_sim <- total_effect_sim %>% rowwise() %>% mutate(
    total_effect = pracma::integral(function(t) {initial + trajectory * t}, 
                                    0, duration),
    total_effect_wellbys = total_effect*WELLBY_conversion_rate
  ) %>% ungroup()
  
  ## Combine the two
  total_effect_combined <- combine_PE_SIM(total_effect_pe, total_effect_sim)
  
  ## Make illustration graph
  
  # Make graph
  total_effect_graph <- total_effect_sim %>% 
    sample_frac(illustration_fraction) %>% 
    mutate(
      id = 1:(simulations*illustration_fraction),
      end_effect = initial+trajectory*duration
    ) %>% ggplot() +
    geom_hline(yintercept = 0, linetype=2)+
    geom_segment(aes(
      x=0, xend=duration,
      y=initial, yend=end_effect
    ), alpha=0.25, color="#3167b4") +
    geom_segment(aes(
      x=0, xend=total_effect_pe$duration,
      y=total_effect_pe$initial, 
      yend=total_effect_pe$initial+
        total_effect_pe$trajectory*total_effect_pe$duration
    ), alpha=1, linewidth=1, color="orange") +
    theme_cowplot() +
    ylab("effect") + xlab("time")
  
  return(
    list(
      total_effect_pe       = total_effect_pe,
      total_effect_sim      = total_effect_sim,
      total_effect_combined = total_effect_combined,
      total_effect_graph    = total_effect_graph
    )
  )
}
