# Special functions related to meta-analysis

## Report lines of meta-analyses in ways easy to copy paste
ma_report <- function(m) {
  sapply(
    1:nrow(m$b),
    function(x){
      sprintf("(%.0f) %.2f (95%% CI: %.2f, %.2f)", x, m$b[x], m$ci.lb[x], m$ci.ub[x])
    })
}

## Get weights from MLM and random effects meta-analyses
# This is how one can extract weights for a multilevel model according to GPT4
# https://chat.openai.com/share/a9644929-8f55-426e-b2ad-5137bc9be63d
ma_extract_weights <- function(m) {
  # Extract the fitted var-cov matrix of the random effects
  V <- m$V
  # Calculate the marginal variances
  margvar <- diag(V)
  # Calculate the weights
  weights <- 1 / margvar
  # Create a data frame with the study labels and weights
  study_weights <- data.frame(
    authors = m$slab,
    weight  = weights
  ) %>% mutate(weight_standardised = weight / sum(weight))

  return(study_weights)
}

## Function to do model comparison easily
# Model comparison is based on:
#  - "Data Analysis: A model Comparison Approach to regression, ANOVA,
# and Beyond" by Judd et al. (2017)
#  - https://github.com/StatQuest/logistic_regression_demo/blob/master/logistic_regression_demo.R
ma_model_comparison <- function (m1, m2){
  
  # Get which model is the reduced (C) or full (A)
  if(m1$parms <= m2$parms) {
    mC <- m1
    mA <- m2
  } else {
    mC <- m2
    mA <- m1
  }
  
  # Calculate the logLik difference
  df <- mA$parms - mC$parms
  logLik_C <- -2*logLik(mC)[1]
  logLik_A <- -2*logLik(mA)[1]
  logLik_diff <- abs(logLik_C - logLik_A) # abs to deal with rare case of positive logLik
  PRD <- abs(logLik_diff / logLik_C)
  
  # Check whether the logLik difference is statistically significant
  chi_crit <- qchisq(1-.05, df = df)
  p_value <- pchisq(logLik_diff, df=df, lower.tail=F)
  
  # Present the models and their values
  models_table <- data.frame(
    n_parameters_C = mC$parms,
    n_parameters_A = mA$parms,
    logLik_C,
    logLik_A,
    AIC_C = AIC(mC),
    AIC_A = AIC(mA)
  ) %>% pivot_longer(
    cols = everything(),
    names_to = c(".value", "model"),
    names_pattern = "(.*)_(.*)"
  ); models_table <- as.data.frame(models_table)
  
  # Present the logLik comparison
  output_table <- data.frame(
    PRD=PRD,
    logLik_diff=logLik_diff,
    df=df,
    chi_crit=chi_crit,
    p=p_value
  )

  # Add synthesis about the model comparison
  synthesis <- data.frame(
    bigger_model_better_logLik = p_value < 0.05,
    bigger_model_better_AIC = AIC(mA) < AIC(mC)
  ) %>% mutate(
    bigger_model_better_both = bigger_model_better_logLik & bigger_model_better_AIC
  )
  
  # Return a list with information calculated here
  # and information about the model specifications
  return(
    list(
      modelC = mC[["call"]],
      modelA = mA[["call"]],
      models_table,
      output_table,
      synthesis
    )
  )
}

# Function to get the SE modified by the Knapp-Hartung t confidence intervals
# So that the monte carlo simulations approach the modified confidence interval
# appropriately, which is usually larger than using the SE naively.
ma_get_adjusted_se <- function(m) {
  # Is it a mlm format model (i.e., rma.mv)
  is_mlm <- !is.null(m$random)
  
  # Get the total number of parameter
  params <- m$p
  
  # Get the total between-effect-sizes variance
  if(is_mlm) (tau2 <- sum(m$sigma2)) else (tau2 <- m$tau2)
  
  # Prepare the data
  adjusted_se <- NULL
  # Loop through and populate the data
  for(i in 1:params){

    # Adjust ddf based on whether it is mlm or not
    if(is_mlm) (ddf <- m$ddf[[i]]) else (ddf <- m$ddf[1])
    
    # Get the prediction interval
    t_crit = qt(0.05/2, df = ddf, lower.tail = FALSE)
    pi_lb = m$b[i] - t_crit * sqrt(tau2 + m$se[i]^2)
    pi_ub = m$b[i] + t_crit * sqrt(tau2 + m$se[i]^2)
    
    new_adjusted_se <- data.frame(
      param = i,
      se    = m$se[i],
      ci_lb = m$ci.lb[i],
      ci_up = m$ci.ub[i],
      # This is how we get a more representation se for the MC
      adjusted_se = (m$ci.ub[i]-m$ci.lb[i])/3.92,
      # Getting some prediction interval corrected se
      pi_lb,
      pi_ub,
      adjusted_se_pi = (pi_ub-pi_lb)/3.92
    )
    
    # Adding the lines
    adjusted_se <- rbind(adjusted_se, new_adjusted_se)
  }
  
  return(adjusted_se)
}