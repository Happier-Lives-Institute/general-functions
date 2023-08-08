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
