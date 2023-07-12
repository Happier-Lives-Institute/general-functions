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
