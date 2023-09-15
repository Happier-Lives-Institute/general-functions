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

## Report a fuller analysis of the RE or MLM model
# Reports information about the moderation if relevant
# Reports information about heterogeneity
# Can flexibly do rma.mv and rma

ma_report_analysis <- function(m) {
  
  #---
  # general variables
  #---
  
  # Is it a mlm format model (i.e., rma.mv)
  is_mlm <- !is.null(m$random)
  
  # Number of effect sizes
  k <- m$k.all
  
  # The weights based on sampling variance
  wi <- 1/m$vi
  
  #---
  # within-effect-sizes variance (nu)
  #---
  # The first step for most heterogeneity calculations is to get the typical
  # within-effect-sizes variance.
  
  # This calculation is based on exploring the calculation in metafor::rma
  # It looks complicated, because it takes into account possible moderators.
  
  # Create a diagonal matrix of the weights (which will weight the contributions
  # of each effect size based on this precision)
  W <- diag(wi, nrow = k, ncol = k)
  
  # Create a design matrix of the model (intercept plus any moderators)
  # If there are no moderators you need to tell it to be an intercept only model
  # First, we need to prepare the data, filtering NAs for any of the moderators
  if(is.null(m$formula.mods)) {
    model_formula <- as.formula("~1")
    m_data <- m$data
  } else {
    model_formula <- m$formula.mods
    m_data <- m$data %>%
    filter(!if_any(all_of(all.vars(model_formula)), is.na))
  }

  # Then we can create the model matrix
  X <- model.matrix(model_formula, data = m_data)
  
  # Copying the rma function to inverse a product of matrices
  .invcalc <- function (X, W, k) {
    sWX <- sqrt(W) %*% X
    res.qrs <- qr.solve(sWX, diag(k))
    return(tcrossprod(res.qrs))
  }
  
  # Get the inverse matrix of the variance weihts by the moderation
  inverse_XWX <- .invcalc(X = X, W = W, k = k)
  
  # Get a projection matrix of the previously calculated matrices
  # (this often used to calculate residuals or other quantities that measure 
  # the fit of a model)
  P <- W - W %*% X %*% inverse_XWX %*% crossprod(X, W)
  
  # Get the number of predictors in the design matrix of this model
  predictors_n <- NCOL(X)
  
  # Copying the rma function to return the trace of a matrix 
  # (sum of diagonal elements)
  .tr <- function (X) {return(sum(diag(X)))}
  
  # Calculate the within-effect-sizes variance (nu)
  nu <- (k - predictors_n)/.tr(P)
  
  #---
  # levels and tau
  #---
  # Prepare the levels structure that will be reported. 
  # Get the tau2s for each level which is important for further calculations.
  if(is_mlm) {
    # Total number of levels set in the 
    max_levels <- length(m$s.names)+1
    
    # Prepare the levels
    model_levels <- data.frame(
      level_number = seq(max_levels,1, -1),
      level_formula = NA,
      level_n = NA,
      variance_symbol = NA, 
      variance = NA,
      nu = NA,
      tau2 = NA
    )
    
    # Add the descriptions of the levels
    # Because the levels are written from the top to the bottom we need to go
    # through these in the opposite direction.
    for (level_i in seq(max_levels,1, -1)) {
      if(level_i == max_levels) {
        # First level is within-effect-size variance
        model_levels$level_formula[level_i] <- "within-effect-size variance"
        model_levels$level_n[level_i] <- k
        model_levels$variance_symbol[level_i] <- "ε"
        model_levels$variance[level_i] <- nu
        model_levels$nu[level_i] <- nu
        model_levels$tau2[level_i] <- NA
      } else {
        # Other levels are dependent on the formulation given in the model
        # (-1 because the formulations don't have the within-effect-size variance)
        model_levels$level_formula[level_i] <- m$s.names[level_i]
        model_levels$level_n[level_i] <- m$s.nlevels[level_i]
        model_levels$variance_symbol[level_i] <- paste0(
          "ζ", model_levels$level_number[level_i]-1
        )
        model_levels$variance[level_i] <- m$sigma2[level_i]
        model_levels$nu[level_i] <- NA
        model_levels$tau2[level_i] <- m$sigma2[level_i]
      }
    }
    
  } else {
    max_levels <- 2
    # Only 2 levels if it is a random effects model
    model_levels <- data.frame(
      level_number = c(2, 1),
      level_formula = c("between-effect-size-variance", "within-effect-size variance"),
      level_n = k,
      variance_symbol = c("ζ","ε"),
      variance = c(m$tau2, nu),
      nu = c(NA, nu),
      tau2 = c(m$tau2, NA)
    )
  }
  
  #---
  # I2
  #---
  # We need the sum of all the variances
  sum_variance <- sum(model_levels$variance, na.rm = T)
  # Get the I2 for each level
  model_levels$I2 <- model_levels$tau2 / sum_variance
  model_levels$`% of total variance` <- round_per(model_levels$I2)
  model_levels$`% of total variance`[max_levels] <- round_per((1 - sum(model_levels$I2, na.rm = T)))
  
  #---
  # ICC
  #---
  # % of the between-effect-size variance that is due to each level
  model_levels$icc <- model_levels$tau2 / sum(model_levels$tau2, na.rm = T)
  
  #---
  # H2
  #---
  # H2 is a ratio of all variance to the within-effect-size variance
  model_levels$H2 <- sum_variance / nu
  
  #---
  # Cochrane's Q
  #---
  
  # Present all the information for Cochrane's Q test
  Q <- list(
    data = data.frame(Q = m$QE, df = k-predictors_n, p = m$QEp),
    summary = paste0("Q(df = ", k-predictors_n, ") = ", round_c(m$QE), ", ", format_p_value(m$QEp))
  )
  
  #---
  # Prediction intervals
  #---
  # Get prediction intervals. It is unclear if this is exactly how it should
  # be done for analyses with moderators.
  
  pi <- ma_get_adjusted_se(m)
  
  #---
  # Moderators information
  #---
  
  ## Calculate R2
  
  # Prepare the null model (model c, the version without moderators)
  # If the model is not mlm, thereby an RE rma, we need to set the random level
  # manually.
  if(is_mlm) {
    model_null_data <- m_data
    model_null_random <- as.formula(m$random[[1]])
  } else {
    model_null_data <- m_data
    model_null_data$es_id <- 1:nrow(model_null_data)
    model_null_random <- as.formula("~1 | es_id")
  }
  model_null <- rma.mv(
    yi = m$yi, V = m$vi, random = model_null_random, data = model_null_data,
    method = m$method, test = m$test
  )
  
  # Calculate R2 for the different models
  model_levels$R2 <- NA
  for (level_i in seq(max_levels,1, -1)) {
    model_levels$R2[level_i] <- positive_or_zero(
      (1 - (model_levels$tau2[level_i]/model_null$sigma2[level_i]))
    )
  }
  
  # Calculate the general R2
  R2 <- positive_or_zero(
    (1 - (sum(model_levels$tau2, na.rm = T)/sum(model_null$sigma2)))
  )
  
  # Calculate overall how much the tau2 reduces from adding moderators
  tau2_reduction_from_moderators <- sum(model_null$sigma2) - sum(model_levels$tau2, na.rm = T)
  
  # present all the information about the moderators
  model <- list(
    moderators_formula = model_formula,
    parameters_n = predictors_n,
    R2 = R2,
    tau2 = sum(model_levels$tau2, na.rm = T),
    tau2_reduction_from_moderators = tau2_reduction_from_moderators
  )
  
  #---
  # Prepare information for return
  #---
  results <- list(
    model = model,
    Q = Q,
    model_levels = model_levels,
    pi = pi
  )
  
  return(results)
}

# Function to produce a ggplot funnel plot
my_funnel_rma.mv <- function(
    m, # the model we are graphing for
    add_contour = T, # Whether there is a contour plot
    se_variable = "function_se", # the variable to use if we want to set a different se, 
    # such as corrected SEs
    point_fill = "" # Variable by which we group the effect sizes in color
) {
  
  # Extract data and estimates from rma.mv model
  dat <- m$data
  dat$yi <- as.numeric(m$yi)
  dat$function_se <- sqrt(m$vi)
  dat$function_se <- dat[, se_variable][[1]]
  
  # Prepare general elements of the graph
  pooled_effect <- m$b[1]
  max_se <- max(dat$function_se)
  
  # Calculate expected (or pseudo) confidence interval for the maximum SE
  expected_ci <- calculate_ci_se(pooled_effect, max_se)
  
  # Generate the funnel plot
  my_funnel_plot <- dat %>% ggplot(aes(x = yi, y = function_se)) +
    # Reverse the standard errors
    scale_y_reverse() +
    
    # Prepare the theme
    theme_cowplot() +
    labs(x = "Effect size",
         y = "Standard Error")
  
  # Conditionally add the contour plot
  if (add_contour) {
    contour_data <- data.frame(
      x = c(
        c(0, -1.64 * max_se, +1.64 * max_se),
        c(0, -1.64 * max_se, -1.96 * max_se), 
        c(0, 1.64 * max_se, 1.96 * max_se), 
        c(0, -2.58 * max_se, -1.96 * max_se),
        c(0, 2.58 * max_se, 1.96 * max_se),
        c(0, -2.58 * max_se, -Inf),
        c(0, 2.58 * max_se, Inf)
      ),
      y = rep(c(0, Inf, Inf), 7),
      fill = c(
        rep("n.s.", 3),
        rep("p < .10", 3),
        rep("p < .10", 3),
        rep("p < .05", 3),
        rep("p < .05", 3),
        rep("p < .001", 3),
        rep("p < .001", 3)
      )
    )
    
    my_funnel_plot <- my_funnel_plot +
      geom_polygon(data = contour_data, aes(x = x, y = y, fill = fill)) +
      scale_fill_manual(
        values = c("n.s." = "gray75", "p < .10" = "gray85", 
                   "p < .05" = "gray95", "p < .001" = "white"),
        breaks = c("p < .001", "p < .05", "p < .10", "n.s."),
        name = "Significance Level",
        drop = F
      )
  }
  
  # Continue with other elements of the funnel plot
  my_funnel_plot <- my_funnel_plot +
    # Create the line for the predicted effect size
    geom_segment(aes(
      x = pooled_effect, xend = pooled_effect, y = Inf, yend = 0
    ), linetype = 3) +
    
    # Create the funnel
    geom_segment(aes(
      x = expected_ci[1], xend = pooled_effect, y = Inf, yend = 0
    ), linetype=3) +
    geom_segment(aes(
      x = pooled_effect, xend = expected_ci[2], y = 0, yend = Inf
    ), linetype=3)
  
  # Dependent on whether a point fill was given, add the effect sizes
  if (point_fill != "") {
    my_funnel_plot <- my_funnel_plot +
      geom_point(alpha = 0.5, 
                 aes(fill = get(point_fill), color = get(point_fill))) +
      # need to modify the legend accordingly
      guides(
        fill = guide_legend(override.aes = list(alpha = 1, shape = NA)),
        color = guide_legend(title = point_fill)
      ) +
      theme(legend.position = "bottom", legend.box = "vertical")
  } else {
    my_funnel_plot <- my_funnel_plot +
      geom_point(alpha = 0.5, shape = 21, color = "black", fill = "blue")
  }
  
  # return the funnel
  return(my_funnel_plot)
  
}