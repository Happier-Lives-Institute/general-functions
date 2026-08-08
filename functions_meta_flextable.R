#~############################################################################~#
# split_soft_returns ----
#~############################################################################~#

# Fix some issues with soft returns
split_soft_returns <- function(ft, data, j) {
  col <- data[[j]]
  rows_with_break <- which(grepl("\n", col, fixed = TRUE))

  for (i in rows_with_break) {
    parts <- strsplit(col[i], "\n", fixed = TRUE)[[1]]
    chunks <- list(parts[1])
    for (p in parts[-1]) {
      chunks <- c(chunks, list(as_chunk("\n")), list(p))
    }
    ft <- mk_par(ft, i = i, j = j, value = do.call(as_paragraph, chunks))
  }

  ft
}

#~############################################################################~#
# Function for Bayesian models ----
#~############################################################################~#

stan_rows <- function(
    fit, m_name,
    bayes_coefficients,
    bayes_n_name,
    bayes_simple_model,
    sigma_pattern = "sigma",
    hdi_prob = 0.95,
    simplified = FALSE
) {

  draws <- fit %>% rstan::extract() %>% as.data.frame()

  ## Coefficients
  fmt <- function(p) {
    x  <- draws[, p]
    ci <- HDInterval::hdi(x, credMass = hdi_prob)
    sprintf("%.2f\n(%.2f, %.2f)", mean(x), ci[["lower"]], ci[["upper"]])
  }
  coef_rows <- do.call(rbind, lapply(names(bayes_coefficients), function(p)
    data.frame(model_name = m_name, variable = bayes_coefficients[[p]], value = fmt(p))))

  ## total tau2
  total_tau2 <- function(d) {
    cols <- grep(sigma_pattern, colnames(d), value = TRUE)
    mean(rowSums(d[, cols, drop = FALSE]^2))
  }
  tau2_full <- total_tau2(draws)

  ## R2 vs baseline
  # bayes_simple_model is optional when simplified = TRUE
  if (!simplified) {
    draws_simple <- bayes_simple_model %>%
      rstan::extract() %>%
      as.data.frame()
    tau2_simple <- total_tau2(draws_simple)
    R2 <- positive_or_zero(
      (1 - (tau2_full/tau2_simple))
    )
  }

  ## Combine the general info
  if (simplified) {
    general_info <- data.frame(model_name = m_name,
                               variable = "Effect sizes", value = round_c(fit@par_dims[[bayes_n_name]],0))
  } else {
    general_info <- data.frame(
      model_name = m_name,
      variable = c("Effect sizes", "Tau²", "R² (Tau² reduction)", "AIC"),
      value    = c(round_c(fit@par_dims[[bayes_n_name]],0), round_c(tau2_full), round_per(R2), NA_character_)
    )
  }

  list(coefficient_rows = coef_rows,
       general_info = general_info,
       coef_vars = coef_rows$variable)
}

#~############################################################################~#
# meta_flextable ----
#~############################################################################~#
# Our own function to present tables of information as we need.

meta_flextable <- function(
    model_list,
    additional_variables = NULL,
    font_size = 10,
    add_spaces = 0,
    my_study_label = "study_label",
    simplified = FALSE,
    knots_list = NULL,
    label_fun = identity,
    bayes_n_name_list = NULL,
    bayes_coefficients_list = NULL,
    bayes_simple_model_list = NULL,
    sigma_pattern = "sigma"
) {

  #~=======================================================~=
  ## Go through each model ----
  #~=======================================================~=

  # preparation
  all_models <- list()
  coef_list <- c()
  other_variables <- c()

  for (i in seq_along(model_list)) {

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Get the model ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Get the model
    m <- model_list[[i]]
    # Get the model name
    m_name <- names(model_list)[i]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Catch if this is a Bayesian model ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (inherits(m, "stanfit")) {
      sr <- stan_rows(
        m, m_name,
        bayes_n_name = bayes_n_name_list[[m_name]],
        bayes_coefficients = bayes_coefficients_list[[m_name]],
        bayes_simple_model = bayes_simple_model_list[[m_name]],
        sigma_pattern = sigma_pattern,
        simplified = simplified
      )
      m_rows          <- sr$coefficient_rows
      new_coef_list   <- sr$coef_vars
      coef_list       <- c(coef_list, setdiff(new_coef_list, coef_list))

      if(!is.null(additional_variables)) {
        for (add_i in seq_along(additional_variables)) {
          # Combine all the value
          m_row <- data.frame(
            model_name = m_name,
            variable = names(additional_variables)[add_i],
            value = (additional_variables[[add_i]][i])

          )

          # Merge
          m_rows <- rbind(m_rows, m_row)
        }
      }

      m_rows <- rbind(m_rows, sr$general_info)
      all_models <- append(all_models, list(m_rows))

      # Prepare for the table (union, so rows only some models have are kept)
      other_variables <- unique(c(
        other_variables, setdiff(m_rows$variable, new_coef_list)
      ))

      next
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Get each coefficient ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Get row per parameter
    m_rows <- NULL

    m_knots    <- if (!is.null(knots_list)) knots_list[[m_name]] else NULL
    coef_names <- coefficient_names(m)

    for (ii in seq_len(nrow(m$b))) {

      variable_name <- coef_names[ii]

      # Combine all the value
      m_row <- data.frame(
        model_name = m_name,
        variable = (variable_name),
        value = (sprintf(
          "%.2f%s\n(%.2f, %.2f)",
          m$b[ii], get_sig_star(m$pval[ii]), m$ci.lb[ii], m$ci.ub[ii]
        ))
      )

      # Merge
      m_rows <- rbind(m_rows, m_row)
    }

    # General coefficient list information to make the table
    new_coef_list <- m_rows$variable
    coef_list <- c(coef_list, setdiff(new_coef_list, coef_list))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Add additional variables ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if(!is.null(additional_variables)) {
      for (add_i in seq_along(additional_variables)) {
        # Combine all the value
        m_row <- data.frame(
          model_name = m_name,
          variable = names(additional_variables)[add_i],
          value = (additional_variables[[add_i]][i])

        )

        # Merge
        m_rows <- rbind(m_rows, m_row)
      }
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Add general information ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    m_report <- ma_report_analysis(m)
    data_single <- m$data %>% filter(!duplicated(get(my_study_label)))

    # Combine all the value
    if(simplified) {
      m_row <- data.frame(
        model_name = m_name,
        variable = c(
          "Effect sizes"
        ),
        value = c(
          m$k.eff
        )
      )
    } else {
      m_row <- data.frame(
        model_name = m_name,
        variable = c(
          "Effect sizes",
          "Tau²",
          "R² (Tau² reduction)",
          "AIC"
        ),
        value = c(
          m$k.eff,
          round_c(m_report$model$tau2),
          round_per(m_report$model$R2),
          round_c(AIC(m), 0)
        )
      )
    }

    # Merge
    m_rows <- rbind(m_rows, m_row)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Add spline knot locations (if provided) ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!is.null(m_knots)) {
      m_spline_var <- spline_var_of(m)
      m_rows <- rbind(m_rows, data.frame(
        model_name = m_name,
        variable   = if (nzchar(m_spline_var)) {
          paste0("Spline knots (", m_spline_var, ")")
        } else "Spline knots",
        value      = paste(round_c(m_knots, 2), collapse = ", ")
      ))
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### Save all models ----
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    all_models <- append(all_models, list(m_rows))

    # Prepare for the table (union, so rows only some models have are kept)
    other_variables <- unique(c(
      other_variables, setdiff(m_rows$variable, new_coef_list)
    ))
  }

  #~=======================================================~=
  ## Merge all the models ----
  #~=======================================================~=

  table_models <- data.frame(Variable = c(coef_list, other_variables))
  for (m in all_models) {

    m_merge <- m %>% select(-model_name)
    colnames(m_merge) <- c("Variable", m$model_name[1])

    table_models <- left_join(
      table_models,
      m_merge,
      by = "Variable"
    )
  }

  #~=======================================================~=
  ## Edit parameter names ----
  #~=======================================================~=

  table_models <- table_models %>%
    mutate(Variable = label_fun(Variable))

  #~=======================================================~=
  ## Make into flextable ----
  #~=======================================================~=

  table <- flextable(table_models) %>%
    set_table_properties(width = 1, layout = "autofit")

  # Align centre some cells
  names_to_align <- colnames(table_models %>% select(-Variable))
  table <- table %>% align(j = names_to_align, align = "center")

  # centre headers for model columns
  header_cols <- setdiff(colnames(table_models), "Variable")
  table <- align(table, j = header_cols, align = "center", part = "header")
  table <- align(table, j = "Variable", align = "left", part = "header")

  # Add separation line
  table <- hline(
    table,
    i = length(coef_list),
    j = 1:ncol(table_models),
    border = officer::fp_border(color = "black", width = 1.5),
    part = "body"
  )

  # Add separation line if there are additional variables
  if(!is.null(additional_variables)) {
    table <- hline(
      table,
      i = length(coef_list) + length(additional_variables),
      j = 1:ncol(table_models),
      border = officer::fp_border(color = "black", width = 1.5),
      part = "body"
    )
  }

  # Fix soft returns
  table <- split_soft_returns(table, table_models, "Variable")

  # Change font size
  table <- fontsize(table, size = font_size, part = "all")

  return(table)
}
