#~############################################################################~#
# Fit splines ----
#~############################################################################~#
get_spline_ma <- function(
    dat,
    spline_var,
    spline_to_try,
    fit_rmv,
    extra_mods = NULL,
    use_local = FALSE
) {

  model_vars <- c(spline_var, extra_mods)
  dat <- dat %>% filter(if_all(all_of(model_vars), ~ !is.na(.x)))

  make_mods <- function(k) {
    terms <- c(sprintf("rcs(%s, %d)", spline_var, k), extra_mods)
    as.formula(paste("~", paste(terms, collapse = " + ")))
  }

  fit_splines <- function(k, method) {
    tryCatch(
      fit_rmv(dat = dat, mods = make_mods(k), method = method),
      error = function(e) {
        if (!use_local) message("        skipped ", k, " knots: ", conditionMessage(e))
        NULL
      }
    )
  }

  # With "ML" the logLik comparison across different fixed effects is valid,
  # and the winner is refitted with REML below.
  m_selection <- NULL
  n_splines   <- NA_integer_

  for (splines in spline_to_try) {
    new_m_selection <- fit_splines(splines, method = "ML")
    if (is.null(new_m_selection)) next

    if (is.null(m_selection)) {
      keep <- TRUE
    } else {
      my_model_comparison <- ma_model_comparison(m_selection, new_m_selection)
      keep <- my_model_comparison[[5]]$bigger_model_better_both
    }

    if (keep) {
      m_selection <- new_m_selection
      n_splines   <- splines
    }
  }

  if (is.null(m_selection)) {
    stop("No spline model converged for ", spline_var,
         " with knots: ", paste(spline_to_try, collapse = ", "))
  }

  m_splines <- fit_splines(n_splines, method = "REML")

  knots <- attr(rcs(dat[[spline_var]], n_splines), "parms")

  return(list(
    m_splines  = m_splines,
    n_splines  = n_splines,
    knots      = knots,
    spline_var = spline_var,
    dat        = dat
  ))
}

#~############################################################################~#
# Predict from splines ----
#~############################################################################~#

# Restricted cubic spline prediction function, for any number of knots.
# Matches rcs() / rcspline.eval(..., inclx = TRUE, norm = 2).
make_g_hat <- function(model, knots, coefs = coef(model)) {
  knots <- sort(knots)
  k     <- length(knots)
  stopifnot(k >= 3)

  if (length(coefs) != k)
    stop(sprintf(
      "Expected %d coefficients (intercept + %d spline terms), got %d",
      k, k - 1, length(coefs)))

  tk  <- knots[k]; tk1 <- knots[k - 1]
  tau <- (tk - knots[1])^2                 # Harrell's norm-2 scaling

  function(X) {
    X  <- as.numeric(X)
    nl <- sapply(seq_len(k - 2), function(j)
      ( pmax(X - knots[j], 0)^3
        - pmax(X - tk1, 0)^3 * (tk  - knots[j]) / (tk - tk1)
        + pmax(X - tk,  0)^3 * (tk1 - knots[j]) / (tk - tk1) ) / tau)
    nl <- matrix(nl, nrow = length(X))     # keep n x (k-2) even when n == 1
    drop(cbind(1, X, nl) %*% coefs)        # intercept + linear + nonlinear
  }
}

# Same, for spline models that carry other moderators too. They are held at the
# values in `at` (named by coefficient, e.g. c(follow_up_years = 0)), or at 0,
# and folded into the intercept before handing over to make_g_hat().
make_g_hat_baseline <- function(model, knots, spline_var,
                                at = NULL, coefs = coef(model)) {

  k <- length(knots)

  coef_names <- names(coefs)
  if (is.null(coef_names)) coef_names <- rownames(model$b)

  is_spline    <- grepl(paste0("^(rcs|ns|bs|pspline)\\(\\s*", spline_var, "\\s*[,)]"),
                        coef_names)
  is_intercept <- coef_names == "intrcpt"

  spline_coefs <- coefs[is_spline]
  if (length(spline_coefs) != k - 1) {
    stop(sprintf(
      "Expected %d spline coefficients for %s (%d knots), got %d: %s",
      k - 1, spline_var, k, length(spline_coefs),
      paste(coef_names[is_spline], collapse = ", ")))
  }
  if (!any(is_intercept)) stop("Model has no intercept ('intrcpt') coefficient.")

  other_coefs <- coefs[!is_spline & !is_intercept]
  intercept   <- unname(coefs[is_intercept])

  if (length(other_coefs) > 0) {
    at_values <- setNames(rep(0, length(other_coefs)), names(other_coefs))

    if (!is.null(at)) {
      unknown <- setdiff(names(at), names(at_values))
      if (length(unknown) > 0) {
        stop("`at` names not found among the model coefficients: ",
             paste(unknown, collapse = ", "), ".\n  Available: ",
             paste(names(at_values), collapse = ", "))
      }
      at_values[names(at)] <- at
    }

    defaulted <- setdiff(names(at_values), names(at))
    if (length(defaulted) > 0) {
      message("make_g_hat_baseline(): holding at 0: ",
              paste(defaulted, collapse = ", "))
    }

    intercept <- intercept + sum(other_coefs * at_values)
  }

  make_g_hat(model, knots, coefs = c(intercept, unname(spline_coefs)))
}
