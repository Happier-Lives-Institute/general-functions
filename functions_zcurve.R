#~############################################################################~#
# z_dens ----
#~############################################################################~#
# modify the zcurve function to get the density so we can do nice ggplots
z_dens <- function (x, annotation = FALSE, CI = FALSE, extrapolate = FALSE,
                    y.anno = c(0.95, 0.88, 0.78, 0.71, 0.61, 0.53, 0.43, 0.35),
                    x.anno = 0.6, cex.anno = 1, ...)
{
  if (is.null(x$boot))
    CI <- FALSE
  additional <- list(...)
  if (is.null(additional$main)) {
    main <- paste(c("z-curve (", ifelse(is.null(x$control$model),
                                        "custom", x$control$model), " via ", x$method, ")"),
                  collapse = "")
  }
  else {
    main <- additional$main
  }
  if (is.null(additional$xlab)) {
    xlab <- "z-scores"
  }
  else {
    xlab <- additional$xlab
  }
  if (is.null(additional$ylab)) {
    ylab <- "Density"
  }
  else {
    ylab <- additional$ylab
  }
  if (is.null(additional$xlim)) {
    xlim <- NULL
  }
  else {
    xlim <- additional$xlim
  }
  if (is.null(additional$ylim)) {
    ylim <- NULL
  }
  else {
    ylim <- additional$ylim
  }
  if (is.null(additional$cex.axis)) {
    cex.axis <- 1
  }
  else {
    cex.axis <- additional$cex.axis
  }
  if (is.null(additional$cex.lab)) {
    cex.lab <- 1
  }
  else {
    cex.lab <- additional$cex.lab
  }
  br1 <- seq(x$control$a, x$control$b, 0.2)
  br2 <- seq(0, x$control$a, 0.2)
  br1[length(br1)] <- x$control$b
  br2[length(br2)] <- x$control$a
  h1 <- graphics::hist(x$data[x$data > x$control$a & x$data <
                                x$control$b], breaks = br1, plot = F)
  if (nrow(x$data_censoring) != 0) {
    cen_counts <- do.call(rbind, lapply(1:nrow(x$data_censoring),
                                        function(i) {
                                          temp_counts <- (x$data_censoring$lb[i] < h1$breaks[-length(h1$breaks)] &
                                                            x$data_censoring$ub[i] > h1$breaks[-length(h1$breaks)])
                                          temp_counts[temp_counts] <- 1/sum(temp_counts)
                                          return(temp_counts)
                                        }))
    cen_counts <- apply(cen_counts, 2, sum)
    h1$counts <- h1$counts + cen_counts
    h1$density <- (h1$counts/sum(h1$counts)) * (length(h1$counts)/(x$control$b -
                                                                     x$control$a))
  }
  if (length(x$data[x$data < x$control$a])) {
    h2 <- graphics::hist(x$data[x$data < x$control$a], breaks = br2,
                         plot = F)
    h2$density <- h2$density * (x$control$a/(x$control$b -
                                               x$control$a))
    h2$density <- h2$density/((length(x$data[x$data > x$control$a &
                                               x$data < x$control$b])/(x$control$b - x$control$a))/(length(x$data[x$data <
                                                                                                                    x$control$a])/(x$control$a)))
  }
  else {
    h2 <- NULL
  }
  x_seq <- seq(0, x$control$b, 0.01)
  y_den <- sapply(1:length(x$fit$mu), function(i) {
    x$fit$weights[i] * exp(zcurve:::.zdist_lpdf(x_seq, x$fit$mu[i],
                                                1, x$control$a, x$control$b))
  })
  y_den <- apply(y_den, 1, sum)
  if (CI & !is.null(x$boot)) {
    # Every bootstrap draw uses the same handful of mu values, so we work
    # out the density for each of those once instead of once per draw.
    # rowSums adds in the same order as apply(, 1, sum), so this gives
    # exactly the same numbers.
    mu_u <- unique(as.vector(x$boot$mu))
    D <- vapply(mu_u, function(mm) {
      exp(zcurve:::.zdist_lpdf(x_seq, mm, 1, x$control$a, x$control$b))
    }, numeric(length(x_seq)))

    idx <- matrix(match(x$boot$mu, mu_u), nrow = nrow(x$boot$mu))
    n_x <- length(x_seq)

    y_den_boot <- vapply(seq_len(nrow(x$boot$mu)), function(b) {
      rowSums(D[, idx[b, ], drop = FALSE] * rep(x$boot$weights[b, ], each = n_x))
    }, numeric(n_x))
    y_den_l.CI <- apply(y_den_boot, 1, stats::quantile,
                        prob = 0.025)
    y_den_u.CI <- apply(y_den_boot, 1, stats::quantile,
                        prob = 0.975)
  }
  x_max <- x$control$b
  x.anno <- x_max * x.anno
  if (extrapolate) {
    y_max <- max(c(y_den, h1$density, h2$density))
  }
  else {
    y_max <- max(c(h1$density, h2$density))
  }
  if (annotation & CI) {
    y_max <- ifelse(max(c(y_den_u.CI[x.anno < x_seq], h1$density[x.anno <
                                                                   h1$breaks[-length(h1$density)]])) > y_max * (y.anno[length(y.anno)] -
                                                                                                                 0.025), max(c(y_den_u.CI[x.anno < x_seq], h1$density[x.anno <
                                                                                                                                                                       h1$breaks[-length(h1$density)]]))/(y.anno[length(y.anno)] -
                                                                                                                                                                                                            0.025), y_max)
  }
  else if (annotation & !CI) {
    y_max <- ifelse(max(c(y_den[x.anno < x_seq], h1$density[x.anno <
                                                              h1$breaks[-length(h1$density)]])) > y_max * (y.anno[length(y.anno)] -
                                                                                                             0.025), max(c(y_den[x.anno < x_seq], h1$density[x.anno <
                                                                                                                                                              h1$breaks[-length(h1$density)]]))/(y.anno[length(y.anno)] -
                                                                                                                                                                                                   0.025), y_max)
  }
  if (!is.null(ylim)) {
    y_min <- ylim[1]
    y_max <- ylim[2]
  }
  else {
    y_min <- 0
  }
  if (!is.null(xlim)) {
    x_min <- xlim[1]
    x_max <- xlim[2]
  }
  else {
    x_min <- 0
  }

  h1_df <- data.frame(
    xmin = h1$breaks[-length(h1$breaks)],
    xmax = h1$breaks[-1],
    density = h1$density
  )
  h2_df <- if (is.null(h2)) {
    data.frame(xmin = numeric(0), xmax = numeric(0), density = numeric(0))
  } else {
    data.frame(
      xmin = h2$breaks[-length(h2$breaks)],
      xmax = h2$breaks[-1],
      density = h2$density
    )
  }

  return(list(
    curve = data.frame(
      x_seq = x_seq,
      y_den = y_den,
      y_den_l.CI = if (CI && !is.null(x$boot)) y_den_l.CI else NA,
      y_den_u.CI = if (CI && !is.null(x$boot)) y_den_u.CI else NA
    ),
    hist1 = h1_df,
    hist2 = h2_df,
    a = x$control$a,
    b = x$control$b
  ))
}

#~############################################################################~#
# zcurve_bias_tests ----
#~############################################################################~#
zcurve_bias_tests <- function(zc, z_values, just_range = c(1.96, 2.58),
                              strong_cutoff = 2.58) {

  stopifnot(inherits(zc, "zcurve"))
  if (missing(z_values))
    stop("You must provide the original z-values used to fit the model.")

  mix_prob <- function(zc, lower, upper) {
    mu <- zc$fit$mu
    w  <- zc$fit$weights
    a  <- zc$control$a
    b  <- zc$control$b

    f <- function(x) {
      rowSums(sapply(seq_along(mu), function(i) {
        w[i] * exp(zcurve:::.zdist_lpdf(x, mu[i], 1, a, b))
      }))
    }
    integrate(f, lower, upper)$value
  }

  k_all <- length(z_values)
  k_sig <- sum(z_values >= zc$control$a)
  odr <- k_sig / k_all
  odr_ci <- prop.test(k_sig, k_all)$conf.int
  odr_lci <- odr_ci[1]; odr_uci <- odr_ci[2]

  edr <- summary(zc)$coefficients["EDR", "Estimate"]
  edr_lci <- summary(zc)$coefficients["EDR", "l.CI"]
  edr_uci <- summary(zc)$coefficients["EDR", "u.CI"]

  odr_exceeds_edr <- if (is.na(edr_uci)) NA else odr > edr_uci

  obs_just <- mean(z_values >= just_range[1] & z_values < just_range[2])
  pred_just <- mix_prob(zc, just_range[1], just_range[2])
  ejs_test <- binom.test(obs_just * k_all, k_all, pred_just, alternative = "greater")

  fit_strong <- zcurve::zcurve(
    z_values[z_values > strong_cutoff],
    bootstrap = FALSE,
    control = list(a = strong_cutoff, b = zc$control$b)
  )
  pred_just_from_strong <- mix_prob(fit_strong, just_range[1], just_range[2])
  phack_test <- binom.test(obs_just * k_all, k_all, pred_just_from_strong, alternative = "greater")

  results <- rbind(
    data.frame(
      Test = "ODR",
      Estimate = round_per(odr),
      LCI = round_per(odr_lci),
      UCI = round_per(odr_uci),
      Observed = round_per(odr),
      Predicted = "",
      p_value = "",
      Notes = sprintf("n_sig=%d of n=%d; binomial CI", k_sig, k_all),
      row.names = NULL
    ),
    data.frame(
      Test = "EDR",
      Estimate = round_per(edr),
      LCI = round_per(edr_lci),
      UCI = round_per(edr_uci),
      Observed = "",
      Predicted = round_per(edr),
      p_value = "",
      Notes = if (is.na(edr_lci)) "no bootstrap CI" else "bootstrap CI",
      row.names = NULL
    ),
    data.frame(
      Test = "ODR > EDR upper CI",
      Estimate = if (is.na(odr_exceeds_edr)) NA else odr_exceeds_edr,
      LCI = "",
      UCI = "",
      Observed = round_per(odr),
      Predicted = if (is.na(edr_uci)) "" else round_per(edr_uci),
      p_value = "",
      Notes = if (is.na(odr_exceeds_edr)) "EDR CI unavailable"
      else if (odr_exceeds_edr) "evidence of selection bias"
      else "no excess over EDR CI",
      row.names = NULL
    ),
    data.frame(
      Test = "Excess Just Significant (binomial)",
      Estimate = round_per(obs_just / pred_just),
      LCI = "",
      UCI = "",
      Observed = round_per(obs_just),
      Predicted = round_per(pred_just),
      p_value = signif(ejs_test$p.value, 3),
      Notes = sprintf("window %.2f to %.2f", just_range[1], just_range[2]),
      row.names = NULL
    ),
    data.frame(
      Test = sprintf("P hacking (binomial; fit z>%s)", strong_cutoff),
      Estimate = round_per(obs_just / pred_just_from_strong),
      LCI = "",
      UCI = "",
      Observed = round_per(obs_just),
      Predicted = round_per(pred_just_from_strong),
      p_value = signif(phack_test$p.value, 3),
      Notes = sprintf("predicting %.2f to %.2f from z>%s",
                      just_range[1], just_range[2], strong_cutoff),
      row.names = NULL
    )
  )

  return(results)
}
