# A version of limitmeta that does not crash
limitmeta_no_crash <- function (x, method.adjust = "beta0", level = x$level, level.ma = x$level.ma,
                                backtransf = x$backtransf, title = x$title, complab = x$complab,
                                outclab = x$outclab)
{
  metasens:::chkclass(x, "meta")
  method.adjust <- metasens:::setchar(method.adjust, c("beta0", "betalim", "mulim"))
  TE <- x$TE
  seTE <- x$seTE
  tau <- x$tau
  w.random <- x$w.random
  k <- x$k
  Q <- x$Q
  sm <- x$sm
  seTE.tau <- sqrt(1 / w.random)
  TE.random <- x$TE.random
  seTE.random <- x$seTE.random
  lower.random <- x$lower.random
  upper.random <- x$upper.random
  statistic.random <- x$statistic.random
  pval.random <- x$pval.random
  reg.f <- metasens:::radialregression(TE, seTE, k)
  reg.r <- metasens:::radialregression(TE, seTE.tau, k)
  alpha.r <- reg.r$intercept
  beta.r <- reg.r$slope
  TE.limit <- beta.r + sqrt(tau^2 / seTE.tau^2) * (TE - beta.r)
  seTE.limit <- seTE / 1
  # m.lim <- metagen(TE.limit, seTE.limit, sm = sm, method.tau.ci = "")
  # reg.l <- metasens:::radialregression(m.lim$TE, m.lim$seTE, k)
  reg.l <- metasens:::radialregression(TE.limit, seTE.limit, k)

  if (method.adjust == "beta0") {
    TE.adjust <- as.vector(beta.r + tau * alpha.r)
    var.beta <- as.vector(1 / var(sqrt(w.random)) / (k - 1))
    var.alpha <- as.vector(mean(w.random) / var(sqrt(w.random)) / (k - 1))
    cov.alpha.beta <- as.vector(-mean(sqrt(w.random)) / var(sqrt(w.random)) / (k - 1))
    seTE.adjust <- sqrt(var.beta + tau^2 * var.alpha + 2 * tau * cov.alpha.beta)
  } else {
    if (method.adjust == "mulim") {
      m.lim <- metagen(TE.limit, seTE.limit, sm = sm, method.tau.ci = "")
      TE.adjust <- m.lim$TE.common
      seTE.adjust <- m.lim$seTE.common
    } else if (method.adjust == "betalim") {
      TE.adjust <- as.vector(reg.l$slope)
      seTE.adjust <- as.vector(reg.l$se.slope)
    }
  }

  ci.adjust <- ci(TE.adjust, seTE.adjust, level = level.ma)
  lower.adjust <- ci.adjust$lower
  upper.adjust <- ci.adjust$upper
  statistic.adjust <- ci.adjust$statistic
  pval.adjust <- ci.adjust$p

  if (inherits(x, c("metaprop"))) {
    statistic.adjust <- NA
    pval.adjust <- NA
  }

  if (!missing(level.ma)) {
    ci.r <- metasens:::ci(TE.random, seTE.random, level = level.ma)
    lower.random <- ci.r$lower
    upper.random <- ci.r$upper
  }

  Q.resid <- reg.f$sigma^2 * (k - 2)
  Q.small <- Q - Q.resid
  G.squared <- 1 - reg.l$r.squared

  res <- list(
    TE = TE, seTE = seTE, TE.limit = TE.limit, seTE.limit = seTE.limit,
    studlab = x$studlab, TE.random = TE.random, seTE.random = seTE.random,
    lower.random = lower.random, upper.random = upper.random,
    statistic.random = statistic.random, pval.random = pval.random,
    w.random = w.random, TE.adjust = TE.adjust, seTE.adjust = seTE.adjust,
    lower.adjust = lower.adjust, upper.adjust = upper.adjust,
    statistic.adjust = statistic.adjust, pval.adjust = pval.adjust,
    alpha.r = alpha.r, beta.r = beta.r, Q = Q, df.Q = x$df.Q,
    Q.small = Q.small, Q.resid = Q.resid, G.squared = G.squared,
    tau = tau, level = level, level.ma = level.ma, k = k,
    sm = sm, method.adjust = method.adjust, title = title,
    complab = complab, outclab = outclab, call = match.call(), x = x
  )

  res$backtransf <- backtransf
  res$version <- utils::packageDescription("metasens")$Version
  class(res) <- "limitmeta"
  res

}
