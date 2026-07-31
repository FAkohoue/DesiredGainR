# A compact anisotropic Gaussian-process regressor used by
# optimize_desired_gains().
#
# This is deliberately self-contained rather than a dependency. The package
# needs a surrogate over at most a few dimensions with a few hundred noisy
# observations, which is well within what a few hundred lines can do reliably,
# and adding a modelling dependency for that would be disproportionate.
#
# The prior mean is constant or linear. An earlier design used a deterministic
# infinitesimal-model recursion as the prior mean; that was rejected, because
# it would embed a second genetic model whose assumptions could be wrong in
# precisely the situations where simulation is most needed.

#' Matern 5/2 correlation matrix
#'
#' @param X1,X2 Numeric matrices with observations in rows.
#' @param lengthscale Positive per-dimension lengthscales.
#'
#' @return The correlation matrix between the rows of `X1` and `X2`.
#' @noRd
.dgr_matern52 <- function(X1, X2, lengthscale) {
  scaled1 <- sweep(X1, 2L, lengthscale, "/")
  scaled2 <- sweep(X2, 2L, lengthscale, "/")
  squared <- outer(rowSums(scaled1^2), rowSums(scaled2^2), "+") -
    2 * tcrossprod(scaled1, scaled2)
  distance <- sqrt(pmax(squared, 0))
  root5 <- sqrt(5)
  (1 + root5 * distance + (5 / 3) * distance^2) * exp(-root5 * distance)
}

#' Design matrix for the prior mean
#'
#' @param X Numeric matrix with observations in rows.
#' @param prior_mean Either `"constant"` or `"linear"`.
#'
#' @return A basis matrix.
#' @noRd
.dgr_mean_basis <- function(X, prior_mean) {
  if (identical(prior_mean, "linear")) cbind(1, X) else matrix(1, nrow(X), 1L)
}

#' Negative log marginal likelihood of the Gaussian process
#'
#' @param parameters Log-scale hyperparameters: lengthscales, then the noise to
#'   signal ratio.
#' @param X,y Training inputs and responses.
#' @param H Prior-mean basis.
#'
#' @return The negative log marginal likelihood, or a large finite penalty when
#'   the covariance is not decomposable.
#' @noRd
.dgr_gp_nlml <- function(parameters, X, y, H) {
  d <- ncol(X)
  lengthscale <- exp(parameters[seq_len(d)])
  nugget_ratio <- exp(parameters[d + 1L])
  n <- nrow(X)

  R <- .dgr_matern52(X, X, lengthscale)
  diag(R) <- diag(R) + nugget_ratio + 1e-8
  chol_R <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(chol_R)) return(1e10)

  solve_R <- function(z) backsolve(chol_R, backsolve(chol_R, z, transpose = TRUE))
  R_inv_y <- solve_R(y)
  R_inv_H <- solve_R(H)
  HtRinvH <- crossprod(H, R_inv_H)
  chol_H <- tryCatch(chol(HtRinvH), error = function(e) NULL)
  if (is.null(chol_H)) return(1e10)

  beta <- backsolve(chol_H, backsolve(chol_H, crossprod(H, R_inv_y),
                                      transpose = TRUE))
  residual <- y - H %*% beta
  R_inv_residual <- solve_R(residual)
  # Profile out the signal variance analytically.
  sigma2 <- as.numeric(crossprod(residual, R_inv_residual)) / n
  if (!is.finite(sigma2) || sigma2 <= 0) return(1e10)

  log_determinant <- 2 * sum(log(diag(chol_R)))
  0.5 * (n * log(sigma2) + log_determinant)
}

#' Fit a Gaussian process
#'
#' @param X Numeric matrix of inputs, observations in rows.
#' @param y Numeric response vector.
#' @param prior_mean Either `"constant"` or `"linear"`.
#' @param n_restarts Number of random restarts of the hyperparameter search.
#'
#' @return A fitted model object used by `.dgr_gp_predict()`.
#' @noRd
.dgr_gp_fit <- function(X, y, prior_mean = "constant", n_restarts = 5L) {
  X <- as.matrix(X)
  d <- ncol(X)
  n <- nrow(X)
  # Standardise the response so that the hyperparameter scales are comparable
  # between problems.
  y_centre <- mean(y)
  y_scale <- stats::sd(y)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  y_standard <- (y - y_centre) / y_scale
  H <- .dgr_mean_basis(X, prior_mean)

  best <- NULL
  best_value <- Inf
  starts <- c(
    list(c(rep(log(0.3), d), log(1e-4))),
    lapply(seq_len(max(0L, n_restarts - 1L)), function(i) {
      c(stats::runif(d, log(0.05), log(2)), stats::runif(1, log(1e-6), log(0.5)))
    })
  )
  for (start in starts) {
    fit <- tryCatch(
      stats::optim(
        start, .dgr_gp_nlml, X = X, y = y_standard, H = H,
        method = "L-BFGS-B",
        lower = c(rep(log(0.02), d), log(1e-8)),
        upper = c(rep(log(10), d), log(2)),
        control = list(maxit = 200L)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < best_value) {
      best_value <- fit$value
      best <- fit
    }
  }
  if (is.null(best)) {
    best <- list(par = c(rep(log(0.3), d), log(1e-3)))
  }

  lengthscale <- exp(best$par[seq_len(d)])
  nugget_ratio <- exp(best$par[d + 1L])
  R <- .dgr_matern52(X, X, lengthscale)
  diag(R) <- diag(R) + nugget_ratio + 1e-8
  chol_R <- chol(R)
  solve_R <- function(z) backsolve(chol_R, backsolve(chol_R, z, transpose = TRUE))
  R_inv_H <- solve_R(H)
  HtRinvH <- crossprod(H, R_inv_H)
  chol_H <- chol(HtRinvH)
  beta <- backsolve(chol_H, backsolve(chol_H, crossprod(H, solve_R(y_standard)),
                                      transpose = TRUE))
  residual <- y_standard - H %*% beta
  R_inv_residual <- solve_R(residual)
  sigma2 <- as.numeric(crossprod(residual, R_inv_residual)) / n

  structure(
    list(
      X = X, y = y, prior_mean = prior_mean,
      y_centre = y_centre, y_scale = y_scale,
      lengthscale = lengthscale, nugget_ratio = nugget_ratio,
      sigma2 = sigma2, beta = beta,
      chol_R = chol_R, chol_H = chol_H,
      R_inv_residual = R_inv_residual, R_inv_H = R_inv_H,
      nlml = best_value
    ),
    class = "desiredgainr_gp"
  )
}

#' Predict from a fitted Gaussian process
#'
#' @param model A model from `.dgr_gp_fit()`.
#' @param X_new Numeric matrix of query points.
#'
#' @return A list with the posterior mean and standard deviation on the
#'   original response scale.
#' @noRd
.dgr_gp_predict <- function(model, X_new) {
  X_new <- as.matrix(X_new)
  k <- .dgr_matern52(X_new, model$X, model$lengthscale)
  H_new <- .dgr_mean_basis(X_new, model$prior_mean)

  mean_standard <- as.numeric(
    H_new %*% model$beta + k %*% model$R_inv_residual
  )
  solve_R <- function(z) {
    backsolve(model$chol_R, backsolve(model$chol_R, z, transpose = TRUE))
  }
  R_inv_k <- solve_R(t(k))
  quadratic <- colSums(t(k) * R_inv_k)
  # Universal-kriging correction for uncertainty in the prior-mean
  # coefficients.
  u <- t(H_new) - crossprod(model$R_inv_H, t(k))
  correction <- colSums(
    u * backsolve(model$chol_H, backsolve(model$chol_H, u, transpose = TRUE))
  )
  variance <- model$sigma2 * pmax(1 - quadratic + correction, 0)

  list(
    mean = mean_standard * model$y_scale + model$y_centre,
    sd = sqrt(variance) * model$y_scale
  )
}

#' Expected improvement for a maximisation problem
#'
#' @param mean,sd Posterior mean and standard deviation.
#' @param best Best value observed so far.
#' @param xi Exploration margin.
#'
#' @return The expected improvement at each point.
#' @noRd
.dgr_expected_improvement <- function(mean, sd, best, xi = 0.01) {
  improvement <- mean - best - xi
  result <- numeric(length(mean))
  positive <- sd > .Machine$double.eps
  z <- improvement[positive] / sd[positive]
  result[positive] <- improvement[positive] * stats::pnorm(z) +
    sd[positive] * stats::dnorm(z)
  pmax(result, 0)
}

#' Probability that a quantity meets or exceeds a floor
#'
#' @param mean,sd Posterior mean and standard deviation.
#' @param floor Required minimum.
#'
#' @return The probability of satisfying the constraint at each point.
#' @noRd
.dgr_probability_feasible <- function(mean, sd, floor) {
  result <- as.numeric(mean >= floor)
  positive <- sd > .Machine$double.eps
  result[positive] <- stats::pnorm((mean[positive] - floor) / sd[positive])
  result
}
