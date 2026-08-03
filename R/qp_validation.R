# Exact numerical solver and independent-validation helpers for population QPs.

# Solve a small, strictly convex box-constrained quadratic programme by exact
# active-set enumeration. Impossible states are excluded before enumeration,
# so a lower-bound-only p-trait problem has 2^p states rather than 3^p.
.dgr_box_qp <- function(H, f, lower, upper, tolerance = 1e-9) {
  H <- as.matrix(H)
  f <- as.numeric(f)
  lower <- as.numeric(lower)
  upper <- as.numeric(upper)
  p <- length(f)
  if (!identical(dim(H), c(p, p)) || length(lower) != p ||
    length(upper) != p || any(!is.finite(H)) || any(!is.finite(f)) ||
    any(is.na(lower)) || any(is.na(upper)) || any(lower > upper)) {
    stop("Invalid box-constrained quadratic programme.", call. = FALSE)
  }
  H <- (H + t(H)) / 2
  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= 0) {
    stop("H must be strictly positive definite.", call. = FALSE)
  }
  fixed <- is.finite(lower) & is.finite(upper) &
    abs(upper - lower) <= tolerance
  allowed_states <- lapply(seq_len(p), function(j) {
    if (fixed[j]) {
      return(2L)
    }
    states <- 0L
    if (is.finite(lower[j])) states <- c(states, -1L)
    if (is.finite(upper[j])) states <- c(states, 1L)
    states
  })
  state_counts <- vapply(allowed_states, length, integer(1L))
  n_states <- prod(state_counts)
  if (!is.finite(n_states) || n_states > 1e6) {
    stop("The exact active-set problem has more than one million states; ",
      "reduce its dimension or use an external convex solver.",
      call. = FALSE
    )
  }
  scale <- max(
    1, max(abs(H)), max(abs(f)),
    max(abs(c(lower[is.finite(lower)], upper[is.finite(upper)])), 0)
  )
  tol <- tolerance * scale
  best <- NULL
  best_objective <- Inf

  for (state_id in 0:(n_states - 1L)) {
    states <- integer(p)
    code <- state_id
    for (j in seq_len(p)) {
      position <- code %% state_counts[j] + 1L
      states[j] <- allowed_states[[j]][position]
      code <- code %/% state_counts[j]
    }
    at_lower <- states == -1L
    free <- states == 0L
    at_upper <- states == 1L
    if (any(at_lower & !is.finite(lower)) ||
      any(at_upper & !is.finite(upper))) {
      next
    }
    x <- numeric(p)
    x[at_lower] <- lower[at_lower]
    x[at_upper] <- upper[at_upper]
    x[fixed] <- lower[fixed]
    active <- !free
    if (any(free)) {
      rhs <- -f[free]
      if (any(active)) {
        rhs <- rhs - H[free, active, drop = FALSE] %*% x[active]
      }
      x[free] <- tryCatch(
        as.numeric(solve(H[free, free, drop = FALSE], rhs)),
        error = function(e) rep(NA_real_, sum(free))
      )
    }
    if (any(!is.finite(x)) || any(x < lower - tol) || any(x > upper + tol)) {
      next
    }
    gradient <- as.numeric(H %*% x + f)
    if (any(free) && max(abs(gradient[free])) > 50 * tol) next
    if (any(at_lower) && min(gradient[at_lower]) < -50 * tol) next
    if (any(at_upper) && max(gradient[at_upper]) > 50 * tol) next
    objective <- as.numeric(
      0.5 * crossprod(x, H %*% x) + crossprod(f, x)
    )
    if (objective < best_objective) {
      best <- x
      best_objective <- objective
    }
  }
  if (is.null(best)) {
    stop("No KKT-feasible active set was found.", call. = FALSE)
  }
  kkt <- .dgr_box_kkt(best, H, f, lower, upper, tolerance)
  list(
    solution = best,
    objective = best_objective,
    kkt = kkt,
    status = if (kkt$overall <= 100 * tolerance * scale) {
      "solved"
    } else {
      "numerical_warning"
    },
    eigenvalues = eigenvalues
  )
}

.dgr_box_kkt <- function(x, H, f, lower, upper, tolerance = 1e-9) {
  x <- as.numeric(x)
  gradient <- as.numeric(H %*% x + f)
  scale <- max(1, max(abs(H)), max(abs(f)), max(abs(x)))
  bound_tolerance <- max(tolerance, 20 * .Machine$double.eps) * scale
  fixed <- is.finite(lower) & is.finite(upper) &
    abs(upper - lower) <= bound_tolerance
  at_lower <- !fixed & is.finite(lower) &
    abs(x - lower) <= bound_tolerance
  at_upper <- !fixed & is.finite(upper) &
    abs(x - upper) <= bound_tolerance
  interior <- !fixed & !at_lower & !at_upper
  stationarity <- numeric(length(x))
  stationarity[interior] <- abs(gradient[interior])
  stationarity[at_lower] <- pmax(0, -gradient[at_lower])
  stationarity[at_upper] <- pmax(0, gradient[at_upper])
  primal <- max(c(
    0,
    lower[is.finite(lower)] - x[is.finite(lower)],
    x[is.finite(upper)] - upper[is.finite(upper)]
  ))
  list(
    stationarity = max(stationarity),
    primal = primal,
    overall = max(max(stationarity), primal),
    gradient = gradient
  )
}
