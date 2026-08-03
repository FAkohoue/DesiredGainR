# Independent validation of DesiredGainR's active-set quadratic solver.
#
# Run from the package root. Optional environment variables:
#   DESIREDGAINR_PYTHON        Python executable
#   DESIREDGAINR_QP_PYTHONPATH Directory containing clarabel/osqp/scipy
#   DESIREDGAINR_QP_REPORT     Output CSV path

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for solver validation.", call. = FALSE)
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Package 'devtools' is required for solver validation.", call. = FALSE)
}
devtools::load_all(".", quiet = TRUE)

encode_bound <- function(x) {
  vapply(x, function(value) {
    if (is.infinite(value) && value > 0) {
      "Inf"
    } else if (is.infinite(value)) {
      "-Inf"
    } else {
      format(value, digits = 17, scientific = TRUE, trim = TRUE)
    }
  }, character(1L))
}

make_problem <- function(id, H, f, lower, upper, class, seed,
                         expected = NULL) {
  list(
    problem_id = id,
    H = unname(H), f = as.numeric(f),
    lower = encode_bound(lower), upper = encode_bound(upper),
    class = class, seed = seed,
    condition_number = kappa(H, exact = TRUE), expected = expected
  )
}

problems <- list()
add_problem <- function(...) {
  problems[[length(problems) + 1L]] <<- make_problem(...)
}

# Exact analytical cases: diagonal clipping and an unconstrained interior.
H <- diag(c(1, 2, 4, 8))
f <- c(-3, 4, -2, 16)
lower <- c(-1, -1, 0, -Inf)
upper <- c(2, 3, 0.4, 1)
expected <- pmin(upper, pmax(lower, -f / diag(H)))
add_problem(
  "analytic_diagonal", H, f, lower, upper, "analytical", 1L,
  expected
)
H <- matrix(c(3, 0.4, 0.4, 2), 2)
f <- c(-0.4, 0.7)
expected <- as.numeric(solve(H, -f))
add_problem(
  "analytic_interior", H, f, rep(-Inf, 2), rep(Inf, 2),
  "analytical", 2L, expected
)
add_problem(
  "analytic_one_dimensional", matrix(3, 1), -12, -1, 2,
  "analytical", 3L, 2
)
add_problem(
  "analytic_narrow_interval", diag(2), c(-1, 1),
  c(0.2, -0.300000000001), c(0.200000000001, -0.3),
  "analytical", 4L, c(0.200000000001, -0.300000000001)
)

# Constructed active sets, including fixed, infinite and narrow bounds.
for (seed in 101:120) {
  set.seed(seed)
  p <- sample(2:6, 1L)
  M <- matrix(stats::rnorm(p^2), p)
  H <- crossprod(M) + diag(p)
  lower <- stats::runif(p, -1.5, -0.2)
  upper <- stats::runif(p, 0.2, 1.5)
  state <- sample(c("lower", "free", "upper"), p, replace = TRUE)
  if (seed %% 4L == 0L) {
    fixed <- stats::runif(1L, lower[1L], upper[1L])
    lower[1L] <- upper[1L] <- fixed
    state[1L] <- "fixed"
  }
  if (seed %% 5L == 0L && p > 2L) {
    lower[2L] <- -Inf
    upper[3L] <- Inf
    if (state[2L] == "lower") state[2L] <- "free"
    if (state[3L] == "upper") state[3L] <- "free"
  }
  x <- vapply(seq_len(p), function(j) {
    lo <- if (is.finite(lower[j])) lower[j] + 0.1 else -0.5
    hi <- if (is.finite(upper[j])) upper[j] - 0.1 else 0.5
    if (lo > hi) lower[j] else stats::runif(1L, lo, hi)
  }, numeric(1L))
  gradient <- numeric(p)
  x[state == "fixed"] <- lower[state == "fixed"]
  gradient[state == "fixed"] <- stats::runif(sum(state == "fixed"), -2, 2)
  x[state == "lower"] <- lower[state == "lower"]
  gradient[state == "lower"] <- stats::runif(sum(state == "lower"), 0.1, 2)
  x[state == "upper"] <- upper[state == "upper"]
  gradient[state == "upper"] <- -stats::runif(sum(state == "upper"), 0.1, 2)
  f <- gradient - as.numeric(H %*% x)
  add_problem(
    sprintf("active_%03d", seed), H, f, lower, upper,
    "constructed_active", seed, x
  )
}

# Controlled conditioning, including cases where objective/KKT agreement is
# more informative than coordinate-wise agreement.
condition_targets <- c(1, 1e2, 1e4, 1e6, 1e8, 1e10)
for (condition_target in condition_targets) {
  for (replicate in 1:3) {
    seed <- as.integer(1000 + log10(condition_target + 1) * 10 + replicate)
    set.seed(seed)
    p <- 5L
    Q <- qr.Q(qr(matrix(stats::rnorm(p^2), p)))
    # Keep the largest eigenvalue at one. This preserves the requested
    # condition number without confounding conditioning with arbitrary
    # objective magnitude. The known optimum is strictly interior.
    eigenvalues <- exp(seq(-log(condition_target), 0, length.out = p))
    H <- Q %*% diag(eigenvalues) %*% t(Q)
    expected <- stats::runif(p, -0.5, 0.5)
    f <- -as.numeric(H %*% expected)
    lower <- rep(-1, p)
    upper <- rep(1, p)
    add_problem(
      sprintf("condition_%010.0f_%d", condition_target, replicate),
      H, f, lower, upper,
      if (condition_target <= 1e4) "well_conditioned" else "conditioned",
      seed, expected
    )
  }
}

# Numerically delicate distances from a bound.
for (distance in c(0, 1e-12, 1e-10, 1e-8)) {
  H <- matrix(c(2, 0.3, 0.3, 1.5), 2)
  lower <- c(0, -1)
  upper <- c(1, 1)
  x <- c(distance, 0.2)
  gradient <- c(if (distance == 0) 0.4 else 0, 0)
  f <- gradient - as.numeric(H %*% x)
  add_problem(
    paste0("near_bound_", format(distance, scientific = TRUE)),
    H, f, lower, upper, "degenerate_active_set", 2001L, x
  )
}

# Invariance cases derived without calling package transformation helpers.
base <- problems[[which(vapply(
  problems, `[[`, character(1L), "problem_id"
) == "active_101")]]
H <- as.matrix(base$H)
f <- unlist(base$f)
lower <- as.numeric(base$lower)
upper <- as.numeric(base$upper)
p <- length(f)
permutation <- rev(seq_len(p))
add_problem(
  "invariance_permutation", H[permutation, permutation],
  f[permutation], lower[permutation], upper[permutation], "invariance", 3001L,
  unlist(base$expected)[permutation]
)
scales <- seq(0.5, 2, length.out = p)
# x = D y, so H_y = D H D and f_y = D f.
D <- diag(scales)
add_problem(
  "invariance_rescaling", D %*% H %*% D, as.numeric(D %*% f),
  lower / scales, upper / scales, "invariance", 3002L,
  unlist(base$expected) / scales
)
signs <- rep(c(-1, 1), length.out = p)
S <- diag(signs)
sign_lower <- ifelse(signs > 0, lower, -upper)
sign_upper <- ifelse(signs > 0, upper, -lower)
add_problem(
  "invariance_sign", S %*% H %*% S, as.numeric(S %*% f),
  sign_lower, sign_upper, "invariance", 3003L,
  signs * unlist(base$expected)
)
add_problem(
  "invariance_objective_scale", 17 * H, 17 * f,
  lower, upper, "invariance", 3004L, unlist(base$expected)
)

input <- tempfile(fileext = ".json")
output <- tempfile(fileext = ".json")
on.exit(unlink(c(input, output), force = TRUE), add = TRUE)
jsonlite::write_json(list(problems = lapply(problems, function(problem) {
  problem[c("problem_id", "H", "f", "lower", "upper")]
})), input, auto_unbox = TRUE, digits = 17, pretty = TRUE)

python <- Sys.getenv("DESIREDGAINR_PYTHON", unset = "")

if (!nzchar(python)) {
  python <- Sys.which("python")
}

if (!nzchar(python)) {
  stop(
    paste(
      "Python was not found.",
      "Set DESIREDGAINR_PYTHON to the interpreter containing",
      "clarabel, osqp, numpy and scipy."
    ),
    call. = FALSE
  )
}

# Preserve the exact executable path supplied by setup-python.
# Do not use normalizePath(), because it dereferences the action-managed
# Python launcher/symlink.
python <- path.expand(python)

if (!file.exists(python)) {
  stop(
    "DESIREDGAINR_PYTHON points to a nonexistent executable: ",
    python,
    call. = FALSE
  )
}

oracle <- file.path(
  "inst",
  "validation",
  "qp_oracle.py"
)

if (!file.exists(oracle)) {
  stop(
    "The Python QP oracle was not found: ",
    oracle,
    call. = FALSE
  )
}

pythonpath <- Sys.getenv(
  "DESIREDGAINR_QP_PYTHONPATH",
  unset = ""
)

# Prevent an unrelated PYTHONHOME from changing this interpreter's prefix.
Sys.unsetenv("PYTHONHOME")

if (nzchar(pythonpath)) {
  Sys.setenv(PYTHONPATH = pythonpath)
}

message("QP oracle Python: ", python)

link_target <- Sys.readlink(python)

if (nzchar(link_target)) {
  message("QP oracle Python link target: ", link_target)
}

required_modules <- c(
  "clarabel",
  "osqp",
  "numpy",
  "scipy"
)

probe_code <- paste(
  "import sys",
  paste(
    sprintf("import %s", required_modules),
    collapse = "; "
  ),
  "print(sys.executable)",
  sep = "; "
)

probe <- system2(
  command = python,
  args = c(
    "-c",
    shQuote(probe_code)
  ),
  stdout = TRUE,
  stderr = TRUE
)

probe_status <- attr(probe, "status")

if (is.null(probe_status)) {
  probe_status <- 0L
}

if (!identical(as.integer(probe_status), 0L)) {
  cat(probe, sep = "\n")
  
  stop(
    paste0(
      "The configured Python interpreter cannot import the ",
      "independent QP solver dependencies. Interpreter: ",
      python
    ),
    call. = FALSE
  )
}

message(
  "QP oracle dependencies confirmed using: ",
  probe[[length(probe)]]
)

status <- system2(
  command = python,
  args = c(
    shQuote(oracle),
    shQuote(input),
    shQuote(output)
  )
)

if (!identical(as.integer(status), 0L)) {
  stop(
    "The Python solver oracle failed with exit status ",
    status,
    ". Interpreter: ",
    python,
    call. = FALSE
  )
}

if (!file.exists(output)) {
  stop(
    "The Python solver oracle completed without producing its output file.",
    call. = FALSE
  )
}



oracle_payload <- jsonlite::read_json(output, simplifyVector = FALSE)
versions <- oracle_payload$versions
external <- oracle_payload$results
names(external) <- vapply(external, `[[`, character(1L), "problem_id")

rows <- lapply(problems, function(problem) {
  H <- as.matrix(problem$H)
  f <- unlist(problem$f)
  lower <- vapply(problem$lower, function(x) {
    if (x == "Inf") Inf else if (x == "-Inf") -Inf else as.numeric(x)
  }, numeric(1L))
  upper <- vapply(problem$upper, function(x) {
    if (x == "Inf") Inf else if (x == "-Inf") -Inf else as.numeric(x)
  }, numeric(1L))
  internal <- DesiredGainR:::.dgr_box_qp(H, f, lower, upper)
  oracle_result <- external[[problem$problem_id]]
  summarise <- function(result, solver) {
    x <- unlist(result$solution)
    dual <- unlist(result$dual)
    lower_dual <- upper_dual <- numeric(length(x))
    if (grepl("clarabel", solver, fixed = TRUE)) {
      if (identical(solver, "clarabel_reordered")) dual <- rev(dual)
      cursor <- 1L
      for (j in seq_along(x)) {
        if (is.finite(upper[j])) {
          upper_dual[j] <- dual[cursor]
          cursor <- cursor + 1L
        }
        if (is.finite(lower[j])) {
          lower_dual[j] <- dual[cursor]
          cursor <- cursor + 1L
        }
      }
    } else {
      lower_dual <- pmax(-dual, 0)
      upper_dual <- pmax(dual, 0)
    }
    list(
      x = x,
      objective = as.numeric(result$objective),
      kkt = DesiredGainR:::.dgr_box_kkt(x, H, f, lower, upper)$overall,
      status = as.character(result$status),
      lower_dual = lower_dual, upper_dual = upper_dual
    )
  }
  clarabel <- summarise(oracle_result$clarabel, "clarabel")
  clarabel_reordered <- summarise(
    oracle_result$clarabel_reordered, "clarabel_reordered"
  )
  osqp_unpolished <- summarise(oracle_result$osqp_unpolished, "osqp")
  osqp_polished <- summarise(oracle_result$osqp_polished, "osqp")
  condition_number <- problem$condition_number
  x_tolerance <- if (condition_number <= 1e4) {
    1e-7
  } else if (condition_number <= 1e6) {
    5e-5
  } else if (condition_number <= 1e8) {
    2e-2
  } else {
    1
  }
  if (problem$class == "degenerate_active_set") x_tolerance <- 1e-5
  objective_tolerance <- if (condition_number <= 1e8) 1e-9 else 1e-6
  x_scale <- max(1, abs(internal$solution), abs(clarabel$x))
  objective_scale <- max(1, abs(internal$objective), abs(clarabel$objective))
  solution_difference <- max(abs(internal$solution - clarabel$x))
  osqp_unpolished_solution_difference <- max(abs(
    internal$solution - osqp_unpolished$x
  ))
  osqp_polished_solution_difference <- max(abs(
    internal$solution - osqp_polished$x
  ))
  objective_difference <- abs(internal$objective - clarabel$objective)
  osqp_unpolished_objective_difference <- abs(
    internal$objective - osqp_unpolished$objective
  )
  osqp_polished_objective_difference <- abs(
    internal$objective - osqp_polished$objective
  )
  analytical_difference <- if (is.null(problem$expected)) {
    NA_real_
  } else {
    max(abs(internal$solution - problem$expected))
  }
  require_solution <- condition_number <= 1e8
  clarabel_kkt_tolerance <- if (
    problem$class == "degenerate_active_set"
  ) {
    1e-5
  } else {
    1e-7
  }
  osqp_kkt_tolerance <- if (condition_number >= 1e6) 2e-6 else 1e-6
  active_tolerance <- max(1e-8, x_tolerance)
  active_lower <- sum(is.finite(lower) &
    abs(internal$solution - lower) <= active_tolerance)
  active_upper <- sum(is.finite(upper) &
    abs(internal$solution - upper) <= active_tolerance)
  fixed_bounds <- sum(is.finite(lower) & is.finite(upper) & lower == upper)
  fixed <- is.finite(lower) & is.finite(upper) & lower == upper
  internal_gradient <- as.numeric(H %*% internal$solution + f)
  internal_lower_dual <- ifelse(
    is.finite(lower) & abs(internal$solution - lower) <= active_tolerance,
    pmax(internal_gradient, 0), 0
  )
  internal_upper_dual <- ifelse(
    is.finite(upper) & abs(internal$solution - upper) <= active_tolerance,
    pmax(-internal_gradient, 0), 0
  )
  dual_difference <- function(result) {
    if (all(fixed)) {
      return(NA_real_)
    }
    max(abs(c(
      internal_lower_dual[!fixed] - result$lower_dual[!fixed],
      internal_upper_dual[!fixed] - result$upper_dual[!fixed]
    )))
  }
  pass <- grepl("Solved", clarabel$status, ignore.case = TRUE) &&
    grepl("Solved", clarabel_reordered$status, ignore.case = TRUE) &&
    grepl("solved", osqp_unpolished$status, ignore.case = TRUE) &&
    grepl("solved", osqp_polished$status, ignore.case = TRUE) &&
    internal$kkt$overall <= 1e-7 && clarabel$kkt <= clarabel_kkt_tolerance &&
    osqp_polished$kkt <= osqp_kkt_tolerance &&
    objective_difference <= objective_tolerance * objective_scale &&
    osqp_unpolished_objective_difference <=
      objective_tolerance * objective_scale &&
    osqp_polished_objective_difference <=
      objective_tolerance * objective_scale &&
    (!require_solution || solution_difference <= x_tolerance * max(x_scale)) &&
    (!require_solution || osqp_polished_solution_difference <=
      x_tolerance * max(x_scale)) &&
    max(abs(clarabel$x - clarabel_reordered$x)) <=
      x_tolerance * max(x_scale) &&
    (is.na(analytical_difference) || analytical_difference <= x_tolerance)
  data.frame(
    problem_id = problem$problem_id, class = problem$class,
    desiredgainr_version = as.character(utils::packageVersion("DesiredGainR")),
    clarabel_version = versions$clarabel, osqp_version = versions$osqp,
    numpy_version = versions$numpy, scipy_version = versions$scipy,
    random_seed = problem$seed, dimension = length(f),
    condition_number = condition_number,
    solution_tolerance = x_tolerance,
    objective_tolerance = objective_tolerance,
    active_lower = active_lower, active_upper = active_upper,
    fixed_bounds = fixed_bounds,
    internal_status = internal$status, clarabel_status = clarabel$status,
    clarabel_reordered_status = clarabel_reordered$status,
    osqp_unpolished_status = osqp_unpolished$status,
    osqp_polished_status = osqp_polished$status,
    internal_objective = internal$objective,
    clarabel_objective = clarabel$objective,
    osqp_unpolished_objective = osqp_unpolished$objective,
    osqp_polished_objective = osqp_polished$objective,
    maximum_solution_difference = solution_difference,
    relative_solution_difference = solution_difference / max(x_scale),
    clarabel_reordered_solution_difference = max(abs(
      clarabel$x - clarabel_reordered$x
    )),
    osqp_unpolished_solution_difference =
      osqp_unpolished_solution_difference,
    osqp_polished_solution_difference = osqp_polished_solution_difference,
    relative_objective_difference = objective_difference / objective_scale,
    osqp_unpolished_relative_objective_difference =
      osqp_unpolished_objective_difference / objective_scale,
    osqp_polished_relative_objective_difference =
      osqp_polished_objective_difference / objective_scale,
    internal_kkt_residual = internal$kkt$overall,
    internal_primal_residual = internal$kkt$primal,
    clarabel_kkt_residual = clarabel$kkt,
    clarabel_primal_residual = DesiredGainR:::.dgr_box_kkt(
      clarabel$x, H, f, lower, upper
    )$primal,
    osqp_unpolished_kkt_residual = osqp_unpolished$kkt,
    osqp_polished_kkt_residual = osqp_polished$kkt,
    osqp_polished_primal_residual = DesiredGainR:::.dgr_box_kkt(
      osqp_polished$x, H, f, lower, upper
    )$primal,
    clarabel_dual_difference = dual_difference(clarabel),
    osqp_unpolished_dual_difference = dual_difference(osqp_unpolished),
    osqp_polished_dual_difference = dual_difference(osqp_polished),
    analytical_difference = analytical_difference,
    pass = pass
  )
})
report <- do.call(rbind, rows)
report_path <- Sys.getenv(
  "DESIREDGAINR_QP_REPORT",
  unset = file.path("artifacts", "qp-validation.csv")
)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(report, report_path, row.names = FALSE)
print(table(report$class, report$pass))
cat("Report:", normalizePath(report_path, mustWork = FALSE), "\n")
if (!all(report$pass)) {
  print(report[!report$pass, ])
  stop(sum(!report$pass), " independent solver validation problem(s) failed.",
    call. = FALSE
  )
}
