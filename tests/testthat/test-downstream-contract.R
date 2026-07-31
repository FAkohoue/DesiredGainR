# Contract relied upon by reverse dependencies.
#
# hapblockr::build_selection_index() delegates method = "dgsi" and
# method = "qgsi" to this package via .hb_build_desiredgain_index(). These
# tests pin exactly the call surface it uses and exactly the result elements
# it reads back, so a future refactor of DesiredGainR cannot silently break
# it. Do not relax these without checking the reverse dependency first.

# Reproduces the argument set that hapblockr passes, including the
# data.frame (not data.table) inputs and the non-default id_col.
downstream_inputs <- function(n = 24L) {
  traits <- c("YLD", "DIS")
  values <- matrix(
    c(stats::rnorm(n), stats::rnorm(n)),
    ncol = 2L, dimnames = list(paste0("cand", seq_len(n)), traits)
  )
  list(
    traits = traits,
    values = values,
    init = data.frame(id = rownames(values), stringsAsFactors = FALSE),
    candidate = data.frame(
      id = rownames(values), values,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  )
}

test_that("run_dgsi() keeps the formals that dgsi_control protects", {
  # hapblockr refuses to let dgsi_control override these, so each must
  # remain a formal argument of run_dgsi() with the same name.
  protected <- c(
    "init_data", "cand_data", "trait_cols", "dg", "P", "G",
    "id_col", "lower_is_better", "n_select"
  )
  expect_true(all(protected %in% names(formals(run_dgsi))))
})

test_that("run_qgsi() keeps the formals the downstream QGSI path passes", {
  required <- c(
    "init_data", "gebv_data", "trait_cols", "linear_weights", "W",
    "id_col", "lower_is_better"
  )
  expect_true(all(required %in% names(formals(run_qgsi))))
  # The downstream call does not pass center_traits, so the default must
  # stay TRUE; otherwise the estimated Gamma warning would start firing.
  expect_true(isTRUE(eval(formals(run_qgsi)$center_traits)))
})

test_that("DGSI result exposes every element the downstream reads", {
  set.seed(404)
  inputs <- downstream_inputs()
  traits <- inputs$traits
  G <- matrix(
    c(0.60, -0.15, -0.15, 0.40), 2,
    dimnames = list(traits, traits)
  )
  P <- matrix(
    c(1.10, -0.20, -0.20, 0.90), 2,
    dimnames = list(traits, traits)
  )

  engine <- run_dgsi(
    init_data = inputs$init,
    cand_data = inputs$candidate,
    trait_cols = traits,
    dg = abs(c(YLD = 0.5, DIS = 0.3)),
    P = P,
    G = G,
    id_col = "id",
    lower_is_better = "DIS",
    n_select = 4L
  )

  # engine$coefficients[traits] and engine$realised_response[traits]: both
  # must stay named numeric vectors indexable by trait name.
  expect_type(engine$coefficients, "double")
  expect_named(engine$coefficients, traits)
  expect_true(all(is.finite(as.numeric(engine$coefficients[traits]))))
  expect_type(engine$realised_response, "double")
  expect_named(engine$realised_response, traits)
  expect_true(all(is.finite(as.numeric(engine$realised_response[traits]))))

  # as.data.frame(engine$ranked_geno) then rename SelectionIndex.
  ranked <- as.data.frame(engine$ranked_geno)
  expect_true(is.data.frame(ranked))
  expect_true("SelectionIndex" %in% names(ranked))
  expect_true("id" %in% names(ranked))
  expect_true(all(is.finite(ranked$SelectionIndex)))
  expect_setequal(ranked$id, rownames(inputs$values))

  # Quality gate: sum(engine$replicate_diagnostics$Chosen) == 1L
  expect_true("Chosen" %in% names(engine$replicate_diagnostics))
  expect_identical(sum(engine$replicate_diagnostics$Chosen), 1L)
  # Used as the uncertainty table, so it must stay rectangular.
  expect_gt(nrow(engine$replicate_diagnostics), 0L)
})

test_that("QGSI result exposes every element the downstream reads", {
  set.seed(405)
  inputs <- downstream_inputs()
  traits <- inputs$traits
  W <- matrix(
    c(0.08, 0.01, 0.01, -0.05), 2,
    dimnames = list(traits, traits)
  )

  # The downstream call passes no center_traits, no Gamma, no n_select.
  # It must not emit warnings under that exact argument set.
  seen <- character()
  withCallingHandlers(
    engine <- run_qgsi(
      init_data = inputs$init,
      gebv_data = inputs$candidate,
      trait_cols = traits,
      linear_weights = abs(c(YLD = 1.0, DIS = 0.4)),
      W = W,
      id_col = "id",
      lower_is_better = "DIS"
    ),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(seen, character())

  ranked <- as.data.frame(engine$ranked_geno)
  expect_true("QGSI" %in% names(ranked))
  expect_true("id" %in% names(ranked))
  expect_true(all(is.finite(ranked$QGSI)))
  expect_setequal(ranked$id, rownames(inputs$values))

  # Used as the uncertainty table.
  expect_true(is.data.frame(engine$component_summary))
  expect_true(all(
    c("Component", "Mean", "SD", "Min", "Max") %in%
      names(engine$component_summary)
  ))
})

test_that("0.3.1 corrections do not touch the fields the downstream reads", {
  # The two numerical corrections in 0.3.1 land in expected_gain_per_trait
  # and theoretical_parameters. The downstream QGSI path reads neither: it
  # sets expected_response to NA and never inspects the theory block. This
  # test documents that separation so it stays true.
  set.seed(406)
  inputs <- downstream_inputs()
  traits <- inputs$traits
  W <- matrix(
    c(0.08, 0.01, 0.01, -0.05), 2,
    dimnames = list(traits, traits)
  )
  engine <- run_qgsi(
    init_data = inputs$init,
    gebv_data = inputs$candidate,
    trait_cols = traits,
    linear_weights = abs(c(YLD = 1.0, DIS = 0.4)),
    W = W,
    id_col = "id",
    lower_is_better = "DIS"
  )
  downstream_fields <- c("ranked_geno", "component_summary")
  expect_true(all(downstream_fields %in% names(engine)))
  # Superseded quantities remain reachable for anyone who did depend on them.
  expect_true(
    "Expected_Genetic_Gain_LinearSD" %in%
      names(engine$expected_gain_per_trait)
  )
  expect_true(
    "variance_ratio_index_to_merit" %in%
      names(engine$theoretical_parameters)
  )
})
