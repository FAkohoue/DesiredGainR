validation_file <- function(name) {
  path <- system.file("validation", name, package = "DesiredGainR")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "validation", name)
  }
  path
}

test_that("empirical validation artifacts are complete and candid", {
  summary <- utils::read.csv(
    validation_file("empirical-validation-summary.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  expect_gte(nrow(summary), 20L)
  expect_equal(anyDuplicated(summary$validation_id), 0L)
  expect_true(all(c(
    "validation_id", "dataset", "evidence_tier", "endpoint", "estimate",
    "lower", "upper", "unit", "n", "status", "interpretation",
    "limitation"
  ) %in% names(summary)))
  expect_false(anyNA(summary$status))
  expect_true(all(nzchar(summary$limitation)))

  boundary <- summary[
    summary$validation_id == "prospective_vector_comparison_design",
  ]
  expect_equal(nrow(boundary), 1L)
  expect_identical(boundary$status, "future-evidence-opportunity")
  expect_true(is.na(boundary$estimate))

  # The report must retain mixed findings; a table containing only supportive
  # endpoints would be evidence filtering, not validation.
  expect_true(any(summary$status == "mixed"))
  expect_true(any(summary$status == "supportive"))
  expect_true(any(summary$status == "programme-concordant"))
})

test_that("real-cycle and temporal anchors remain numerically plausible", {
  summary <- utils::read.csv(
    validation_file("empirical-validation-summary.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  row <- function(id) summary[summary$validation_id == id, , drop = FALSE]

  cna <- row("cna6_yld_check_contrast")
  expect_equal(nrow(cna), 1L)
  expect_gt(cna$estimate, 150)
  expect_lt(cna$estimate, 320)
  expect_gt(cna$lower, 0)

  hir <- row("hir_hif_connected_trend")
  expect_equal(nrow(hir), 1L)
  expect_gt(hir$estimate, 0)
  expect_gt(hir$lower, 0)

  inia <- row("inia_e1_advancement_auc")
  expect_equal(nrow(inia), 1L)
  expect_gt(inia$estimate, 0.7)
  expect_gt(inia$lower, 0.5)

  increment <- row("inia_e1_auc_increment_over_yield")
  expect_equal(nrow(increment), 1L)
  expect_gt(increment$lower, 0)
  expect_lt(increment$estimate, 0.1)

  g2f <- row("g2f_index_transport_rho")
  expect_equal(nrow(g2f), 1L)
  expect_gt(g2f$lower, 0)
  expect_lte(g2f$n, 100L)
})

test_that("empirical audit distinguishes trial keys from plot rows", {
  audit <- utils::read.csv(
    validation_file("empirical-data-audit.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  expect_gte(nrow(audit), 8L)
  expect_equal(anyDuplicated(audit$dataset), 0L)
  expect_true(all(nzchar(audit$principal_limitation)))
  expect_true(all(c(
    "genotyped_candidates", "markers", "founder_simulation_compatibility"
  ) %in% names(audit)))

  inia <- audit[audit$dataset == "INIA historical rice breeding programme", ]
  expect_equal(inia$rows, 91636L)
  expect_equal(inia$environments, 996L)
  expect_equal(inia$genotyped_candidates, 965L)
  expect_equal(inia$markers, 61260L)
  expect_match(inia$founder_simulation_compatibility, "39 of 965")

  dataverse <- audit[grepl("Ceron-Rojas", audit$dataset, fixed = TRUE), ]
  expect_equal(dataverse$rows, 239L)
  expect_match(dataverse$principal_limitation, "not present")
})

test_that("empirical report states the causal boundary", {
  report <- paste(
    readLines(validation_file("empirical-validation-report.md"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(report, "validate the package's", fixed = TRUE)
  expect_match(report, "evidence extension", fixed = TRUE)
  expect_match(report, "selected-only", fixed = TRUE)
})

test_that("CIMMYT RCGS validation reproduces gain and genomic signal", {
  summary <- utils::read.csv(
    validation_file("empirical-validation-summary.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  row <- function(id) summary[summary$validation_id == id, , drop = FALSE]

  c0_c4 <- row("rcgs_gain_c0_c4_slope")
  expect_equal(nrow(c0_c4), 1L)
  expect_gt(c0_c4$estimate, 0.12)
  expect_lt(c0_c4$estimate, 0.23)
  expect_identical(c0_c4$status, "programme-concordant")

  genomic <- row("rcgs_raw_subsampled_gblup_cv_accuracy")
  expect_equal(nrow(genomic), 1L)
  expect_gt(genomic$estimate, 0.40)
  expect_lt(genomic$estimate, 0.60)

  negative <- row("rcgs_raw_gblup_permutation_control")
  expect_equal(nrow(negative), 1L)
  expect_lt(abs(negative$estimate), 0.10)

  boundary <- row("rcgs_actual_founder_simulation_calibration")
  expect_identical(boundary$status, "not-compatible")
  expect_true(is.na(boundary$estimate))
})

test_that("CIMMYT RCGS environment validation is nested and auditable", {
  result <- utils::read.csv(
    validation_file("cimmyt-rcgs-environment-validation.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  expect_equal(nrow(result), 4L)
  direction_columns <- c(
    "direction_GY", "direction_EarRot", "direction_EarAspect",
    "direction_Moisture"
  )
  expect_equal(rowSums(result[, direction_columns]), rep(1, 4), tolerance = 1e-12)
  expect_true(all(result$predicted_joint_success_probability > 0 &
    result$predicted_joint_success_probability < 1))
  expect_true(all(c(
    "chosen_GY", "chosen_EarRot", "chosen_EarAspect", "chosen_Moisture",
    "equal_GY", "yield_only_GY"
  ) %in% names(result)))

  integrity <- utils::read.csv(
    validation_file("cimmyt-rcgs-deposit-integrity.csv"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  duplicate <- integrity[integrity$item == "C3/C4 selected HapMap identity", ]
  expect_equal(nrow(duplicate), 1L)
  expect_identical(duplicate$status, "fail-do-not-use-as-C4")
  expect_true(any(integrity$item == "training pedigree-to-sample linkage" &
    integrity$status == "pass"))
  expect_true(any(integrity$item == "C3 gain-trial pedigree linkage" &
    integrity$status == "pass"))
})
