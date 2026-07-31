# Single registry of data.table non-standard-evaluation symbols.
#
# These are column names referenced inside data.table expressions and `..`
# selections. Declaring them here suppresses spurious "no visible binding for
# global variable" NOTEs from R CMD check. This is the only place in the
# package that calls globalVariables(); keep it in sync when a new column is
# created by reference.

utils::globalVariables(c(
  # data.table `..` column selection
  "..trait_cols",

  # run_dgsi()
  "Chosen",
  "Eligible",
  "Replicate",
  "SelectionIndex",
  "Selected",

  # run_qgsi()
  "LinearPart",
  "QGSI",
  "QuadraticPart",
  "Rank",

  # run_qgsi_desired_gain() compatibility aliases
  "LinearDGPart",
  "QGSI_DG",
  "QuadraticDGPart",
  "Rank_QGSI_DG",

  # compare_dg_and_qgsi()
  "AbsRankDiff_DG_vs_QGSI",
  "DG_Rank",
  "DG_SelectionIndex",
  "QGSI_Rank",
  "RankDiff_DG_minus_QGSI",

  # selection_index()
  "id",
  "score",
  "selected",

  # simulate_selection_cycles()
  "cycle",
  "trait",
  "genetic_variance",

  # dosage_diagnostics()
  "variant",
  "heterozygosity",

  # optimize_desired_gains()
  "pareto_optimal"
))
