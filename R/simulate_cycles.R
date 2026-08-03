# Multi-cycle recurrent-selection simulation for evaluating a desired-gain
# vector.
#
# The purpose of this simulation is narrow and should stay narrow: it exists to
# compare desired-gain directions over several cycles, not to design a crossing
# plan. Parent, cross and mating decisions belong to a dedicated mating-design
# tool. Consequently the selection rule here is deliberately simple.

#' Pull the genotypes used to measure diversity
#'
#' Diversity is measured on a **genome-wide marker panel**, not on the
#' quantitative trait loci.
#'
#' Measuring it on the QTL confounds two different things. Selection is
#' supposed to change QTL frequencies: that is what genetic gain *is*. A
#' relationship matrix built from QTL therefore rises whenever selection
#' succeeds, so using it as a diversity axis makes the optimiser trade gain
#' against a quantity that is partly gain itself, and a direction is penalised
#' for working. Coancestry is a property of the genome as a whole, and the
#' markers that are not under selection are what carry it.
#'
#' The panel is fixed at [founder_population()] time and stored in the setup,
#' so every measurement uses the same loci.
#'
#' @param pop An AlphaSimR population.
#' @param SP AlphaSimR simulation parameters.
#' @param panel Optional marker locus identifiers held in the setup. When
#'   absent, all segregating sites are used, which is still genome-wide.
#' @noRd
.dgr_diversity_geno <- function(pop, SP, panel = NULL) {
  geno <- tryCatch(
    AlphaSimR::pullSnpGeno(pop, simParam = SP),
    error = function(e) NULL
  )
  if (is.null(geno) || !length(geno)) {
    # No SNP chip was declared on this SimParam, so fall back to every
    # segregating site. That is still genome-wide and still excludes nothing
    # under selection preferentially.
    geno <- AlphaSimR::pullSegSiteGeno(pop, simParam = SP)
  }
  if (!is.null(panel) && length(panel)) {
    # The panel is required exactly. Intersecting silently, and falling back to
    # the full genotype matrix when nothing matched, would change the basis on
    # which diversity is measured without saying so -- and the fallback would
    # reintroduce the quantitative trait loci the panel exists to exclude.
    available <- colnames(geno)
    if (is.null(available)) {
      stop("The genotype matrix carries no marker names, so the stored ",
        "diversity panel cannot be matched against it.",
        call. = FALSE
      )
    }
    missing_loci <- setdiff(panel, available)
    if (length(missing_loci)) {
      stop(
        length(missing_loci), " of ", length(panel), " markers in the stored ",
        "diversity panel are absent from this population, including ",
        paste(utils::head(missing_loci, 3L), collapse = ", "),
        ". Diversity would otherwise be measured on a different basis from ",
        "the founders, making the two incomparable.",
        call. = FALSE
      )
    }
    # Ordered by identifier so the column order cannot vary between calls.
    geno <- geno[, sort(panel), drop = FALSE]
  }
  geno
}

#' Allele frequencies of the founder generation on the marker panel
#'
#' Relatedness accumulates relative to a fixed base. Recomputing allele
#' frequencies from each generation would destroy exactly the signal being
#' measured, so the base is established once and reused.
#'
#' @param pop The founder population.
#' @param SP AlphaSimR simulation parameters.
#' @param panel Optional marker identifiers.
#'
#' @return A list holding the retained loci, their founder frequencies, and the
#'   VanRaden scaling constant.
#' @noRd
.dgr_base_frequency <- function(pop, SP, panel = NULL) {
  geno <- .dgr_diversity_geno(pop, SP, panel)
  frequency <- colMeans(geno) / 2
  polymorphic <- frequency > 0 & frequency < 1
  if (!any(polymorphic)) {
    stop("The founder population is monomorphic at every marker, so ",
      "relatedness cannot be measured.",
      call. = FALSE
    )
  }
  # The retained loci are identified by name rather than by position. A
  # logical vector is only meaningful against the exact column order it was
  # built from, so any later population with a different order would be
  # silently measured on the wrong loci.
  retained <- colnames(geno)[polymorphic]
  if (is.null(retained) || anyNA(retained)) {
    stop("Diversity markers must carry names so that the same loci can be ",
      "measured in every generation.",
      call. = FALSE
    )
  }
  frequency <- frequency[polymorphic]
  names(frequency) <- retained
  list(
    retained = retained,
    frequency = frequency,
    scaling = 2 * sum(frequency * (1 - frequency)),
    n_markers = length(frequency),
    panel = panel,
    basis = if (is.null(panel)) {
      "all segregating sites (no marker panel was declared)"
    } else {
      "genome-wide markers"
    }
  )
}

#' Founder frequencies at the quantitative trait loci
#'
#' Reported separately from marker-based coancestry, because the two answer
#' different questions: this one measures how much of the *functional*
#' variation has been used up, which is informative but is not diversity.
#'
#' @noRd
.dgr_base_qtl_frequency <- function(pop, SP) {
  geno <- AlphaSimR::pullQtlGeno(pop, simParam = SP)
  frequency <- colMeans(geno) / 2
  retained <- frequency > 0 & frequency < 1
  if (!any(retained)) {
    return(NULL)
  }
  list(retained = retained, frequency = frequency[retained])
}

#' Proportion of founder QTL variation still segregating
#'
#' @return The fraction of initially polymorphic QTL still segregating, and the
#'   genic variance relative to the founders.
#' @noRd
.dgr_qtl_diversity <- function(pop, SP, qtl_base) {
  if (is.null(qtl_base)) {
    return(list(segregating = NA_real_, genic_variance_ratio = NA_real_))
  }
  geno <- AlphaSimR::pullQtlGeno(pop, simParam = SP)
  frequency <- colMeans(geno[, qtl_base$retained, drop = FALSE]) / 2
  heterozygosity <- 2 * frequency * (1 - frequency)
  founder_heterozygosity <- 2 * qtl_base$frequency * (1 - qtl_base$frequency)
  list(
    segregating = mean(frequency > 0 & frequency < 1),
    genic_variance_ratio = sum(heterozygosity) / sum(founder_heterozygosity)
  )
}

#' Mean pairwise genomic relationship relative to the founder base
#'
#' Uses the VanRaden (2008) construction with the allele frequencies **of the
#' founder generation**, not of the population being measured. Centring on a
#' population's own frequencies forces every column of the centred genotype
#' matrix to sum to zero, so the whole relationship matrix sums to zero and the
#' mean off-diagonal element collapses to \eqn{-\overline{G}_{ii}/(n-1)}: a
#' function of the number of individuals and nothing else. Relatedness measured
#' against a moving base cannot accumulate, which defeats the purpose.
#'
#' @param pop An AlphaSimR population.
#' @param SP AlphaSimR simulation parameters.
#' @param base A list from `.dgr_base_frequency()`.
#'
#' @return A list with the mean off-diagonal relationship and the mean genomic
#'   inbreeding coefficient, both relative to the founders.
#' @noRd
.dgr_mean_relationship <- function(pop, SP, base) {
  geno <- .dgr_diversity_geno(pop, SP, base$panel)
  if (is.null(dim(geno)) || nrow(geno) < 2L) {
    return(list(relationship = NA_real_, inbreeding = NA_real_))
  }
  absent <- setdiff(base$retained, colnames(geno))
  if (length(absent)) {
    stop(length(absent), " of the ", length(base$retained), " founder ",
      "diversity markers are absent from this population, so relatedness ",
      "would be measured against a different base.",
      call. = FALSE
    )
  }
  geno <- geno[, base$retained, drop = FALSE]
  centred <- sweep(geno, 2L, 2 * base$frequency, "-")
  relationship <- tcrossprod(centred) / base$scaling
  list(
    relationship = mean(relationship[upper.tri(relationship)]),
    # The diagonal of a VanRaden matrix is 1 + F, so the mean genomic
    # inbreeding coefficient relative to the founders is the mean diagonal
    # less one.
    inbreeding = mean(diag(relationship)) - 1
  )
}

#' Effective population size implied by a change in inbreeding
#'
#' @param before,after Mean genomic inbreeding before and after one cycle.
#'
#' @return An approximate effective population size, or `NA_real_` when
#'   inbreeding did not increase, which happens through sampling noise at small
#'   population sizes and does not mean the effective size is infinite.
#' @noRd
.dgr_effective_size <- function(before, after) {
  if (!is.finite(before) || !is.finite(after)) {
    return(NA_real_)
  }
  delta <- (after - before) / (1 - before)
  if (!is.finite(delta) || delta <= 0) {
    return(NA_real_)
  }
  1 / (2 * delta)
}

#' Index coefficients for a desired-gain direction
#'
#' Uses the Yamada formulation, which is the one [run_dgsi()] implements.
#'
#' @param d Desired-gain direction.
#' @param G,P Covariance matrices in the analysis space.
#' @param ridge Small ridge added to both solves.
#'
#' @return Named index coefficients.
#' @noRd
.dgr_yamada_coefficients <- function(d, G, P, ridge = 1e-8) {
  p <- length(d)
  P_ridged <- P + diag(ridge, p)
  Xpg <- solve(P_ridged, G)
  middle <- crossprod(G, Xpg)
  middle <- (middle + t(middle)) / 2 + diag(ridge, p)
  coefficients <- as.numeric(Xpg %*% solve(middle, d))
  names(coefficients) <- names(d)
  coefficients
}

#' Simulate recurrent selection under a fixed desired-gain direction
#'
#' Runs a recurrent-selection programme forward for a stated number of cycles,
#' selecting each cycle on a desired-gain index built from the supplied
#' direction, and records the cumulative genetic gain and the loss of genetic
#' variation. The purpose is to compare desired-gain directions, because a
#' direction that maximises response in the first cycle may exhaust the
#' variance that response depends on within a few more.
#'
#' @details
#' # What the simulation adds beyond the single-cycle formula
#'
#' A single-cycle response prediction cannot distinguish desired-gain
#' directions beyond what the achievable-response ellipsoid already states, and
#' [gain_feasibility()] gives that answer exactly and without simulation.
#' Simulation is informative only because it represents the processes that make
#' cycle five differ from cycle one: (i) the reduction and reshaping of genetic
#' variance caused by truncation selection, (ii) drift and the accumulation of
#' relatedness in a finite population, (iii) the point at which an antagonistic
#' correlation drives a secondary trait past an unacceptable level, and (iv)
#' the compounding of estimation error when the index is rebuilt each cycle.
#'
#' # What a cycle means
#'
#' This is the definition to check before interpreting any number below.
#'
#' **Cycle 0** is the founder population, before any selection. It exists so
#' that the cycle 1 row has something to be a response *to*.
#'
#' **Cycle \eqn{t}** records the population produced by selecting parents from
#' the cycle \eqn{t-1} candidates and crossing them. Every row from cycle 1
#' onward is therefore a transmitted selection response, and `n_cycles = 1`
#' gives exactly one such response.
#'
#' Earlier versions measured the candidates rather than their selected
#' descendants. Because random mating does not shift a population mean, that
#' made the cycle 1 gain zero in expectation *for every desired-gain
#' direction*, so a one-cycle run could not distinguish directions at all.
#' Results from before this change should be re-run.
#'
#' Genetic mean, variance, relatedness and inbreeding are all measured on the
#' same population, the cycle's response. `parent_inbreeding` additionally
#' reports the selected parents, which is a different and smaller group.
#'
#' # Mating systems
#'
#' \describe{
#'   \item{`"self"`}{Self-pollinated line development. Selected parents are
#'     intercrossed, and the resulting families are advanced by selfing or
#'     doubled haploidy before evaluation. Recombination therefore releases
#'     variation slowly.}
#'   \item{`"outcross"`}{Recurrent selection in a random-mating population.
#'     Selected parents are intercrossed each cycle, so half the disequilibrium
#'     generated by selection decays per generation.}
#'   \item{`"clonal"`}{Clonally propagated crops, with the sexual and clonal
#'     phases kept distinct. Recombination happens once per cycle, when
#'     selected parents are crossed to raise a seedling generation; every
#'     evaluation stage after that is a copy of a seedling, so no further
#'     meiosis occurs and dominance is transmitted intact. Selection therefore
#'     acts on total genotypic value rather than breeding value, and response
#'     is a broad-sense quantity. Requires a setup built with
#'     `dominance_degree`, and `G` should be the genotypic covariance. See
#'     [founder_population()] for how the additive and dominance components are
#'     calibrated to that target.}
#' }
#'
#' # Index re-estimation
#'
#' With `reestimate_index = TRUE`, the covariance of the observed selection
#' criterion is re-estimated each cycle. Phenotypic selection retains the
#' breeder-supplied `G_target`; it never estimates genetic covariance from
#' hidden simulated breeding values. RR-BLUP uses the calibrated-prediction
#' identity `Cov(A, Ahat) = Var(Ahat)` and reports held-out prediction accuracy
#' as a diagnostic. With `FALSE`, cycle-one coefficients are reused. Thus no
#' selection decision has access to simulation truth.
#'
#' @param setup An object from [founder_population()].
#' @param desired_gains Named desired-gain direction. Only the direction
#'   matters, because the attainable magnitude is fixed by selection intensity.
#' @param n_cycles Number of selection cycles.
#' @param mating_system One of `"self"`, `"outcross"`, or `"clonal"`.
#' @param n_parents Number of parents recycled each cycle. Fewer parents
#'   increase short-term gain and reduce long-term potential.
#' @param n_crosses Number of crosses made each cycle.
#' @param n_progeny_per_cross Progeny evaluated per cross.
#' @param n_selfing_generations Selfing generations before evaluation, used by
#'   `"self"` when `use_doubled_haploids` is `FALSE`.
#' @param use_doubled_haploids Whether `"self"` produces doubled haploids
#'   instead of advancing by selfing.
#' @param lower_is_better Traits for which smaller values are favourable.
#' @param reestimate_index Whether to rebuild the index from each cycle's own
#'   simulated data. See Details.
#' @param n_clonal_replicates Ramets evaluated per genotype in a clonal
#'   programme. A clonal trial phenotypes several copies of each genotype and
#'   averages them, so its selection is more accurate than a single plot;
#'   leaving this at 1 understates the response a clonal programme achieves.
#'   Only meaningful with `mating_system = "clonal"`.
#' @param n_threads Threads AlphaSimR may use. The default of 1 is deliberate:
#'   above one, the order in which random numbers are consumed is not
#'   guaranteed, so a run cannot be reproduced exactly from its seed. Values
#'   above 1 warn.
#' @param seed Random seed. The caller's random number generator state is
#'   restored on exit, and `setup$SP` is deep-cloned so that the caller's
#'   `SimParam` is not advanced by the simulation.
#' @param prediction A named list describing the selection criterion. The
#'   default `list(method = "phenotype")` preserves phenotypic selection.
#'   `method = "rrblup"` uses leakage-free cross-fitted RR-BLUP predictions and
#'   requires a marker panel; optional fields are `folds` (5), `max_iter`
#'   (10000), `update_training` (`TRUE`) and `max_training` (`NULL`).
#'
#' @return An object of class `desiredgainr_simulation` containing a per-cycle
#'   table of genetic means, genetic variances, mean relationship, effective
#'   population size, prediction accuracy and prediction-calibration slope,
#'   together with the cumulative gain per trait.
#'
#' @seealso [founder_population()], [gain_feasibility()]
#' @export
simulate_selection_cycles <- function(
  setup,
  desired_gains,
  n_cycles = 5L,
  mating_system = c("self", "outcross", "clonal"),
  n_parents = 20L,
  n_crosses = 50L,
  n_progeny_per_cross = 10L,
  n_selfing_generations = 3L,
  use_doubled_haploids = FALSE,
  lower_is_better = NULL,
  reestimate_index = TRUE,
  n_clonal_replicates = 1L,
  n_threads = 1L,
  seed = 42L,
  prediction = list(method = "phenotype")
) {
  .dgr_require_alphasimr("simulate_selection_cycles()")
  mating_system <- match.arg(mating_system)
  if (!inherits(setup, "desiredgainr_sim_setup")) {
    stop("setup must be created by founder_population().", call. = FALSE)
  }
  if (mating_system == "clonal" && !isTRUE(setup$dominance)) {
    stop(
      "mating_system = 'clonal' requires a setup built with ",
      "dominance_degree, because selection acts on total genetic value.",
      call. = FALSE
    )
  }
  prediction <- .dgr_merge_named_list(
    list(
      method = "phenotype", folds = 5L, update_training = TRUE,
      max_training = NULL, max_iter = 10000L
    ),
    prediction, "prediction"
  )
  prediction$method <- match.arg(
    prediction$method, c("phenotype", "rrblup")
  )
  prediction$folds <- .dgr_positive_integer(
    prediction$folds, "prediction$folds"
  )
  prediction$max_iter <- .dgr_positive_integer(
    prediction$max_iter, "prediction$max_iter"
  )
  if (!is.logical(prediction$update_training) ||
    length(prediction$update_training) != 1L ||
    is.na(prediction$update_training)) {
    stop("prediction$update_training must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(prediction$max_training)) {
    prediction$max_training <- .dgr_positive_integer(
      prediction$max_training, "prediction$max_training"
    )
  }
  if (identical(prediction$method, "rrblup")) {
    if (is.null(setup$marker_panel) || !length(setup$marker_panel)) {
      stop("prediction$method = 'rrblup' requires a setup with a marker panel; ",
        "supply n_markers_per_chromosome to founder_population().",
        call. = FALSE
      )
    }
    if (identical(mating_system, "clonal") || isTRUE(setup$dominance)) {
      stop("The RR-BLUP pathway currently supports additive, non-clonal ",
        "programmes only; it predicts breeding value rather than total ",
        "genotypic value.",
        call. = FALSE
      )
    }
  }
  trait_cols <- setup$trait_cols
  d <- .dgr_named_vector(desired_gains, trait_cols, "desired_gains")
  if (all(d == 0)) {
    stop("desired_gains must contain at least one non-zero element.",
      call. = FALSE
    )
  }
  n_cycles <- .dgr_positive_integer(n_cycles, "n_cycles")
  n_parents <- .dgr_positive_integer(n_parents, "n_parents")
  n_crosses <- .dgr_positive_integer(n_crosses, "n_crosses")
  n_progeny_per_cross <- .dgr_positive_integer(
    n_progeny_per_cross, "n_progeny_per_cross"
  )
  # Zero selfing generations is meaningful: it evaluates the F1.
  n_selfing_generations <- .dgr_non_negative_integer(
    n_selfing_generations, "n_selfing_generations"
  )
  n_threads <- .dgr_positive_integer(n_threads, "n_threads")
  n_clonal_replicates <- .dgr_positive_integer(
    n_clonal_replicates, "n_clonal_replicates"
  )
  if (n_clonal_replicates > 1L && !identical(mating_system, "clonal")) {
    stop("n_clonal_replicates > 1 describes replicated clonal evaluation and ",
      "is only meaningful with mating_system = 'clonal'.",
      call. = FALSE
    )
  }
  seed <- .dgr_seed(seed)
  if (n_parents > n_crosses * n_progeny_per_cross) {
    stop("n_parents (", n_parents, ") exceeds the number of individuals ",
      "produced each cycle (", n_crosses * n_progeny_per_cross,
      "). Selection cannot retain more parents than exist.",
      call. = FALSE
    )
  }

  direction <- rep(1, length(trait_cols))
  names(direction) <- trait_cols
  if (length(lower_is_better)) {
    unknown <- setdiff(lower_is_better, trait_cols)
    if (length(unknown)) {
      stop("Unknown lower_is_better traits: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    direction[lower_is_better] <- -1
  }

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)

  # AlphaSimR's SimParam is an R6 object, so every crossing call mutates it by
  # reference: lastId advances and the caller's setup is silently changed.
  # Running the same setup twice would then not reproduce, and an optimiser
  # that evaluates many directions from one setup would have each evaluation
  # depend on the ones before it. Work on a deep clone instead.
  SP <- setup$SP$clone(deep = TRUE)
  caller_last_id <- setup$SP$lastId
  thread_note <- NULL
  if (n_threads > 1L) {
    thread_note <- paste(
      "n_threads =", n_threads, "was requested. AlphaSimR parallelises some",
      "operations across threads, and the order in which random numbers are",
      "consumed is then not guaranteed, so results may not reproduce exactly",
      "from the same seed. Use n_threads = 1 for an auditable result."
    )
    warning(thread_note, call. = FALSE)
  }
  # Record the thread count that took effect, not the one requested. If the
  # assignment fails on an older AlphaSimR the run can remain multithreaded
  # while the provenance claims a single thread, which is worse than not
  # recording it at all: the result would carry a false reproducibility claim.
  assigned <- tryCatch(
    {
      SP$nThreads <- n_threads
      TRUE
    },
    error = function(e) FALSE
  )
  effective_threads <- tryCatch(as.integer(SP$nThreads),
    error = function(e) NA_integer_
  )
  if (!isTRUE(assigned) || !identical(effective_threads, n_threads)) {
    thread_note <- paste(
      "The thread count could not be set to", n_threads,
      "on AlphaSimR", as.character(utils::packageVersion("AlphaSimR")),
      "- the effective count is",
      if (is.na(effective_threads)) "unknown" else effective_threads,
      ". This run's reproducibility from its seed is not guaranteed."
    )
    warning(thread_note, call. = FALSE)
  }

  # Genetic gain must be tracked on an absolute scale. AlphaSimR's bv() returns
  # breeding values expressed relative to the supplied population, so its
  # column means are zero for every population and a difference between cycles
  # would always be zero. gv() is anchored by the trait mean fixed when the
  # trait was added, and is therefore the quantity whose change across cycles
  # is genetic gain.
  merit_of <- function(pop) {
    values <- AlphaSimR::gv(pop)
    colnames(values) <- trait_cols
    values
  }

  # Estimating the covariance used to rebuild the index is a separate question.
  # For an additive programme the index should be built on additive covariance,
  # and centring does not affect a covariance, so bv() is appropriate there.
  # For a clonal programme the clone inherits dominance intact, so the
  # genotypic covariance is the relevant one.
  genetic_for_covariance <- function(pop) {
    values <- if (mating_system == "clonal") {
      AlphaSimR::gv(pop)
    } else {
      AlphaSimR::bv(pop, simParam = SP)
    }
    colnames(values) <- trait_cols
    values
  }

  # The base against which relatedness is measured is fixed at the founders and
  # never recomputed, and it is measured on genome-wide markers rather than on
  # the quantitative trait loci. See .dgr_diversity_geno() for why both matter.
  base <- .dgr_base_frequency(setup$founder_pop, SP, setup$marker_panel)
  qtl_base <- .dgr_base_qtl_frequency(setup$founder_pop, SP)

  # Develop the response population from a set of selected parents. This is the
  # only place recombination happens, and for the clonal system it is the
  # explicit sexual phase: clones do not recombine, so a clonal programme
  # crosses its selected parents to raise a new seedling generation and then
  # propagates those seedlings without further meiosis.
  develop_progeny <- function(parents) {
    seedlings <- AlphaSimR::randCross(
      parents,
      nCrosses = n_crosses, nProgeny = n_progeny_per_cross,
      simParam = SP
    )
    if (mating_system == "self") {
      if (isTRUE(use_doubled_haploids)) {
        return(AlphaSimR::makeDH(seedlings, nDH = 1L, simParam = SP))
      }
      advanced <- seedlings
      for (generation in seq_len(n_selfing_generations)) {
        advanced <- AlphaSimR::self(advanced, nProgeny = 1L, simParam = SP)
      }
      return(advanced)
    }
    # Outcrossing evaluates the F1 directly. The clonal system also evaluates
    # the seedling generation, but every subsequent stage is a copy of it, so
    # no further recombination occurs and dominance is transmitted intact.
    seedlings
  }

  measure <- function(pop) {
    values <- merit_of(pop)
    diversity <- .dgr_mean_relationship(pop, SP, base)
    qtl <- .dgr_qtl_diversity(pop, SP, qtl_base)
    list(
      mean = colMeans(values),
      variance = apply(values, 2L, stats::var),
      relationship = diversity$relationship,
      inbreeding = diversity$inbreeding,
      qtl_segregating = qtl$segregating,
      genic_variance_ratio = qtl$genic_variance_ratio,
      n = AlphaSimR::nInd(pop)
    )
  }

  founders <- measure(setup$founder_pop)
  baseline_merit <- founders$mean
  names(baseline_merit) <- trait_cols

  # Cycle 0 is the founder population itself: the state before any selection
  # has been applied. Recording it separately is what makes the cycle 1 row a
  # transmitted response rather than a baseline.
  records <- vector("list", n_cycles + 1L)
  records[[1L]] <- data.table::data.table(
    cycle = 0L,
    trait = trait_cols,
    genetic_mean = as.numeric(founders$mean),
    cumulative_gain = 0,
    genetic_variance = as.numeric(founders$variance),
    mean_relationship = founders$relationship,
    mean_inbreeding = founders$inbreeding,
    parent_inbreeding = NA_real_,
    qtl_segregating = founders$qtl_segregating,
    genic_variance_ratio = founders$genic_variance_ratio,
    effective_size = NA_real_,
    selection_intensity = NA_real_,
    n_evaluated = founders$n,
    n_parents_selected = NA_integer_,
    prediction_accuracy = NA_real_,
    prediction_calibration_slope = NA_real_
  )

  current <- setup$founder_pop
  previous_inbreeding <- founders$inbreeding
  fixed_coefficients <- NULL
  training_history <- NULL

  for (cycle in seq_len(n_cycles)) {
    # Selection in cycle t acts on the candidates available at the start of
    # cycle t, and the population it produces is the response for cycle t.
    # Measuring the candidates instead, as this loop previously did, records
    # an unselected random cross whose mean equals the previous generation's
    # in expectation, so cycle 1 showed no gain for any desired-gain direction
    # and n_cycles = 1 could not distinguish one direction from another.
    # A clonal programme evaluates each genotype as several ramets, so its
    # phenotypes average over replicates and are more precise than a single
    # plot. Ignoring that overstates the error variance and understates the
    # response, because selection is more accurate than the simulation shows.
    candidates <- if (n_clonal_replicates > 1L) {
      AlphaSimR::setPheno(
        current,
        reps = n_clonal_replicates, simParam = SP
      )
    } else {
      AlphaSimR::setPheno(current, simParam = SP)
    }
    n_available <- AlphaSimR::nInd(candidates)
    if (n_parents > n_available) {
      stop("n_parents (", n_parents, ") exceeds the ", n_available,
        " candidates available in cycle ", cycle, ".",
        call. = FALSE
      )
    }

    phenotypes <- AlphaSimR::pheno(candidates)
    colnames(phenotypes) <- trait_cols
    oriented_pheno <- sweep(phenotypes, 2L, direction, "*")
    oriented_genetic <- sweep(
      genetic_for_covariance(candidates), 2L, direction, "*"
    )

    selection_values <- oriented_pheno
    if (identical(prediction$method, "rrblup")) {
      n_folds <- min(prediction$folds, n_available - 1L)
      if (n_folds < 2L || n_available < 4L) {
        stop("RR-BLUP cross-fitting requires at least four candidates per cycle.",
          call. = FALSE
        )
      }
      fold <- integer(n_available)
      fold[order(as.character(candidates@id))] <- rep(
        seq_len(n_folds),
        length.out = n_available
      )
      predicted <- matrix(
        NA_real_,
        nrow = n_available, ncol = length(trait_cols),
        dimnames = list(NULL, trait_cols)
      )
      for (held_out in seq_len(n_folds)) {
        test <- which(fold == held_out)
        training <- candidates[-test]
        if (isTRUE(prediction$update_training) && !is.null(training_history)) {
          training <- AlphaSimR::mergePops(list(training_history, training))
        }
        solution <- AlphaSimR::RRBLUP(
          training,
          traits = seq_along(trait_cols), use = "pheno",
          snpChip = 1L, maxIter = prediction$max_iter, simParam = SP
        )
        held_out_pop <- AlphaSimR::setEBV(
          candidates[test], solution,
          value = "bv", simParam = SP
        )
        predicted[test, ] <- AlphaSimR::ebv(held_out_pop)
      }
      selection_values <- sweep(predicted, 2L, direction, "*")
      if (isTRUE(prediction$update_training)) {
        training_history <- if (is.null(training_history)) {
          candidates
        } else {
          AlphaSimR::mergePops(list(training_history, candidates))
        }
        if (!is.null(prediction$max_training) &&
          AlphaSimR::nInd(training_history) > prediction$max_training) {
          retained <- utils::tail(
            seq_len(AlphaSimR::nInd(training_history)),
            prediction$max_training
          )
          training_history <- training_history[retained]
        }
      }
    }
    prediction_accuracy <- vapply(seq_along(trait_cols), function(trait) {
      if (stats::sd(selection_values[, trait]) == 0 ||
        stats::sd(oriented_genetic[, trait]) == 0) {
        return(NA_real_)
      }
      stats::cor(selection_values[, trait], oriented_genetic[, trait])
    }, numeric(1L))
    prediction_calibration_slope <- vapply(
      seq_along(trait_cols),
      function(trait) {
        variance <- stats::var(selection_values[, trait])
        if (!is.finite(variance) || variance <= 0) {
          return(NA_real_)
        }
        stats::cov(
          oriented_genetic[, trait], selection_values[, trait]
        ) / variance
      },
      numeric(1L)
    )

    coefficients <- if (!is.null(fixed_coefficients) &&
      !isTRUE(reestimate_index)) {
      fixed_coefficients
    } else {
      P_hat <- stats::cov(selection_values)
      if (identical(prediction$method, "rrblup")) {
        # With calibrated GEBV, Cov(A, A-hat) = Var(A-hat). Cross-fitting keeps
        # this working assumption from being inflated by training-set reuse;
        # the realised accuracies below diagnose it using simulation truth.
        G_hat <- P_hat
      } else {
        signs <- diag(direction, nrow = length(direction))
        G_hat <- signs %*% setup$G_target %*% signs
      }
      dimnames(P_hat) <- dimnames(G_hat) <- list(trait_cols, trait_cols)
      .dgr_yamada_coefficients(d, G_hat, P_hat)
    }
    if (is.null(fixed_coefficients)) fixed_coefficients <- coefficients

    score <- as.numeric(selection_values %*% coefficients)
    # Ties are broken by AlphaSimR identifier rather than by row order, so the
    # selected set does not depend on the order individuals happen to arrive
    # in. The identifier is the population's `id` slot; AlphaSimR exports no
    # accessor for it.
    keep <- order(-score, as.character(candidates@id))[seq_len(n_parents)]
    parents <- candidates[keep]
    parent_diversity <- .dgr_mean_relationship(parents, SP, base)

    response <- develop_progeny(parents)
    observed <- measure(response)
    names(observed$mean) <- trait_cols

    records[[cycle + 1L]] <- data.table::data.table(
      cycle = cycle,
      trait = trait_cols,
      genetic_mean = as.numeric(observed$mean),
      cumulative_gain = as.numeric(observed$mean - baseline_merit),
      genetic_variance = as.numeric(observed$variance),
      mean_relationship = observed$relationship,
      mean_inbreeding = observed$inbreeding,
      parent_inbreeding = parent_diversity$inbreeding,
      # Reported beside coancestry but never substituted for it: QTL frequency
      # change is what selection is supposed to cause.
      qtl_segregating = observed$qtl_segregating,
      genic_variance_ratio = observed$genic_variance_ratio,
      effective_size = .dgr_effective_size(
        previous_inbreeding, observed$inbreeding
      ),
      selection_intensity = .dgr_intensity(n_parents / n_available),
      n_evaluated = observed$n,
      n_parents_selected = n_parents,
      prediction_accuracy = prediction_accuracy,
      prediction_calibration_slope = prediction_calibration_slope
    )
    previous_inbreeding <- observed$inbreeding
    current <- response
  }

  cycle_table <- data.table::rbindlist(records)
  final <- cycle_table[cycle == max(cycle)]
  cumulative <- stats::setNames(final$cumulative_gain, final$trait)

  result <- list(
    desired_gains = d,
    direction = direction,
    mating_system = mating_system,
    n_cycles = n_cycles,
    n_parents = n_parents,
    reestimate_index = reestimate_index,
    prediction = prediction,
    cycle_table = cycle_table,
    cumulative_gain = cumulative,
    final_variance = stats::setNames(final$genetic_variance, final$trait),
    final_relationship = final$mean_relationship[1L],
    final_inbreeding = final$mean_inbreeding[1L],
    diversity_basis = paste(
      "Relatedness and inbreeding are measured on", base$n_markers,
      "genome-wide markers against the founder allele frequencies, which are",
      "fixed at the first generation. Two choices matter here. Recomputing the",
      "frequencies each cycle would pin the mean off-diagonal at",
      "-mean(diagonal)/(n-1) and prevent any accumulation from being detected.",
      "Measuring on quantitative trait loci instead of markers would make the",
      "metric rise whenever selection succeeded, since changing QTL",
      "frequencies is what genetic gain is, so a direction would be penalised",
      "for working. QTL frequency change is reported separately as",
      "qtl_segregating and genic_variance_ratio."
    ),
    coefficients_cycle_one = fixed_coefficients,
    seed = seed,
    # Everything needed to reproduce this run, and to detect that a stored
    # result was produced under different conditions from the current session.
    provenance = list(
      seed = seed,
      n_threads_requested = n_threads,
      n_threads_effective = effective_threads,
      # Retained under the original name so that anything reading it keeps
      # working, but it now holds the count that actually applied.
      n_threads = effective_threads,
      thread_warning = thread_note,
      marker_basis = base$basis,
      n_diversity_markers = base$n_markers,
      alphasimr_version = as.character(utils::packageVersion("AlphaSimR")),
      desiredgainr_version = as.character(
        utils::packageVersion("DesiredGainR")
      ),
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      n_crosses = n_crosses,
      n_progeny_per_cross = n_progeny_per_cross,
      n_selfing_generations = n_selfing_generations,
      n_clonal_replicates = n_clonal_replicates,
      prediction = prediction,
      use_doubled_haploids = isTRUE(use_doubled_haploids),
      n_founders = AlphaSimR::nInd(setup$founder_pop),
      n_qtl = ncol(AlphaSimR::pullQtlGeno(setup$founder_pop, simParam = SP)),
      simparam_isolated = TRUE,
      caller_last_id_unchanged = identical(caller_last_id, setup$SP$lastId)
    ),
    cycle_semantics = paste(
      "Cycle 0 is the founder population before any selection. Cycle t records",
      "the population produced by selecting parents from the cycle t-1",
      "candidates and crossing them, so every row from cycle 1 onward is a",
      "transmitted selection response rather than an unselected sample."
    ),
    scope = paste(
      "This simulation compares desired-gain directions. It is not a crossing",
      "plan, and its selection rule is deliberately simple truncation on the",
      "index. Parent, cross and mating allocation require a dedicated",
      "mating-design tool."
    )
  )
  class(result) <- c("desiredgainr_simulation", "list")
  result
}

#' @export
print.desiredgainr_simulation <- function(x, ...) {
  cat("<desiredgainr_simulation>\n")
  cat(sprintf(
    "  %s system, %d cycles, %d parents recycled\n",
    x$mating_system, x$n_cycles, x$n_parents
  ))
  cat(
    "  Index re-estimated each cycle:",
    if (isTRUE(x$reestimate_index)) "yes" else "no", "\n"
  )
  cat("  Selection criterion:", x$prediction$method, "\n")
  cat("  Cumulative genetic gain:\n")
  print(round(x$cumulative_gain, 4L))
  cat(sprintf(
    "  Final relatedness: %.4f   inbreeding: %.4f\n",
    x$final_relationship, x$final_inbreeding
  ))
  cat("  (Both measured against the fixed founder allele frequencies.)\n")
  cat("  Cycle 0 is the founder population; cycle 1 onward are transmitted\n")
  cat("  selection responses.\n")
  if (!is.null(x$provenance$thread_warning)) {
    cat(
      "  WARNING: run with", x$provenance$n_threads,
      "threads and may not reproduce exactly.\n"
    )
  }
  invisible(x)
}
