# Construction of an AlphaSimR founder population from real phased haplotypes.
#
# DesiredGainR never simulates founder genomes from a coalescent model. The
# founders must come from the breeder's own phased marker data, so that the
# linkage disequilibrium structure, allele frequency spectrum and population
# structure of the simulation are those of the target programme rather than
# those of an assumed demography.
#
# This file deliberately does not depend on, import, or call HapBlockR. It
# accepts the data shape that HapBlockR's read_phased_vcf() returns, which
# keeps the dependency arrow one-directional: HapBlockR may use DesiredGainR,
# and DesiredGainR must therefore never use HapBlockR.

#' Require AlphaSimR
#'
#' @param context Short description of the calling operation.
#'
#' @noRd
.dgr_require_alphasimr <- function(context) {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop(
      context, " requires the AlphaSimR package. Install it with ",
      "install.packages(\"AlphaSimR\").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Validate one homologue matrix
#'
#' Missing values are permitted at this stage, because phased variant call
#' format output legitimately carries them. They are resolved later under an
#' explicit policy.
#'
#' @param x Candidate matrix.
#' @param name Object name used in messages.
#'
#' @return An integer matrix of 0, 1 and `NA`.
#' @noRd
.dgr_check_homologue <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(name, " must have variant row names and individual column names.",
      call. = FALSE
    )
  }
  observed <- unique(as.vector(x))
  unexpected <- setdiff(observed[!is.na(observed)], c(0L, 1L))
  if (length(unexpected)) {
    stop(
      name, " must contain only 0, 1 and NA, but also contains: ",
      paste(sort(unexpected), collapse = ", "),
      ". A haplotype is a single chromosome copy, so at a biallelic locus it ",
      "carries one allele rather than an allele count. A 0/1/2 matrix is ",
      "genotype dosage, which is the sum of the homologues and therefore ",
      "discards phase; see haplotypes_from_inbred_dosage() when the material ",
      "is highly inbred.",
      call. = FALSE
    )
  }
  x
}

#' Summarise heterozygosity and missingness in a dosage matrix
#'
#' No universally appropriate residual-heterozygosity threshold exists for
#' self-pollinated material, because the level depends on the generation,
#' mating history, crop, population type, genotyping error rate and
#' quality-control procedure. It must therefore be measured from each dataset
#' rather than assumed. Run this function before
#' [haplotypes_from_inbred_dosage()] and inspect the per-individual and
#' per-marker distributions, since a low overall rate can still conceal a small
#' number of badly affected individuals or markers.
#'
#' @param dosage Variant-by-individual matrix coded 0, 1, 2 with `NA` for
#'   missing calls.
#'
#' @return An object of class `desiredgainr_dosage_diagnostics` giving the
#'   overall heterozygous and missing rates together with per-individual and
#'   per-marker tables.
#'
#' @examples
#' dosage <- matrix(
#'   c(0, 2, 1, 0, 2, 0, 2, NA, 2),
#'   nrow = 3,
#'   dimnames = list(c("v1", "v2", "v3"), c("l1", "l2", "l3"))
#' )
#' dosage_diagnostics(dosage)
#'
#' @seealso [haplotypes_from_inbred_dosage()]
#' @export
dosage_diagnostics <- function(dosage) {
  dosage <- as.matrix(dosage)
  storage.mode(dosage) <- "integer"
  if (is.null(rownames(dosage)) || is.null(colnames(dosage))) {
    stop("dosage must have variant row names and individual column names.",
      call. = FALSE
    )
  }
  observed <- unique(as.vector(dosage))
  unexpected <- setdiff(observed[!is.na(observed)], c(0L, 1L, 2L))
  if (length(unexpected)) {
    stop("dosage must contain only 0, 1, 2 and NA, but also contains: ",
      paste(sort(unexpected), collapse = ", "),
      call. = FALSE
    )
  }

  called <- !is.na(dosage)
  n_called <- sum(called)
  heterozygous <- dosage == 1L & called
  safe_mean <- function(numerator, denominator) {
    ifelse(denominator > 0, numerator / denominator, NA_real_)
  }

  result <- list(
    n_variants = nrow(dosage),
    n_individuals = ncol(dosage),
    n_calls = n_called,
    overall_heterozygosity = if (n_called > 0) {
      sum(heterozygous) / n_called
    } else {
      NA_real_
    },
    overall_missing = mean(!called),
    by_individual = data.table::data.table(
      individual = colnames(dosage),
      heterozygosity = safe_mean(colSums(heterozygous), colSums(called)),
      missing = colMeans(!called)
    ),
    by_marker = data.table::data.table(
      variant = rownames(dosage),
      heterozygosity = safe_mean(rowSums(heterozygous), rowSums(called)),
      missing = rowMeans(!called)
    )
  )
  class(result) <- c("desiredgainr_dosage_diagnostics", "list")
  result
}

#' @export
print.desiredgainr_dosage_diagnostics <- function(x, ...) {
  cat("<desiredgainr_dosage_diagnostics>\n")
  cat(sprintf(
    "  %d variants x %d individuals\n", x$n_variants, x$n_individuals
  ))
  cat(sprintf(
    "  Heterozygous calls: %.3f%%   Missing calls: %.3f%%\n",
    100 * x$overall_heterozygosity, 100 * x$overall_missing
  ))
  cat(sprintf(
    "  Per-individual heterozygosity, maximum %.3f%% (%s)\n",
    100 * max(x$by_individual$heterozygosity, na.rm = TRUE),
    x$by_individual$individual[which.max(x$by_individual$heterozygosity)]
  ))
  cat(sprintf(
    "  Per-marker heterozygosity, maximum %.3f%%\n",
    100 * max(x$by_marker$heterozygosity, na.rm = TRUE)
  ))
  invisible(x)
}

#' Apply a drop-or-mask policy to a set of aligned matrices
#'
#' @param flagged Logical variant-by-individual matrix marking affected calls.
#' @param matrices List of aligned matrices to be filtered.
#' @param policy One of `"error"`, `"drop_variant"`, `"drop_individual"`,
#'   `"mask"`.
#' @param label Description of the affected calls, used in messages.
#'
#' @return A list with the filtered matrices and a record of what was removed.
#' @noRd
.dgr_apply_call_policy <- function(flagged, matrices, policy, label) {
  removed_variants <- character(0)
  removed_individuals <- character(0)
  n_flagged <- sum(flagged)

  if (n_flagged == 0L) {
    return(list(
      matrices = matrices,
      removed_variants = removed_variants,
      removed_individuals = removed_individuals,
      n_flagged = 0L
    ))
  }
  if (policy == "error") {
    stop(
      sprintf(
        paste(
          "%d %s cannot be resolved automatically. Choose a policy of",
          "'drop_variant', 'drop_individual' or 'mask', or resolve them",
          "upstream. Run dosage_diagnostics() to see whether the affected",
          "calls concentrate in a few markers or a few individuals."
        ),
        n_flagged, label
      ),
      call. = FALSE
    )
  }
  if (policy == "drop_variant") {
    keep <- rowSums(flagged) == 0L
    if (!any(keep)) {
      stop("Every variant carries ", label, "; none can be retained. ",
        "Consider policy = 'drop_individual'.",
        call. = FALSE
      )
    }
    removed_variants <- rownames(flagged)[!keep]
    matrices <- lapply(matrices, function(m) m[keep, , drop = FALSE])
  } else if (policy == "drop_individual") {
    keep <- colSums(flagged) == 0L
    if (sum(keep) < 2L) {
      stop("Fewer than two individuals remain after dropping those carrying ",
        label, ". Consider policy = 'drop_variant'.",
        call. = FALSE
      )
    }
    removed_individuals <- colnames(flagged)[!keep]
    matrices <- lapply(matrices, function(m) m[, keep, drop = FALSE])
  } else if (policy == "mask") {
    matrices <- lapply(matrices, function(m) {
      m[flagged] <- NA_integer_
      m
    })
  }
  list(
    matrices = matrices,
    removed_variants = removed_variants,
    removed_individuals = removed_individuals,
    n_flagged = n_flagged
  )
}

#' Derive phased haplotypes from inbred-line dosage
#'
#' Phase is required for simulation because it determines which favourable
#' alleles segregate together. However, in a diploid inbred line a dosage of 0
#' or 2 is already unambiguous: both homologues carry the same allele.
#' Therefore dosage from doubled haploids, recombinant inbred lines or advanced
#' selfed generations can be converted without external phasing, and only the
#' heterozygous and missing calls require a decision.
#'
#' Hence this helper serves self-pollinated diploid programmes, such as wheat,
#' rice and common bean, that hold dosage rather than phased data. It must not
#' be used for outcrossing or clonal material, where heterozygosity is
#' pervasive. Where heterozygosity is appreciable, or where cis and trans
#' relationships matter, phase the data externally instead.
#'
#' @details
#' No threshold is imposed on residual heterozygosity, because no universally
#' appropriate level exists. The rates are measured, reported and returned
#' instead, and heterozygous calls are never resolved silently: an explicit
#' policy is always required. Inspect [dosage_diagnostics()] first, since a low
#' overall rate can conceal a few badly affected individuals or markers, and
#' the cheaper of `"drop_variant"` and `"drop_individual"` depends on which.
#'
#' @param dosage Variant-by-individual matrix coded 0, 1, 2 with `NA` for
#'   missing calls.
#' @param heterozygous_policy Treatment of dosage 1. `"error"` refuses to
#'   proceed, `"drop_variant"` removes affected variants, `"drop_individual"`
#'   removes affected individuals, and `"mask"` sets the calls to `NA` for
#'   resolution by [founder_haplotypes()].
#' @param missing_policy Treatment of missing calls, with the same options.
#'
#' @return A list with `hap1` and `hap2` suitable for [founder_haplotypes()],
#'   carrying a `conversion` attribute that records the observed rates and
#'   everything removed or masked.
#'
#' @examples
#' dosage <- matrix(
#'   c(0, 2, 2, 0, 2, 0),
#'   nrow = 3,
#'   dimnames = list(c("v1", "v2", "v3"), c("line1", "line2"))
#' )
#' haplotypes_from_inbred_dosage(dosage)
#'
#' @seealso [dosage_diagnostics()], [founder_haplotypes()]
#' @export
haplotypes_from_inbred_dosage <- function(
  dosage,
  heterozygous_policy = c(
    "error", "drop_variant", "drop_individual", "mask"
  ),
  missing_policy = c(
    "error", "drop_variant", "drop_individual", "mask"
  )
) {
  heterozygous_policy <- match.arg(heterozygous_policy)
  missing_policy <- match.arg(missing_policy)
  diagnostics <- dosage_diagnostics(dosage)

  dosage <- as.matrix(dosage)
  storage.mode(dosage) <- "integer"

  missing_step <- .dgr_apply_call_policy(
    flagged = is.na(dosage),
    matrices = list(dosage = dosage),
    policy = missing_policy,
    label = "missing call(s)"
  )
  dosage <- missing_step$matrices$dosage

  heterozygous_step <- .dgr_apply_call_policy(
    flagged = !is.na(dosage) & dosage == 1L,
    matrices = list(dosage = dosage),
    policy = heterozygous_policy,
    label = "heterozygous call(s), which dosage alone cannot phase"
  )
  dosage <- heterozygous_step$matrices$dosage

  hap1 <- matrix(
    NA_integer_, nrow(dosage), ncol(dosage),
    dimnames = dimnames(dosage)
  )
  hap2 <- hap1
  homozygous_reference <- !is.na(dosage) & dosage == 0L
  homozygous_alternative <- !is.na(dosage) & dosage == 2L
  hap1[homozygous_reference] <- 0L
  hap2[homozygous_reference] <- 0L
  hap1[homozygous_alternative] <- 1L
  hap2[homozygous_alternative] <- 1L

  result <- list(hap1 = hap1, hap2 = hap2)
  attr(result, "conversion") <- list(
    diagnostics = diagnostics,
    heterozygous_policy = heterozygous_policy,
    missing_policy = missing_policy,
    n_heterozygous_calls = heterozygous_step$n_flagged,
    n_missing_calls = missing_step$n_flagged,
    removed_variants = unique(c(
      missing_step$removed_variants, heterozygous_step$removed_variants
    )),
    removed_individuals = unique(c(
      missing_step$removed_individuals, heterozygous_step$removed_individuals
    )),
    basis = paste(
      "Phase was determined from homozygous dosage only. No heterozygous call",
      "was assigned to a homologue, because dosage does not identify which",
      "alleles are in cis across loci."
    )
  )
  result
}

#' Validate and package phased founder haplotypes
#'
#' Multi-cycle simulation in DesiredGainR is calibrated to the breeder's own
#' germplasm rather than to an assumed genetic architecture. Hence the founder
#' population is always built from real phased marker data.
#'
#' @details
#' # Input shape
#'
#' `hap1` and `hap2` are variant-by-individual matrices coded 0, 1 or `NA`,
#' giving the allele carried on the first and second homologue. Row names must
#' be variant identifiers and column names individual identifiers, and both
#' matrices must carry identical dimnames in identical order. This is the shape
#' returned by `HapBlockR::read_phased_vcf()`, whose `hap1` and `hap2` can be
#' passed directly, with its `snp_info` supplying `map`.
#'
#' # Why a haplotype is coded 0 or 1, and dosage is not accepted
#'
#' A haplotype is one chromosome copy, so at a biallelic locus it carries a
#' single allele rather than an allele count. Genotype dosage coded 0, 1, 2 is
#' the sum of the homologues; it records how many alternative alleles an
#' individual carries but not which copy carries them, and therefore does not
#' identify which alleles lie in cis across loci. Since that co-occurrence is
#' the linkage disequilibrium the simulation exists to represent, dosage is
#' rejected. Where the material is highly inbred diploid, use
#' [haplotypes_from_inbred_dosage()], which derives phase from homozygous calls
#' alone.
#'
#' # Diploid only
#'
#' The contract is deliberately diploid and biallelic. Extending `hap1` and
#' `hap2` to further homologues would not make the simulation polyploid,
#' because autopolyploid and allopolyploid meiosis require an explicit pairing
#' and recombination model, including the treatment of multivalents and double
#' reduction, which differs by crop. Supplying a polyploid dataset therefore
#' raises an error rather than silently applying a diploid meiotic model.
#'
#' # Per-marker haplotypes, not block haplotype states
#'
#' Two haplotype representations are common in haplotype-aware pipelines, and
#' only the first belongs here.
#'
#' \describe{
#'   \item{Per-marker phased haplotypes}{One row per marker, coded 0 or 1,
#'     stating the allele carried on a given chromosome copy. A genome
#'     simulator models recombination between individual loci and therefore
#'     needs their separate positions.}
#'   \item{Block haplotype states}{One entry per linkage-disequilibrium block,
#'     recording which distinct haplotype an individual carries in that block.
#'     These are the correct representation for a haplotype relationship
#'     matrix, block association testing and local breeding values, but they
#'     have aggregated across the constituent markers and no longer provide a
#'     chromosome-wide homologue representation or a marker-level recombination
#'     map.}
#' }
#'
#' Read a phased variant call format file with a reader that preserves phase.
#' A general genotype reader may return dosage even for a phased file, in which
#' case the phase is discarded on import.
#'
#' @param hap1,hap2 Variant-by-individual homologue matrices coded 0, 1 or
#'   `NA`.
#' @param map Data frame describing the variants, requiring the variant
#'   identifier, the chromosome, and either a genetic or a physical position.
#' @param missing_policy Treatment of missing calls, which AlphaSimR cannot
#'   accept: `"error"`, `"drop_variant"` or `"drop_individual"`.
#' @param variant_col,chromosome_col Column names in `map`.
#' @param position_morgan_col Column giving genetic position in Morgans. When
#'   absent, `position_bp_col` and `recombination_rate` are used instead.
#' @param position_bp_col Column giving physical position in base pairs.
#' @param recombination_rate Morgans per base pair, used only when genetic
#'   positions are absent. The default corresponds to one centimorgan per
#'   megabase.
#' @param haplotypes Optional list of homologue matrices. Only a list of length
#'   two is accepted; see Details.
#'
#' @return An object of class `desiredgainr_founders` holding the
#'   per-chromosome genetic map and interleaved haplotype matrices required by
#'   AlphaSimR, together with the individual and variant identifiers and a
#'   record of any calls removed.
#'
#' @examples
#' variants <- c("v1", "v2", "v3", "v4")
#' individuals <- c("a", "b")
#' hap1 <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0),
#'   nrow = 4,
#'   dimnames = list(variants, individuals)
#' )
#' hap2 <- matrix(c(0, 0, 1, 1, 1, 1, 0, 0),
#'   nrow = 4,
#'   dimnames = list(variants, individuals)
#' )
#' map <- data.frame(
#'   variant_id = variants,
#'   chromosome = c(1, 1, 2, 2),
#'   position_bp = c(1e6, 5e6, 1e6, 4e6)
#' )
#' founder_haplotypes(hap1, hap2, map)
#'
#' @seealso [haplotypes_from_inbred_dosage()], [founder_population()]
#' @export
founder_haplotypes <- function(
  hap1 = NULL,
  hap2 = NULL,
  map,
  missing_policy = c("error", "drop_variant", "drop_individual"),
  variant_col = "variant_id",
  chromosome_col = "chromosome",
  position_morgan_col = "position_morgan",
  position_bp_col = "position_bp",
  recombination_rate = 1e-8,
  haplotypes = NULL
) {
  missing_policy <- match.arg(missing_policy)
  if (is.null(haplotypes)) {
    if (is.null(hap1) || is.null(hap2)) {
      stop("Supply hap1 and hap2, the two homologue matrices.", call. = FALSE)
    }
    haplotypes <- list(hap1 = hap1, hap2 = hap2)
  }
  if (!is.list(haplotypes)) {
    stop("haplotypes must be a list of homologue matrices.", call. = FALSE)
  }
  if (length(haplotypes) != 2L) {
    stop(
      "Only diploid material is supported: exactly two homologue matrices are ",
      "required, but ", length(haplotypes), " were supplied. Polyploid ",
      "simulation is not a matter of adding homologues, because autopolyploid ",
      "and allopolyploid meiosis require an explicit crop-specific pairing and ",
      "recombination model, including multivalent formation and double ",
      "reduction. Applying a diploid meiotic model to polyploid homologues ",
      "would produce confidently wrong results.",
      call. = FALSE
    )
  }
  homologue_names <- names(haplotypes)
  if (is.null(homologue_names)) homologue_names <- c("hap1", "hap2")
  haplotypes <- lapply(seq_len(2L), function(k) {
    .dgr_check_homologue(haplotypes[[k]], homologue_names[k])
  })
  if (!identical(dimnames(haplotypes[[1L]]), dimnames(haplotypes[[2L]]))) {
    stop(
      "Both homologue matrices must have identical dimnames in identical ",
      "order.",
      call. = FALSE
    )
  }

  # AlphaSimR requires complete haplotypes, so missing calls must be resolved
  # here under an explicit policy rather than imputed silently.
  missing_flag <- is.na(haplotypes[[1L]]) | is.na(haplotypes[[2L]])
  missing_rate <- mean(missing_flag)
  resolved <- .dgr_apply_call_policy(
    flagged = missing_flag,
    matrices = list(hap1 = haplotypes[[1L]], hap2 = haplotypes[[2L]]),
    policy = missing_policy,
    label = "missing haplotype call(s), which AlphaSimR cannot accept"
  )
  haplotypes <- list(resolved$matrices$hap1, resolved$matrices$hap2)

  variant_id <- rownames(haplotypes[[1L]])
  individual_id <- colnames(haplotypes[[1L]])
  if (anyDuplicated(variant_id) || anyDuplicated(individual_id)) {
    stop("Variant and individual identifiers must be unique.", call. = FALSE)
  }
  if (length(individual_id) < 2L) {
    stop("At least two founder individuals are required.", call. = FALSE)
  }

  map <- as.data.frame(map, stringsAsFactors = FALSE)
  absent <- setdiff(c(variant_col, chromosome_col), names(map))
  if (length(absent)) {
    stop(
      "map is missing columns: ", paste(absent, collapse = ", "),
      ". When using HapBlockR::read_phased_vcf() output, pass its snp_info ",
      "table and set variant_col, chromosome_col and position_bp_col to the ",
      "corresponding column names.",
      call. = FALSE
    )
  }
  missing_variants <- setdiff(variant_id, map[[variant_col]])
  if (length(missing_variants)) {
    stop(
      "map is missing ", length(missing_variants), " variant(s), including: ",
      paste(utils::head(missing_variants, 5L), collapse = ", "),
      call. = FALSE
    )
  }
  map <- map[match(variant_id, map[[variant_col]]), , drop = FALSE]

  has_morgan <- !is.null(position_morgan_col) &&
    position_morgan_col %in% names(map)
  if (has_morgan) {
    position <- as.numeric(map[[position_morgan_col]])
    position_source <- "supplied genetic position in Morgans"
  } else {
    if (!position_bp_col %in% names(map)) {
      stop("map must supply either ", position_morgan_col, " or ",
        position_bp_col, ".",
        call. = FALSE
      )
    }
    if (!is.numeric(recombination_rate) ||
      length(recombination_rate) != 1L ||
      !is.finite(recombination_rate) || recombination_rate <= 0) {
      stop("recombination_rate must be a positive finite scalar.",
        call. = FALSE
      )
    }
    position <- as.numeric(map[[position_bp_col]]) * recombination_rate
    position_source <- sprintf(
      "converted from physical position at %.3g Morgan per base pair",
      recombination_rate
    )
  }
  if (any(!is.finite(position))) {
    stop("Variant positions must be finite.", call. = FALSE)
  }
  chromosome <- as.character(map[[chromosome_col]])
  if (anyNA(chromosome) || any(!nzchar(chromosome))) {
    stop("Every variant must carry a non-missing chromosome label.",
      call. = FALSE
    )
  }

  chromosomes <- unique(chromosome)
  n_individuals <- length(individual_id)
  gen_map <- vector("list", length(chromosomes))
  interleaved <- vector("list", length(chromosomes))
  names(gen_map) <- names(interleaved) <- chromosomes

  for (chr in chromosomes) {
    rows <- which(chromosome == chr)
    rows <- rows[order(position[rows])]
    if (length(rows) < 2L) {
      stop("Chromosome ", chr, " carries fewer than two variants; AlphaSimR ",
        "requires at least two loci per chromosome.",
        call. = FALSE
      )
    }
    chromosome_position <- position[rows] - min(position[rows])
    if (anyDuplicated(chromosome_position)) {
      # Break exact ties by a negligible offset, because AlphaSimR requires a
      # strictly increasing map.
      chromosome_position <- chromosome_position +
        seq_along(chromosome_position) * 1e-12
    }
    names(chromosome_position) <- variant_id[rows]
    gen_map[[chr]] <- chromosome_position

    # AlphaSimR expects homologues grouped by individual: individual 1
    # homologue 1, individual 1 homologue 2, individual 2 homologue 1, and so
    # on.
    block <- matrix(
      0L,
      nrow = 2L * n_individuals, ncol = length(rows),
      dimnames = list(NULL, variant_id[rows])
    )
    block[seq(1L, 2L * n_individuals, by = 2L), ] <-
      t(haplotypes[[1L]][rows, , drop = FALSE])
    block[seq(2L, 2L * n_individuals, by = 2L), ] <-
      t(haplotypes[[2L]][rows, , drop = FALSE])
    interleaved[[chr]] <- block
  }

  result <- list(
    gen_map = gen_map,
    haplotypes = interleaved,
    ploidy = 2L,
    individual_id = individual_id,
    variant_id = variant_id,
    n_individuals = n_individuals,
    n_variants = length(variant_id),
    n_chromosomes = length(chromosomes),
    chromosomes = chromosomes,
    missing_data = list(
      policy = missing_policy,
      rate_before_resolution = missing_rate,
      n_calls_affected = resolved$n_flagged,
      removed_variants = resolved$removed_variants,
      removed_individuals = resolved$removed_individuals
    ),
    provenance = list(
      position_source = position_source,
      note = paste(
        "Founder haplotypes supplied by the user. DesiredGainR does not",
        "simulate founder genomes, so the linkage disequilibrium and allele",
        "frequency structure of the simulation are those of the supplied",
        "germplasm."
      )
    )
  )
  class(result) <- c("desiredgainr_founders", "list")
  result
}

#' @export
print.desiredgainr_founders <- function(x, ...) {
  cat("<desiredgainr_founders>\n")
  cat(sprintf(
    "  %d individuals, %d variants, %d chromosomes, diploid\n",
    x$n_individuals, x$n_variants, x$n_chromosomes
  ))
  cat("  Map:", x$provenance$position_source, "\n")
  cat(sprintf(
    "  Missing calls before resolution: %.3f%% (policy: %s)\n",
    100 * x$missing_data$rate_before_resolution, x$missing_data$policy
  ))
  if (length(x$missing_data$removed_variants)) {
    cat(sprintf(
      "  Variants removed: %d\n",
      length(x$missing_data$removed_variants)
    ))
  }
  if (length(x$missing_data$removed_individuals)) {
    cat(sprintf(
      "  Individuals removed: %d\n",
      length(x$missing_data$removed_individuals)
    ))
  }
  invisible(x)
}

#' Build an AlphaSimR founder population with a target genetic covariance
#'
#' The founder genomes come from the breeder's phased marker data through
#' [founder_haplotypes()], whereas the trait architecture is calibrated to the
#' genetic covariance matrix the breeder has already estimated. Therefore the
#' simulation reproduces both the observed germplasm structure and the observed
#' trait variances and correlations, rather than an assumed demography.
#'
#' @details
#' Trait variances are taken from the diagonal of `G` and trait correlations
#' from the corresponding correlation matrix.
#'
#' # Additive and clonal programmes calibrate differently
#'
#' For a self-pollinated or cross-pollinated programme, `G` is the additive
#' genetic covariance and is passed straight to AlphaSimR.
#'
#' For a clonal programme, supply `dominance_degree`. `G` is then the total
#' **genotypic** covariance, because the unit of selection is the clone and
#' dominance is inherited intact. This needs a calibration step that earlier
#' versions omitted. AlphaSimR's `var` argument sets the *additive* variance,
#' so with dominance present the dominance variance is added on top and the
#' realised genotypic variance exceeds the supplied target. Since scaling every
#' quantitative trait locus effect by \eqn{c} scales both the additive and the
#' dominance variance by \eqn{c^2}, the correction is exact in one pass: the
#' trait is built, the realised genotypic variance measured, and the additive
#' target rescaled by the ratio.
#'
#' The correlations are calibrated too, by the same fixed point. Rescaling the
#' variances alone would leave the realised genotypic correlations emergent, so
#' the setup would accept a full covariance matrix as its target and reproduce
#' only its diagonal. Dominance perturbs the correlation structure as well, and
#' that perturbation is a smooth, near-identity function of the additive
#' correlations supplied, so correcting both and projecting the correction back
#' onto the set of valid correlation matrices recovers the whole of `G`.
#'
#' `setup$G_realised` records what was achieved and `setup$calibration_error`
#' the largest remaining deviation in each of the variances and correlations.
#' Convergence is not guaranteed for every target: a strongly antagonistic
#' correlation structure at a high dominance degree may not be attainable by
#' any additive-plus-dominance architecture. Failure warns and is recorded in
#' `setup$calibration_converged` rather than passing silently.
#'
#' # Which heritability
#'
#' `heritability = "narrow"` sets the error variance so that
#' \eqn{V_A / V_P} equals `h2`; `"broad"` targets \eqn{V_G / V_P}. They
#' coincide without dominance and diverge with it, so a clonal programme must
#' say which is meant. Selection in a clonal programme acts on genotypic value,
#' making broad sense usually the relevant one.
#'
#' # Markers for coancestry
#'
#' `n_markers_per_chromosome` creates a neutral marker panel, held disjoint
#' from the quantitative trait loci where the installed AlphaSimR supports it.
#' Diversity should not be measured on the loci under selection: changing their
#' frequencies is what genetic gain is, so a relationship matrix built from
#' them rises whenever selection succeeds, and a desired-gain direction would
#' be penalised for working.
#'
#' @param founders An object from [founder_haplotypes()].
#' @param G Genetic variance-covariance matrix, named by trait. Without
#'   dominance this is the additive covariance. With dominance it is the total
#'   **genotypic** covariance, because that is what a clonal programme selects
#'   on; see Details.
#' @param h2 Named heritabilities, or a single value applied to every trait.
#'   Narrow- or broad-sense according to `heritability`.
#' @param residual_covariance Optional named residual covariance matrix. When
#'   supplied it is passed directly to AlphaSimR and takes precedence over the
#'   marginal `h2` values. This is the route used when propagating a sampled
#'   residual covariance; `h2` remains recorded as the nominal setup input.
#' @param n_qtl_per_chromosome Number of quantitative trait loci simulated per
#'   chromosome.
#' @param n_markers_per_chromosome Optional number of neutral markers per
#'   chromosome, used for coancestry. When supplied, the markers are held
#'   disjoint from the quantitative trait loci where the installed AlphaSimR
#'   supports it. Strongly recommended: without a panel, diversity is measured
#'   on all segregating sites.
#' @param dominance_degree Optional named mean degree of dominance per trait.
#'   Supplying it adds dominance to the simulated traits, which the clonal
#'   mating system requires.
#' @param dominance_variance Optional named variance of the **degree of
#'   dominance** across loci. This is AlphaSimR's `varDD` and is not the
#'   dominance genetic variance; the two are different quantities.
#' @param heritability Whether `h2` is narrow-sense (the default, a statement
#'   about breeding values) or broad-sense (a statement about genotypic values).
#'   They differ exactly when dominance is simulated, which is when a clonal
#'   programme is being represented, so it must be stated rather than assumed.
#' @param seed Random seed. The caller's random number generator state is
#'   restored on exit.
#'
#' @return A list of class `desiredgainr_sim_setup` holding the AlphaSimR
#'   simulation parameters, the founder population, and the calibration
#'   targets.
#'
#' @seealso [founder_haplotypes()], [simulate_selection_cycles()]
#' @export
founder_population <- function(
  founders,
  G,
  h2,
  residual_covariance = NULL,
  n_qtl_per_chromosome = 100L,
  n_markers_per_chromosome = NULL,
  dominance_degree = NULL,
  dominance_variance = NULL,
  heritability = c("narrow", "broad"),
  seed = 42L
) {
  .dgr_require_alphasimr("founder_population()")
  if (!inherits(founders, "desiredgainr_founders")) {
    stop("founders must be created by founder_haplotypes().", call. = FALSE)
  }
  trait_cols <- colnames(G)
  if (is.null(trait_cols) || anyDuplicated(trait_cols)) {
    stop("G must have unique trait dimnames.", call. = FALSE)
  }
  G <- .dgr_covariance(G, trait_cols, "G")
  .dgr_check_psd(G, "G")
  if (!is.null(residual_covariance)) {
    residual_covariance <- .dgr_covariance(
      residual_covariance, trait_cols, "residual_covariance"
    )
    .dgr_check_psd(residual_covariance, "residual_covariance")
  }
  n_traits <- length(trait_cols)

  if (length(h2) == 1L) {
    h2 <- stats::setNames(rep(as.numeric(h2), n_traits), trait_cols)
  }
  h2 <- .dgr_named_vector(h2, trait_cols, "h2")
  if (any(h2 <= 0 | h2 > 1)) {
    stop("Every heritability must lie in (0, 1].", call. = FALSE)
  }
  n_qtl_per_chromosome <- .dgr_positive_integer(
    n_qtl_per_chromosome, "n_qtl_per_chromosome"
  )
  if (!is.null(n_markers_per_chromosome)) {
    n_markers_per_chromosome <- .dgr_positive_integer(
      n_markers_per_chromosome, "n_markers_per_chromosome"
    )
  }
  seed <- .dgr_seed(seed)

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  entry_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", entry_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)

  heritability <- match.arg(heritability)
  use_dominance <- !is.null(dominance_degree)
  if (!use_dominance && identical(heritability, "broad")) {
    stop("heritability = 'broad' is only meaningful when dominance is ",
      "simulated. Without dominance the two coincide, so use 'narrow'.",
      call. = FALSE
    )
  }

  build <- function(additive_target, correlation_target = correlations) {
    # The seed is reset on every call so that each build draws the SAME
    # quantitative trait loci, positions and dominance degrees, and differs
    # only in the additive scale. Without this the calibration below would
    # measure one trait architecture and apply the correction to a different
    # one, chasing a target that moves each pass.
    set.seed(seed)
    map_pop <- AlphaSimR::newMapPop(
      genMap = founders$gen_map,
      haplotypes = founders$haplotypes
    )
    SP <- AlphaSimR::SimParam$new(map_pop)
    SP$setTrackPed(TRUE)

    # Markers used for coancestry must not be the loci under selection, or the
    # diversity metric rises whenever selection succeeds. Ask AlphaSimR to
    # partition the segregating sites so that the SNP chip and the QTL are
    # disjoint. Older versions lack the argument, in which case some overlap
    # is possible and is recorded in the setup.
    marker_overlap <- TRUE
    if (!is.null(n_markers_per_chromosome)) {
      restricted <- tryCatch(
        {
          SP$restrSegSites(
            minQtlPerChr = n_qtl_per_chromosome,
            minSnpPerChr = n_markers_per_chromosome,
            overlap = FALSE
          )
          TRUE
        },
        error = function(e) FALSE
      )
      marker_overlap <- !restricted
      tryCatch(
        SP$addSnpChip(nSnpPerChr = n_markers_per_chromosome),
        error = function(e) {
          stop("Could not create a marker panel of ",
            n_markers_per_chromosome, " markers per chromosome: ",
            conditionMessage(e),
            "\n  Reduce n_markers_per_chromosome or supply denser ",
            "founder haplotypes.",
            call. = FALSE
          )
        }
      )
    }

    if (use_dominance) {
      SP$addTraitAD(
        nQtlPerChr = n_qtl_per_chromosome,
        mean = rep(0, n_traits),
        var = as.numeric(additive_target),
        meanDD = as.numeric(dominance_degree),
        varDD = as.numeric(dominance_variance),
        corA = correlation_target
      )
    } else {
      SP$addTraitA(
        nQtlPerChr = n_qtl_per_chromosome,
        mean = rep(0, n_traits),
        var = as.numeric(additive_target),
        corA = correlation_target
      )
    }
    population <- AlphaSimR::newPop(map_pop, simParam = SP)
    list(SP = SP, pop = population, marker_overlap = marker_overlap)
  }

  variances <- diag(G)
  correlations <- stats::cov2cor(G)
  if (use_dominance) {
    dominance_degree <- .dgr_named_vector(
      dominance_degree, trait_cols, "dominance_degree"
    )
    # AlphaSimR's varDD is the variance of the dominance *degree* across loci,
    # not the dominance genetic variance. The two are different quantities and
    # conflating them was the source of the miscalibration below.
    dominance_variance <- if (is.null(dominance_variance)) {
      stats::setNames(rep(0.2, n_traits), trait_cols)
    } else {
      .dgr_named_vector(dominance_variance, trait_cols, "dominance_variance")
    }
  }

  built <- build(variances)
  realised <- stats::cov(AlphaSimR::gv(built$pop))
  dimnames(realised) <- list(trait_cols, trait_cols)

  calibration_note <- NULL
  calibration_iterations <- 0L
  calibration_converged <- NA
  if (use_dominance) {
    # AlphaSimR's `var` argument sets the ADDITIVE genetic variance. With
    # dominance present the dominance variance is added on top, so the total
    # genotypic variance exceeds the supplied target -- which is the matrix a
    # clonal programme selects on.
    #
    # Because build() reseeds, the trait architecture is held fixed and only
    # the additive scale varies, so the ratio correction converges quickly.
    # It is applied iteratively rather than once because the relationship
    # between the additive target and the realised genotypic variance is not
    # exactly quadratic: AlphaSimR rescales to hit the requested additive
    # variance after the dominance deviations are formed.
    # Both the variances and the correlations are calibrated. Rescaling the
    # additive variances alone leaves the realised genotypic CORRELATIONS
    # emergent, so the setup would accept a full covariance matrix as its
    # target and reproduce only its diagonal. Dominance perturbs the
    # correlation structure as well, and that perturbation is itself a smooth,
    # near-identity function of the additive correlations supplied, so a
    # fixed-point iteration on both recovers the whole matrix.
    variance_target <- variances
    correlation_target <- correlations
    additive_variance <- variances
    additive_correlation <- correlations
    tolerance <- 1e-3
    trace <- numeric(0)
    for (iteration in seq_len(40L)) {
      realised_variance <- diag(realised)
      if (any(!is.finite(realised_variance)) || any(realised_variance <= 0)) {
        stop("A trait had non-positive realised genotypic variance, so the ",
          "dominance calibration cannot be scaled.",
          call. = FALSE
        )
      }
      realised_correlation <- stats::cov2cor(realised)
      variance_error <- max(
        abs(realised_variance - variance_target) / variance_target
      )
      correlation_error <- max(abs(realised_correlation - correlation_target))
      trace <- c(trace, max(variance_error, correlation_error))
      calibration_iterations <- iteration - 1L
      if (variance_error < tolerance && correlation_error < tolerance) break

      additive_variance <- additive_variance *
        (variance_target / realised_variance)
      # Newton-style correction on the correlations, projected back onto the
      # set of valid correlation matrices so the next AlphaSimR call receives
      # an admissible corA.
      additive_correlation <- .dgr_project_correlation(
        additive_correlation + (correlation_target - realised_correlation)
      )
      built <- build(additive_variance, additive_correlation)
      realised <- stats::cov(AlphaSimR::gv(built$pop))
      dimnames(realised) <- list(trait_cols, trait_cols)
    }
    final_variance_error <- max(
      abs(diag(realised) - variance_target) / variance_target
    )
    final_correlation_error <- max(abs(
      stats::cov2cor(realised) - correlation_target
    ))
    calibration_converged <- final_variance_error < tolerance &&
      final_correlation_error < tolerance
    if (!isTRUE(calibration_converged)) {
      warning(
        "The genotypic covariance could not be calibrated to G within ",
        format(100 * tolerance, digits = 2), "%. Largest remaining deviation: ",
        format(100 * final_variance_error, digits = 3), "% on the variances, ",
        format(final_correlation_error, digits = 3), " on the correlations. ",
        "The usual cause is too few quantitative trait loci for the dominance ",
        "model to be scaled smoothly, or a target correlation structure that ",
        "the additive-plus-dominance architecture cannot attain at this ",
        "dominance degree. Increase n_qtl_per_chromosome, reduce ",
        "dominance_degree, and use setup$G_realised rather than G as the ",
        "truth for this simulation.",
        call. = FALSE
      )
    }
    calibration_note <- paste(
      "With dominance simulated, G is interpreted as the total GENOTYPIC",
      "covariance, because that is what a clonal programme selects on. Both",
      "the additive variances and the additive correlations were calibrated",
      "by fixed-point iteration over", calibration_iterations, "passes, so the",
      "realised genotypic covariance reproduces the whole of G rather than",
      "only its diagonal. G_realised records what was achieved."
    )
  }

  SP <- built$SP
  population <- built$pop
  # Narrow-sense heritability is a statement about breeding values and broad
  # sense about genotypic values. They differ exactly when dominance is
  # present, which is when a clonal programme is being represented, so the
  # caller must say which one h2 means rather than have it assumed.
  if (!is.null(residual_covariance)) {
    SP$setVarE(varE = residual_covariance)
  } else if (identical(heritability, "broad")) {
    SP$setVarE(H2 = as.numeric(h2))
  } else {
    SP$setVarE(h2 = as.numeric(h2))
  }
  # Re-phenotype the existing founders rather than creating a second
  # population, which would duplicate the identifiers already issued and
  # advance lastId for no reason.
  population <- AlphaSimR::setPheno(population, simParam = SP)

  marker_panel <- tryCatch(
    colnames(AlphaSimR::pullSnpGeno(population, simParam = SP)),
    error = function(e) NULL
  )

  result <- list(
    SP = SP,
    founder_pop = population,
    trait_cols = trait_cols,
    G_target = G,
    G_realised = realised,
    calibration_iterations = calibration_iterations,
    calibration_converged = calibration_converged,
    calibration_error = if (use_dominance) {
      list(
        variance = max(abs(diag(realised) - variances) / variances),
        correlation = max(abs(
          stats::cov2cor(realised) - stats::cov2cor(G)
        )),
        trace = trace
      )
    } else {
      NULL
    },
    h2 = h2,
    residual_covariance = residual_covariance,
    heritability_type = heritability,
    ploidy = 2L,
    n_qtl_per_chromosome = n_qtl_per_chromosome,
    n_markers_per_chromosome = n_markers_per_chromosome,
    marker_panel = marker_panel,
    marker_qtl_overlap = built$marker_overlap,
    dominance = use_dominance,
    dominance_degree = if (use_dominance) dominance_degree else NULL,
    dominance_degree_variance = if (use_dominance) dominance_variance else NULL,
    seed = seed,
    founders = founders,
    calibration = paste(
      c(
        paste(
          "Trait variances and correlations were calibrated to the supplied G.",
          "Genome structure came from the supplied phased haplotypes and was",
          "not simulated."
        ),
        calibration_note,
        if (is.null(marker_panel)) {
          paste(
            "No marker panel was requested, so diversity will be measured on",
            "all segregating sites. Supply n_markers_per_chromosome for a",
            "panel disjoint from the quantitative trait loci."
          )
        } else {
          paste(
            "Diversity is measured on a panel of", length(marker_panel),
            "markers",
            if (isTRUE(built$marker_overlap)) {
              "which may overlap the quantitative trait loci."
            } else {
              "held disjoint from the quantitative trait loci."
            }
          )
        }
      ),
      collapse = " "
    )
  )
  class(result) <- c("desiredgainr_sim_setup", "list")
  result
}

# Project a matrix onto the set of valid correlation matrices.
#
# Used inside the clonal calibration, where a Newton-style correction to the
# additive correlations can leave the matrix indefinite or off the unit
# diagonal. AlphaSimR requires an admissible corA, so the correction is
# projected before it is used.
#
# This is emphatically NOT bend_covariance(). That function repairs a user's
# estimate and reports the repair, because the adjustment changes what the data
# are taken to say. This one keeps an internal iterate inside its feasible set;
# nothing about the user's inputs is altered, and the quantity that matters,
# the realised genotypic covariance, is measured afterwards and reported as
# G_realised.
#
# Written as plain comments rather than roxygen deliberately: an internal
# helper needs no manual page, and a roxygen block placed immediately after
# another block's @export tag merges with it, so the prose becomes export
# names. That is how this file failed to install.
.dgr_project_correlation <- function(x, min_eigenvalue = 1e-6) {
  x <- (x + t(x)) / 2
  diag(x) <- 1
  # Alternate between flooring the eigenvalues and restoring the unit
  # diagonal. Each step alone can undo the other, so a few passes are needed;
  # this converges quickly because the input is already close to admissible.
  for (pass in seq_len(20L)) {
    decomposition <- eigen(x, symmetric = TRUE)
    if (min(decomposition$values) >= min_eigenvalue &&
      max(abs(diag(x) - 1)) < 1e-10) {
      break
    }
    values <- pmax(decomposition$values, min_eigenvalue)
    x <- decomposition$vectors %*% diag(values, nrow(x)) %*%
      t(decomposition$vectors)
    x <- (x + t(x)) / 2
    scaling <- sqrt(diag(x))
    scaling[!is.finite(scaling) | scaling <= 0] <- 1
    x <- sweep(sweep(x, 1L, scaling, "/"), 2L, scaling, "/")
    diag(x) <- 1
  }
  # Numerical drift can leave an entry marginally outside [-1, 1].
  x[x > 1] <- 1
  x[x < -1] <- -1
  diag(x) <- 1
  x
}

#' @export
print.desiredgainr_sim_setup <- function(x, ...) {
  cat("<desiredgainr_sim_setup>\n")
  cat(sprintf(
    "  Founders: %d individuals, %d chromosomes, %d QTL per chromosome\n",
    x$founders$n_individuals, x$founders$n_chromosomes,
    x$n_qtl_per_chromosome
  ))
  cat("  Traits:", paste(x$trait_cols, collapse = ", "), "\n")
  cat("  Dominance simulated:", if (isTRUE(x$dominance)) "yes" else "no", "\n")
  cat("  Heritability supplied as:", x$heritability_type, "sense\n")
  if (is.null(x$marker_panel)) {
    cat("  Marker panel: none (diversity uses all segregating sites)\n")
  } else {
    cat(sprintf(
      "  Marker panel: %d markers%s\n", length(x$marker_panel),
      if (isTRUE(x$marker_qtl_overlap)) " (may overlap QTL)" else ", QTL-free"
    ))
  }
  if (isTRUE(x$dominance)) {
    cat(sprintf(
      "  Genotypic covariance calibrated in %d passes (%s)\n",
      x$calibration_iterations,
      if (isTRUE(x$calibration_converged)) "converged" else "NOT CONVERGED"
    ))
    cat(sprintf(
      "    largest deviation: %.3g on variances, %.3g on correlations\n",
      x$calibration_error$variance, x$calibration_error$correlation
    ))
  }
  invisible(x)
}
