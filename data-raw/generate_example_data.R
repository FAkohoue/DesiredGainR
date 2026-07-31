## data-raw/generate_example_data.R
## ---------------------------------------------------------------------------
## Generates the example datasets shipped with DesiredGainR.
## Run once from the package root with:
##   source("data-raw/generate_example_data.R")
##
## Why these data exist
## --------------------
## Earlier releases demonstrated the package on uncorrelated standard normal
## variates, so the genetic covariance matrix was effectively an identity
## matrix. That cannot illustrate a multi-trait selection index, because the
## entire purpose of an index is to resolve the tension between correlated
## traits measured on different scales.
##
## These data therefore have three properties the previous examples lacked:
##
##   (i)   genuinely antagonistic genetic correlations, notably a positive
##         correlation between grain yield and anthesis date when the breeding
##         objective requires yield to rise and anthesis date to fall;
##   (ii)  trait scales spanning two orders of magnitude, from ears per plant
##         near 1.0 to plant height near 205 cm, which is what makes
##         standardisation and matrix conditioning matter; and
##   (iii) one internally consistent population, so that the markers actually
##         explain the trait values rather than being generated separately.
##
## The material is simulated. It is patterned on a tropical maize breeding
## programme but represents no real germplasm and no real trial.
##
## Datasets produced
## -----------------
##   dgr_traits      trait metadata: unit, direction, heritability
##   dgr_G           6 x 6 additive genetic covariance matrix
##   dgr_P           6 x 6 phenotypic covariance matrix, equal to G + E
##   dgr_candidates  200 candidates x 6 adjusted means, with family labels
##   dgr_gebv        200 candidates x 6 genomic estimated breeding values
##   dgr_history     a previous cycle's selection decision, for
##                   retrospective_weights()
##   dgr_hap1        450 markers x 200 candidates, first homologue, 0 or 1
##   dgr_hap2        450 markers x 200 candidates, second homologue, 0 or 1
##   dgr_map         marker map: identifier, chromosome, physical position
## ---------------------------------------------------------------------------

set.seed(20260730)

# -- Programme definition -----------------------------------------------------

TRAITS <- c("GY", "PHT", "AD", "ASI", "EPP", "GLS")

dgr_traits <- data.frame(
  trait = TRAITS,
  description = c(
    "Grain yield",
    "Plant height",
    "Anthesis date",
    "Anthesis-silking interval",
    "Ears per plant",
    "Grey leaf spot severity"
  ),
  unit = c("t/ha", "cm", "days", "days", "count", "score 1-9"),
  direction = c(
    "increase", "decrease", "decrease", "decrease", "increase", "decrease"
  ),
  heritability = c(0.35, 0.60, 0.70, 0.30, 0.40, 0.45),
  mean = c(5.20, 205.0, 66.0, 1.80, 0.98, 3.90),
  genetic_sd = c(0.75, 12.0, 2.50, 0.60, 0.10, 0.85),
  stringsAsFactors = FALSE
)
rownames(dgr_traits) <- NULL
N_TRAITS <- length(TRAITS)

## Additive genetic correlations. The structure is deliberately awkward:
## yield correlates positively with anthesis date, yet the objective is to
## raise yield while shortening the cycle, so the two cannot both be pushed
## freely. Yield also correlates strongly and negatively with the
## anthesis-silking interval, which is the standard drought-stress indicator.
genetic_correlation <- matrix(
  c(
    #  GY    PHT     AD    ASI    EPP    GLS
     1.00,  0.25,  0.30, -0.55,  0.45, -0.35,
     0.25,  1.00,  0.35, -0.10,  0.15, -0.20,
     0.30,  0.35,  1.00,  0.20,  0.05, -0.15,
    -0.55, -0.10,  0.20,  1.00, -0.30,  0.25,
     0.45,  0.15,  0.05, -0.30,  1.00, -0.20,
    -0.35, -0.20, -0.15,  0.25, -0.20,  1.00
  ),
  nrow = N_TRAITS, byrow = TRUE,
  dimnames = list(TRAITS, TRAITS)
)

## Residual correlations are weaker and do not mirror the genetic ones, which
## is what makes the phenotypic and genetic covariance matrices genuinely
## different objects rather than one rescaled into the other.
residual_correlation <- diag(N_TRAITS)
dimnames(residual_correlation) <- list(TRAITS, TRAITS)
residual_pairs <- list(
  c("GY", "EPP", 0.20), c("GY", "PHT", 0.10),
  c("PHT", "AD", 0.15), c("AD", "ASI", 0.10)
)
for (pair in residual_pairs) {
  residual_correlation[pair[1], pair[2]] <- as.numeric(pair[3])
  residual_correlation[pair[2], pair[1]] <- as.numeric(pair[3])
}

# -- Helper: project a correlation matrix to the nearest positive definite ----
# A correlation matrix written down by hand need not be positive definite, and
# a covariance matrix that is not positive definite cannot be inverted. The
# projection is reported so that any adjustment is visible rather than silent.
nearest_positive_definite <- function(M, minimum_eigenvalue = 1e-6,
                                      label = "matrix") {
  decomposition <- eigen((M + t(M)) / 2, symmetric = TRUE)
  if (min(decomposition$values) >= minimum_eigenvalue) {
    message(label, ": already positive definite (smallest eigenvalue ",
            format(min(decomposition$values), digits = 3), ")")
    return(M)
  }
  message(label, ": adjusted to positive definite (smallest eigenvalue was ",
          format(min(decomposition$values), digits = 3), ")")
  values <- pmax(decomposition$values, minimum_eigenvalue)
  adjusted <- decomposition$vectors %*% diag(values) %*% t(decomposition$vectors)
  adjusted <- stats::cov2cor(adjusted)
  dimnames(adjusted) <- dimnames(M)
  adjusted
}

genetic_correlation <- nearest_positive_definite(
  genetic_correlation, label = "Genetic correlation"
)
residual_correlation <- nearest_positive_definite(
  residual_correlation, label = "Residual correlation"
)

## Covariance matrices. The residual variance follows from the declared
## heritability, so the three quantities stay mutually consistent.
genetic_sd <- dgr_traits$genetic_sd
names(genetic_sd) <- TRAITS
residual_sd <- genetic_sd * sqrt((1 - dgr_traits$heritability) /
                                   dgr_traits$heritability)
names(residual_sd) <- TRAITS

dgr_G <- diag(genetic_sd) %*% genetic_correlation %*% diag(genetic_sd)
dimnames(dgr_G) <- list(TRAITS, TRAITS)
residual_covariance <- diag(residual_sd) %*% residual_correlation %*%
  diag(residual_sd)
dimnames(residual_covariance) <- list(TRAITS, TRAITS)
dgr_P <- dgr_G + residual_covariance
dimnames(dgr_P) <- list(TRAITS, TRAITS)

stopifnot(
  all(abs(diag(dgr_G) / diag(dgr_P) - dgr_traits$heritability) < 1e-8)
)
message("Heritabilities reproduced from G and P: ",
        paste(sprintf("%s=%.2f", TRAITS, diag(dgr_G) / diag(dgr_P)),
              collapse = ", "))

# -- Marker panel with linkage disequilibrium ---------------------------------

N_IND <- 200L
N_CHR <- 3L
MARKERS_PER_BLOCK <- 15L
BLOCKS_PER_CHR <- 10L
N_FOUNDER_HAPLOTYPES <- 8L

markers_per_chr <- MARKERS_PER_BLOCK * BLOCKS_PER_CHR
n_markers <- markers_per_chr * N_CHR

## Each block carries a small pool of founder haplotypes, and every individual
## draws one pool member per homologue. Markers within a block therefore travel
## together, which is what creates linkage disequilibrium; markers in different
## blocks segregate independently.
draw_homologue <- function() {
  homologue <- matrix(0L, nrow = n_markers, ncol = N_IND)
  marker_offset <- 0L
  for (chromosome in seq_len(N_CHR)) {
    for (block in seq_len(BLOCKS_PER_CHR)) {
      pool <- matrix(
        rbinom(N_FOUNDER_HAPLOTYPES * MARKERS_PER_BLOCK, 1L,
               runif(1L, 0.25, 0.75)),
        nrow = N_FOUNDER_HAPLOTYPES
      )
      chosen <- sample.int(N_FOUNDER_HAPLOTYPES, N_IND, replace = TRUE)
      rows <- marker_offset + seq_len(MARKERS_PER_BLOCK)
      homologue[rows, ] <- t(pool[chosen, , drop = FALSE])
      marker_offset <- marker_offset + MARKERS_PER_BLOCK
    }
  }
  homologue
}

candidate_id <- sprintf("CAND%03d", seq_len(N_IND))
marker_id <- sprintf("M%04d", seq_len(n_markers))

dgr_hap1 <- draw_homologue()
dgr_hap2 <- draw_homologue()
dimnames(dgr_hap1) <- dimnames(dgr_hap2) <- list(marker_id, candidate_id)

## Discard markers that came out monomorphic, since they carry no information
## and would only complicate the worked examples.
dosage <- dgr_hap1 + dgr_hap2
allele_frequency <- rowMeans(dosage) / 2
informative <- allele_frequency > 0.05 & allele_frequency < 0.95
dgr_hap1 <- dgr_hap1[informative, , drop = FALSE]
dgr_hap2 <- dgr_hap2[informative, , drop = FALSE]
dosage <- dosage[informative, , drop = FALSE]
marker_id <- marker_id[informative]
n_markers <- length(marker_id)
message("Retained ", n_markers, " polymorphic markers of ",
        length(informative))

## Physical map. Markers sit about one megabase apart within a block, with a
## larger gap between blocks.
position <- integer(0)
chromosome_label <- character(0)
retained_index <- which(informative)
original_chromosome <- rep(seq_len(N_CHR), each = markers_per_chr)
original_block <- rep(
  rep(seq_len(BLOCKS_PER_CHR), each = MARKERS_PER_BLOCK), times = N_CHR
)
original_within <- rep(
  rep(seq_len(MARKERS_PER_BLOCK), times = BLOCKS_PER_CHR), times = N_CHR
)
position_all <- (original_block - 1L) * 25e6 + (original_within - 1L) * 1e6 + 1e6

dgr_map <- data.frame(
  variant_id = marker_id,
  chromosome = original_chromosome[retained_index],
  position_bp = as.numeric(position_all[retained_index]),
  stringsAsFactors = FALSE
)

# -- Genetic values generated from the markers --------------------------------

N_QTL <- 90L
qtl_index <- sort(sample.int(n_markers, N_QTL))
qtl_dosage <- t(dosage[qtl_index, , drop = FALSE])

## Correlated quantitative-trait-locus effects, so the traits share a genuine
## genetic basis rather than being drawn independently.
effect_root <- chol(genetic_correlation)
qtl_effects <- matrix(rnorm(N_QTL * N_TRAITS), nrow = N_QTL) %*% effect_root
raw_genetic <- scale(qtl_dosage, center = TRUE, scale = FALSE) %*% qtl_effects

## Whiten, then recolour to the target covariance, so the realised genetic
## covariance of these candidates equals dgr_G rather than merely resembling it.
raw_centred <- scale(raw_genetic, center = TRUE, scale = FALSE)
whitening <- solve(chol(stats::cov(raw_centred)))
genetic_values <- raw_centred %*% whitening %*% chol(dgr_G)
colnames(genetic_values) <- TRAITS
stopifnot(max(abs(stats::cov(genetic_values) - dgr_G)) < 1e-8)
message("Realised genetic covariance matches dgr_G to within 1e-8")

# -- Adjusted means and genomic estimated breeding values ---------------------

multivariate_normal <- function(n, Sigma) {
  decomposition <- eigen(Sigma, symmetric = TRUE)
  root <- decomposition$vectors %*%
    diag(sqrt(pmax(decomposition$values, 0)), nrow(Sigma))
  matrix(rnorm(n * nrow(Sigma)), n) %*% t(root)
}

residuals_matrix <- multivariate_normal(N_IND, residual_covariance)
colnames(residuals_matrix) <- TRAITS
adjusted_means <- sweep(
  genetic_values + residuals_matrix, 2L, dgr_traits$mean, "+"
)

dgr_candidates <- data.frame(
  GenoID = candidate_id,
  Family = sprintf("F%02d", rep(seq_len(20L), each = N_IND / 20L)),
  round(as.data.frame(adjusted_means), 4L),
  stringsAsFactors = FALSE
)

## Genomic estimated breeding values follow the prediction identity: the
## predictor and the prediction error are orthogonal, so the true genetic value
## is the prediction plus an independent error, and the reliability determines
## how much of the genetic variance the prediction retains.
reliability <- c(GY = 0.35, PHT = 0.65, AD = 0.70, ASI = 0.30,
                 EPP = 0.40, GLS = 0.50)[TRAITS]
reliability_root <- diag(sqrt(reliability))
prediction_covariance <- reliability_root %*% dgr_G %*% reliability_root
dimnames(prediction_covariance) <- list(TRAITS, TRAITS)
error_covariance <- dgr_G - prediction_covariance
error_covariance <- nearest_positive_definite(
  stats::cov2cor(error_covariance), label = "Prediction error correlation"
)
error_sd <- sqrt(pmax(diag(dgr_G - prediction_covariance), 1e-10))
error_covariance <- diag(error_sd) %*% error_covariance %*% diag(error_sd)

gebv <- multivariate_normal(N_IND, prediction_covariance)
colnames(gebv) <- TRAITS
dgr_gebv <- data.frame(
  GenoID = candidate_id,
  round(as.data.frame(gebv), 4L),
  stringsAsFactors = FALSE
)

# -- A previous cycle's selection decision ------------------------------------

## The historical decision is generated from a known weight vector, which is
## attached to the dataset so that the recovery achieved by
## retrospective_weights() can be checked rather than merely asserted. An
## example whose answer cannot be verified demonstrates nothing.
##
## The weights act on the favourable-direction, standardised scale, so a
## positive value always means improvement.
generating_weights <- c(GY = 1.00, PHT = 0.25, AD = 0.60, ASI = 0.45,
                        EPP = 0.30, GLS = 0.55)[TRAITS]
hidden_weights <- generating_weights
direction_sign <- ifelse(dgr_traits$direction == "increase", 1, -1)
names(direction_sign) <- TRAITS
standardised <- scale(adjusted_means)
oriented <- sweep(standardised, 2L, direction_sign, "*")
historical_merit <- as.numeric(oriented %*% hidden_weights)
selected_rows <- order(-historical_merit)[seq_len(40L)]

dgr_history <- data.frame(
  GenoID = candidate_id,
  selected = seq_len(N_IND) %in% selected_rows,
  stringsAsFactors = FALSE
)
attr(dgr_history, "generating_weights") <- generating_weights
attr(dgr_history, "generating_scale") <- paste(
  "Standardised, favourable-direction trait space: each trait was centred and",
  "scaled, then oriented so that larger is better, before the weights were",
  "applied."
)
message("Historical decision: ", sum(dgr_history$selected), " of ", N_IND,
        " candidates selected (", round(100 * 40 / N_IND, 1), "%)")
message("Generating weights attached to dgr_history for verification.")

# -- Save ---------------------------------------------------------------------

usethis::use_data(dgr_traits,     overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_G,          overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_P,          overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_candidates, overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_gebv,       overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_history,    overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_hap1,       overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_hap2,       overwrite = TRUE, compress = "xz")
usethis::use_data(dgr_map,        overwrite = TRUE, compress = "xz")

message("Saved nine datasets to data/")

# -- Sanity checks that the data demonstrate what they are meant to -----------

condition_number <- function(M) {
  values <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  max(abs(values)) / min(abs(values))
}
scale_ratio <- max(genetic_sd) / min(genetic_sd)

message("Trait standard deviations span a factor of ", round(scale_ratio, 0),
        ", so the standardisation diagnostics have something to detect")
message("Condition number of G: ", format(condition_number(dgr_G), digits = 4))
message("Condition number of P: ", format(condition_number(dgr_P), digits = 4))
message("Strongest antagonism: GY with AD, genetic correlation ",
        genetic_correlation["GY", "AD"],
        " while the objective requires GY up and AD down")
message("Data generation complete.")
