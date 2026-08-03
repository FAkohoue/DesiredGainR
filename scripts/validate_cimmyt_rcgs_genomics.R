# Marker-subsampled raw-genomic validation for the CIMMYT RCGS deposit.
#
# This is intentionally a separate, slower script. It streams the 1.9 GB C0
# HapMap instead of loading it into memory, constructs a VanRaden relationship
# matrix from a deterministic genome-wide marker sample, and performs repeated
# five-fold kernel-GBLUP validation. Run from the package root:
#
#   Rscript scripts/validate_cimmyt_rcgs_genomics.R

required <- c("data.table", "readxl")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) {
  stop("Install validation dependencies: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

data_root <- Sys.getenv("DESIREDGAINR_PUBLIC_DATA", unset = "Public data")
root <- file.path(data_root, "dataverse_files_New")
out_dir <- file.path("inst", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
target_markers <- as.integer(Sys.getenv("DESIREDGAINR_RCGS_MARKERS", unset = "6000"))
repeats <- as.integer(Sys.getenv("DESIREDGAINR_RCGS_CV_REPEATS", unset = "10"))
seed <- 20260802L

normalise_sample <- function(x) {
  x <- trimws(as.character(x))
  numeric_id <- grepl("^[0-9]+:", x)
  x[numeric_id] <- sub(":.*$", "", x[numeric_id])
  x <- sub(":[A-Z0-9]+:[12]:[0-9]+$", "", x)
  toupper(gsub("[^A-Z0-9]", "", x))
}
normalise_pedigree <- function(x) {
  x <- sub("//\\([^)]*\\)$", "", trimws(as.character(x)))
  toupper(gsub("[^A-Z0-9]", "", x))
}
as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

iupac <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")
decode_marker <- function(parts) {
  alleles <- strsplit(toupper(parts[2L]), "[/]", fixed = FALSE)[[1L]]
  if (length(alleles) != 2L || any(nchar(alleles) != 1L)) {
    return(rep(NA_real_, length(parts) - 11L))
  }
  a1 <- alleles[1L]
  a2 <- alleles[2L]
  calls <- toupper(parts[-seq_len(11L)])
  out <- rep(NA_real_, length(calls))
  out[calls %in% c(a1, paste0(a1, a1))] <- 0
  out[calls %in% c(a2, paste0(a2, a2))] <- 2
  pair <- paste(sort(c(a1, a2)), collapse = "")
  hetero <- c(paste0(a1, a2), paste0(a2, a1), unname(iupac[pair]))
  out[calls %in% hetero] <- 1
  out
}

hapmap <- list.files(
  file.path(root, "RCGS Genetic gains_Genotypic data", "Genotypic data_C0"),
  pattern = "C0_1000 samples[.]hmp[.]txt$", full.names = TRUE
)
if (length(hapmap) != 1L) stop("C0 1000-sample HapMap not found.", call. = FALSE)
marker_info <- data.table::fread(file.path(root, "Marker information file.csv"))
n_deposit_markers <- nrow(marker_info)
step <- max(1L, floor(n_deposit_markers / target_markers))

con <- file(hapmap, open = "rt")
header <- strsplit(readLines(con, n = 1L, warn = FALSE), "\\t")[[1L]]
sample_ids <- header[-seq_len(11L)]
sample_keys <- normalise_sample(sample_ids)
capacity <- ceiling(n_deposit_markers / step) + 2L
dosage <- matrix(NA_real_, nrow = capacity, ncol = length(sample_ids))
marker_names <- character(capacity)
row_index <- 0L
kept <- 0L
repeat {
  lines <- readLines(con, n = 5000L, warn = FALSE)
  if (!length(lines)) break
  line_numbers <- row_index + seq_along(lines)
  take <- which((line_numbers - 1L) %% step == 0L)
  if (length(take)) {
    parsed <- strsplit(lines[take], "\\t", fixed = FALSE)
    for (j in seq_along(parsed)) {
      parts <- parsed[[j]]
      if (length(parts) != length(sample_ids) + 11L) next
      kept <- kept + 1L
      marker_names[kept] <- parts[1L]
      dosage[kept, ] <- decode_marker(parts)
    }
  }
  row_index <- row_index + length(lines)
}
close(con)
dosage <- dosage[seq_len(kept), , drop = FALSE]
marker_names <- marker_names[seq_len(kept)]

# Candidate phenotype: mean of four within-environment standardised grain-yield
# records. This removes environment scale without borrowing test-fold outcomes
# during genomic CV; all CV splits occur between candidate families.
training_dir <- file.path(root, "RCGS_Training population phenotypic data")
training_paths <- sort(list.files(training_dir, pattern = "[.]xls$", full.names = TRUE))
env_yield <- data.table::rbindlist(lapply(seq_along(training_paths), function(i) {
  raw <- suppressMessages(readxl::read_excel(training_paths[i], sheet = "Master",
                                              .name_repair = "unique_quiet"))
  x <- data.table::as.data.table(raw[-1L, c("Entry", "mGrainYieldTons_FieldWt")])
  x[, Entry := as_number(Entry)]
  x[, yield := as_number(mGrainYieldTons_FieldWt)]
  x <- x[Entry >= 1 & Entry <= 1000 & is.finite(yield)]
  x <- x[, .(yield = mean(yield)), by = Entry]
  x[, y_z := as.numeric(scale(yield))]
  x[, environment := i]
  x[, .(Entry, environment, y_z)]
}))
phenotype <- env_yield[, .(y = mean(y_z)), by = Entry]

fieldbook <- suppressMessages(readxl::read_excel(training_paths[1L], sheet = "Fieldbook",
                                                  .name_repair = "unique_quiet"))
fieldbook <- data.table::as.data.table(fieldbook[, c("Entry", "BreedersPedigree1")])
fieldbook[, Entry := as_number(Entry)]
fieldbook <- unique(fieldbook[Entry >= 1 & Entry <= 1000,
                              .(Entry, pedigree_key = normalise_pedigree(BreedersPedigree1))])
material <- data.table::fread(file.path(root, "Genetic Material.csv"),
                              na.strings = c("", "NA"))
material <- material[Cycles == "Training Population"]
material[, pedigree_key := normalise_pedigree(Pedigree)]
material[, sample_key := normalise_sample(`Sample Name`)]
key <- merge(fieldbook, material[, .(pedigree_key, sample_key)], by = "pedigree_key")
phenotype <- merge(phenotype, key[, .(Entry, sample_key)], by = "Entry")
y <- phenotype$y[match(sample_keys, phenotype$sample_key)]

keep_samples <- is.finite(y) & !duplicated(sample_keys)
y <- y[keep_samples]
sample_keys <- sample_keys[keep_samples]
M <- t(dosage[, keep_samples, drop = FALSE])
call_rate <- colMeans(is.finite(M))
p <- colMeans(M, na.rm = TRUE) / 2
maf <- pmin(p, 1 - p)
keep_markers <- is.finite(p) & call_rate >= 0.90 & maf > 0.05
M <- M[, keep_markers, drop = FALSE]
p <- p[keep_markers]
Z <- sweep(M, 2L, 2 * p, "-")
Z[!is.finite(Z)] <- 0
denominator <- sum(2 * p * (1 - p))
K <- tcrossprod(Z) / denominator
rownames(K) <- colnames(K) <- sample_keys

gcv_lambda <- function(K_train, y_train) {
  eig <- eigen(K_train, symmetric = TRUE)
  d <- pmax(eig$values, 0)
  centred <- y_train - mean(y_train)
  projection <- crossprod(eig$vectors, centred)
  grid <- 10^seq(-3, 3, length.out = 49L)
  score <- vapply(grid, function(lambda) {
    shrink <- d / (d + lambda)
    fitted <- eig$vectors %*% (shrink * projection)
    denom <- 1 - mean(shrink)
    mean((centred - fitted)^2) / denom^2
  }, numeric(1L))
  grid[which.min(score)]
}

run_cv <- function(response, repeats, seed) {
  set.seed(seed)
  out <- data.table::rbindlist(lapply(seq_len(repeats), function(rep) {
    fold <- sample(rep(seq_len(5L), length.out = length(response)))
    prediction <- rep(NA_real_, length(response))
    lambdas <- numeric(5L)
    for (k in seq_len(5L)) {
      test <- which(fold == k)
      train <- which(fold != k)
      lambda <- gcv_lambda(K[train, train, drop = FALSE], response[train])
      alpha <- solve(K[train, train, drop = FALSE] + diag(lambda, length(train)),
                     response[train] - mean(response[train]))
      prediction[test] <- mean(response[train]) +
        K[test, train, drop = FALSE] %*% alpha
      lambdas[k] <- lambda
    }
    data.table::data.table(
      iteration = rep, correlation = stats::cor(prediction, response),
      median_lambda = stats::median(lambdas)
    )
  }))
  out
}

observed <- run_cv(y, repeats, seed)
set.seed(seed + 1L)
permuted <- run_cv(sample(y), max(5L, floor(repeats / 2L)), seed + 2L)
detail <- data.table::rbindlist(list(
  observed[, scenario := "observed"],
  permuted[, scenario := "permuted-negative-control"]
), use.names = TRUE, fill = TRUE)
summary <- data.table::data.table(
  sampled_marker_rows = kept,
  filtered_markers = ncol(M), candidates = length(y),
  repeats = repeats, folds = 5L,
  mean_correlation = mean(observed$correlation),
  sd_correlation = stats::sd(observed$correlation),
  lower_correlation = unname(stats::quantile(observed$correlation, 0.025, type = 8)),
  upper_correlation = unname(stats::quantile(observed$correlation, 0.975, type = 8)),
  published_correlation = 0.55,
  negative_control_correlation = mean(permuted$correlation),
  marker_filter = "C0 call rate >= 0.90 and MAF > 0.05",
  relationship = "VanRaden-I from deterministic genome-wide marker sample",
  validation = "repeated five-fold family CV; GCV-tuned kernel ridge/GBLUP lambda"
)
data.table::fwrite(detail, file.path(out_dir, "cimmyt-rcgs-genomic-cv-detail.csv"), na = "")
data.table::fwrite(summary, file.path(out_dir, "cimmyt-rcgs-genomic-cv-summary.csv"), na = "")
message("CIMMYT RCGS raw-genomic validation written to ", out_dir)
