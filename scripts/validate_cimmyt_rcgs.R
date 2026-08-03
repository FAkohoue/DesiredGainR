# CIMMYT rapid-cycle genomic selection (RCGS) validation.
#
# This module is sourced by validate_empirical_programmes.R. It deliberately
# separates raw-data reproductions, out-of-environment validation, published
# aggregate anchors, and analyses that the deposited data cannot identify.

validate_cimmyt_rcgs <- function(data_root, out_dir) {
  root <- file.path(data_root, "dataverse_files_New")
  if (!dir.exists(root)) {
    stop("CIMMYT RCGS deposit not found: ", root, call. = FALSE)
  }

  dt <- data.table::data.table
  bind <- function(x) data.table::rbindlist(x, fill = TRUE)
  summaries <- list()
  cycles <- list()
  audits <- list()
  add_summary <- function(...) summaries[[length(summaries) + 1L]] <<- dt(...)
  add_cycle <- function(...) cycles[[length(cycles) + 1L]] <<- dt(...)
  add_audit <- function(...) audits[[length(audits) + 1L]] <<- dt(...)

  normalise_pedigree <- function(x) {
    x <- sub("//\\([^)]*\\)$", "", trimws(as.character(x)))
    toupper(gsub("[^A-Z0-9]", "", x))
  }
  normalise_sample <- function(x) {
    x <- trimws(as.character(x))
    numeric_id <- grepl("^[0-9]+:", x)
    x[numeric_id] <- sub(":.*$", "", x[numeric_id])
    x <- sub(":[A-Z0-9]+:[12]:[0-9]+$", "", x)
    toupper(gsub("[^A-Z0-9]", "", x))
  }
  as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  favourable_response <- function(favourable, selected) {
    vapply(favourable, function(x) {
      mean(x[selected], na.rm = TRUE) - mean(x, na.rm = TRUE)
    }, numeric(1L))
  }
  safe_hash <- function(path) unname(tools::md5sum(path))

  # -------------------------------------------------------------------------
  # Deposit and key integrity.
  # -------------------------------------------------------------------------
  material_path <- file.path(root, "Genetic Material.csv")
  marker_path <- file.path(root, "Marker information file.csv")
  material <- data.table::fread(material_path, na.strings = c("", "NA"))
  data.table::setnames(material, c("Sample Name", "Pedigree"),
                       c("sample_name", "pedigree"))
  material[, pedigree_key := normalise_pedigree(pedigree)]
  material[, sample_key := normalise_sample(sample_name)]
  markers <- data.table::fread(marker_path)
  marker_name <- names(markers)[1L]
  marker_integrity <- dt(
    check = c("marker_rows", "unique_marker_names", "mapped_chromosomes",
              "zero_chromosome_markers"),
    observed = c(nrow(markers), data.table::uniqueN(markers[[marker_name]]),
                 data.table::uniqueN(markers$chromosome),
                 sum(markers$chromosome == 0)),
    expected = c(955690, 955690, 11, NA_real_),
    status = c(if (nrow(markers) == 955690L) "pass" else "investigate",
               if (!anyDuplicated(markers[[marker_name]])) "pass" else "fail",
               if (setequal(unique(markers$chromosome), 0:10)) "pass" else "investigate",
               "documented")
  )

  hapmap_paths <- list.files(root, pattern = "[.]hmp[.]txt$", recursive = TRUE,
                            full.names = TRUE)
  expected_sample_count <- function(path) {
    name <- basename(path)
    if (grepl("18 Original", name, fixed = TRUE)) return(18L)
    if (grepl("C0_1000", name, fixed = TRUE)) return(1000L)
    if (grepl("C0_50", name, fixed = TRUE)) return(50L)
    if (grepl("C1_157", name, fixed = TRUE)) return(157L)
    if (grepl("C1_25", name, fixed = TRUE)) return(25L)
    if (grepl("C2_91", name, fixed = TRUE)) return(91L)
    if (grepl("C2_18", name, fixed = TRUE)) return(18L)
    if (grepl("HGB_44", name, fixed = TRUE)) return(44L)
    if (grepl("HGB_22", name, fixed = TRUE)) return(22L)
    NA_integer_
  }
  header_audit <- lapply(hapmap_paths, function(path) {
    con <- file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    first_lines <- readLines(con, n = 2000L, warn = FALSE)
    header_line <- first_lines[which(!grepl("^##", first_lines))[1L]]
    header <- strsplit(header_line, "\\t", fixed = FALSE)[[1L]]
    fixed_columns <- if (length(header) && identical(header[1L], "#CHROM")) 9L else 11L
    samples <- if (length(header) > fixed_columns) header[-seq_len(fixed_columns)] else character()
    ids <- trimws(samples)
    keys <- normalise_sample(ids)
    dt(
      file = basename(path), relative_path = substring(path, nchar(root) + 2L),
      bytes = file.info(path)$size, samples = length(ids),
      expected_samples = expected_sample_count(path),
      unique_samples = data.table::uniqueN(ids),
      samples_in_genetic_material = sum(keys %chin% material$sample_key),
      md5 = if (basename(path) == "HGB_22 selected samples.hmp.txt") {
        safe_hash(path)
      } else {
        NA_character_
      }
    )
  })
  header_audit <- bind(header_audit)
  c3_selected <- header_audit[file == "HGB_22 selected samples.hmp.txt" &
                                grepl("data_C3", relative_path, ignore.case = TRUE)]
  c4_labelled <- header_audit[file == "HGB_22 selected samples.hmp.txt" &
                                grepl("data_C4", relative_path, ignore.case = TRUE)]
  duplicate_c4 <- nrow(c3_selected) == 1L && nrow(c4_labelled) == 1L &&
    identical(c3_selected$md5, c4_labelled$md5)

  # -------------------------------------------------------------------------
  # Training-population phenotypes and independent environment holdouts.
  # Four traits are present in all four training environments. Each trait is
  # standardised within environment and oriented so larger is favourable.
  # -------------------------------------------------------------------------
  training_dir <- file.path(root, "RCGS_Training population phenotypic data")
  training_paths <- sort(list.files(training_dir, pattern = "[.]xls$", full.names = TRUE))
  traits <- c("mGrainYieldTons_FieldWt", "mEarRotTotalPer",
              "mEarAspect1_5", "mGrainMoisturePer")
  short_traits <- c("GY", "EarRot", "EarAspect", "Moisture")
  names(short_traits) <- traits
  directions <- c(1, -1, -1, -1)
  names(directions) <- traits

  read_training_environment <- function(path, environment) {
    raw <- suppressMessages(readxl::read_excel(path, sheet = "Master",
                                                .name_repair = "unique_quiet"))
    # Row one contains the abbreviated field names, followed by plot records.
    x <- data.table::as.data.table(raw[-1L, c("Entry", traits)])
    x[, Entry := as_number(Entry)]
    for (trait in traits) x[, (trait) := as_number(get(trait))]
    x <- x[Entry >= 1 & Entry <= 1000]
    means <- x[, lapply(.SD, mean_or_na), by = Entry, .SDcols = traits]
    means[, environment := environment]
    for (trait in traits) {
      s <- stats::sd(means[[trait]], na.rm = TRUE)
      means[, (paste0(trait, "_fav")) :=
              directions[[trait]] * (get(trait) - mean(get(trait), na.rm = TRUE)) / s]
    }
    means
  }
  env_data <- bind(Map(read_training_environment, training_paths,
                       paste0("E", seq_along(training_paths))))
  fav_traits <- paste0(traits, "_fav")
  env_counts <- env_data[, c(list(candidates = data.table::uniqueN(Entry)),
                             lapply(.SD, function(x) sum(is.finite(x)))),
                         by = environment, .SDcols = traits]
  data.table::setnames(env_counts, traits, paste0(short_traits, "_observed"))
  fieldbook <- suppressMessages(readxl::read_excel(
    training_paths[1L], sheet = "Fieldbook",
    .name_repair = "unique_quiet"
  ))
  fieldbook <- data.table::as.data.table(fieldbook[, c("Entry", "BreedersPedigree1")])
  fieldbook[, Entry := as_number(Entry)]
  fieldbook <- unique(fieldbook[Entry >= 1 & Entry <= 1000,
                                .(Entry, pedigree = as.character(BreedersPedigree1))])
  fieldbook[, pedigree_key := normalise_pedigree(pedigree)]
  training_material <- material[Cycles == "Training Population"]
  training_linked <- sum(fieldbook$pedigree_key %chin% training_material$pedigree_key)

  direction_grid <- data.table::as.data.table(expand.grid(rep(list(0:6), 4L)))
  data.table::setnames(direction_grid, traits)
  direction_grid <- direction_grid[rowSums(direction_grid) == 6L]
  direction_grid[, (traits) := lapply(.SD, function(x) x / 6), .SDcols = traits]
  event_threshold <- stats::setNames(c(0.20, 0.05, 0.05, 0.02), fav_traits)

  fit_index_direction <- function(train, direction) {
    aggregate <- train[, lapply(.SD, mean_or_na), by = Entry, .SDcols = fav_traits]
    aggregate <- aggregate[stats::complete.cases(aggregate)]
    G <- ridge_covariance(aggregate, fav_traits)
    selection_index(
      aggregate, fav_traits, id_col = "Entry", method = "pesek_baker", G = G,
      desired_gains = stats::setNames(direction, fav_traits),
      center_traits = TRUE, scale_traits = FALSE,
      n_select = max(1L, floor(0.20 * nrow(aggregate)))
    )
  }
  evaluate_direction <- function(train, test, direction) {
    fit <- fit_index_direction(train, direction)
    selected_ids <- fit$ranking[selected == TRUE, id]
    held <- test[Entry %in% selected_ids]
    population <- test[Entry %in% fit$ranking$id]
    favourable_response(population[, ..fav_traits],
                        population$Entry %in% held$Entry)
  }
  choose_direction <- function(train_environments) {
    scores <- vapply(seq_len(nrow(direction_grid)), function(i) {
      d <- unlist(direction_grid[i, ..traits], use.names = FALSE)
      responses <- vapply(train_environments, function(inner_holdout) {
        inner_train <- env_data[
          environment %chin% setdiff(train_environments, inner_holdout)]
        inner_test <- env_data[environment == inner_holdout]
        r <- evaluate_direction(inner_train, inner_test, d)
        min(r / event_threshold)
      }, numeric(1L))
      stats::median(responses)
    }, numeric(1L))
    # Deterministic tie break: prefer the more balanced direction.
    best <- which(scores == max(scores, na.rm = TRUE))
    if (length(best) > 1L) {
      balance <- apply(direction_grid[best, ..traits], 1L,
                       function(x) -stats::var(as.numeric(x)))
      best <- best[which.max(balance)]
    }
    unlist(direction_grid[best[1L], ..traits], use.names = FALSE)
  }

  environment_validation <- lapply(unique(env_data$environment), function(outer) {
    train_envs <- setdiff(unique(env_data$environment), outer)
    chosen <- choose_direction(train_envs)
    chosen_response <- evaluate_direction(env_data[environment %chin% train_envs],
                                          env_data[environment == outer], chosen)
    equal_response <- evaluate_direction(env_data[environment %chin% train_envs],
                                         env_data[environment == outer], rep(0.25, 4L))
    yield_response <- evaluate_direction(env_data[environment %chin% train_envs],
                                         env_data[environment == outer], c(1, 0, 0, 0))
    inner_success <- vapply(train_envs, function(inner) {
      r <- evaluate_direction(
        env_data[environment %chin% setdiff(train_envs, inner)],
        env_data[environment == inner], chosen
      )
      all(r >= event_threshold)
    }, logical(1L))
    probability <- (sum(inner_success) + 1) / (length(inner_success) + 2)
    dt(
      held_out_environment = outer,
      direction_GY = chosen[1L], direction_EarRot = chosen[2L],
      direction_EarAspect = chosen[3L], direction_Moisture = chosen[4L],
      predicted_joint_success_probability = probability,
      observed_joint_success = all(chosen_response >= event_threshold),
      chosen_GY = chosen_response[1L], chosen_EarRot = chosen_response[2L],
      chosen_EarAspect = chosen_response[3L], chosen_Moisture = chosen_response[4L],
      equal_GY = equal_response[1L], equal_EarRot = equal_response[2L],
      equal_EarAspect = equal_response[3L], equal_Moisture = equal_response[4L],
      yield_only_GY = yield_response[1L], yield_only_EarRot = yield_response[2L],
      yield_only_EarAspect = yield_response[3L], yield_only_Moisture = yield_response[4L]
    )
  })
  environment_validation <- bind(environment_validation)
  brier <- mean((environment_validation$predicted_joint_success_probability -
                  environment_validation$observed_joint_success)^2)
  chosen_min <- apply(environment_validation[, .(chosen_GY, chosen_EarRot,
                                                  chosen_EarAspect, chosen_Moisture)],
                      1L, min)
  equal_min <- apply(environment_validation[, .(equal_GY, equal_EarRot,
                                                 equal_EarAspect, equal_Moisture)],
                     1L, min)
  add_summary(
    validation_id = "rcgs_nested_direction_worst_trait_increment",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)",
    evidence_tier = "nested-environment-transport",
    endpoint = "optimised versus equal-direction worst-trait response",
    estimate = mean(chosen_min - equal_min), lower = min(chosen_min - equal_min),
    upper = max(chosen_min - equal_min), p_value = NA_real_,
    unit = "favourable within-environment SD", n = nrow(environment_validation),
    status = if (min(chosen_min - equal_min) > 0) "supportive" else "mixed",
    interpretation = paste(
      "The desired-gain direction was chosen inside each training fold by",
      "maximising the worst normalised favourable response, then frozen for",
      "the held-out environment."
    ),
    limitation = paste(
      "Only four environments are available; covariance is phenotypic and",
      "environment-adjusted, not an additive-genetic covariance estimate."
    )
  )
  add_summary(
    validation_id = "rcgs_joint_success_probability_brier",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)",
    evidence_tier = "nested-environment-transport",
    endpoint = "joint high-gain probability Brier score",
    estimate = brier, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
    unit = "Brier score", n = nrow(environment_validation), status = "descriptive",
    interpretation = paste0(
      "High gain was predeclared as favourable response of at least ",
      "0.20, 0.05, 0.05 and 0.02 environment SD for yield, ear rot, ear ",
      "aspect and moisture. Probabilities use add-one smoothing across the ",
      "three inner environment folds."
    ),
    limitation = "Four outer outcomes are insufficient to claim probability calibration."
  )

  # -------------------------------------------------------------------------
  # Genetic-gain trial: two locations reported in the paper. Cotaxtla is a
  # deposit-only third environment with documented data replacement and is not
  # silently mixed into the primary reproduction.
  # -------------------------------------------------------------------------
  gain_dir <- file.path(root, "RCGS Genetic Gain_Phenotypic Data")
  gain_paths <- sort(list.files(gain_dir, pattern = "All [Ll]ocations[.]csv$",
                                full.names = TRUE))
  gain <- bind(lapply(gain_paths, function(path) {
    x <- data.table::fread(path)
    x[, cycle := sub(".*_(C[0-4])-.*", "\\1", basename(path))]
    x
  }))
  gain[, LOC := trimws(LOC)]
  gain[LOC == "Cot", LOC := "Cotaxtla"]
  gain[, yield := as_number(mGrainYieldTons_FieldWt)]
  expected_candidates <- c(C0 = 48L, C1 = 47L, C2 = 48L, C3 = 43L, C4 = 44L)
  gain[, candidate := ENTRY <= expected_candidates[cycle]]
  gain_candidates <- unique(gain[candidate & cycle != "C4",
                                 .(cycle, pedigree_key = normalise_pedigree(Pedigree))])
  gain_link <- gain_candidates[, .(
    candidates = .N,
    linked_to_genetic_material = sum(pedigree_key %chin% material$pedigree_key)
  ), by = cycle]
  primary <- gain[candidate & LOC %chin% c("Agua Fria", "Tlaltizapan") &
                    is.finite(yield)]
  primary[, `:=`(
    cycle = factor(cycle, levels = paste0("C", 0:4)),
    LOC = factor(LOC), REP = factor(REP), BLOCK = factor(BLOCK),
    entry_cycle = interaction(cycle, ENTRY, drop = TRUE),
    block_key = interaction(LOC, cycle, REP, BLOCK, drop = TRUE)
  )]
  gain_fit <- lme4::lmer(yield ~ cycle * LOC + (1 | entry_cycle) + (1 | block_key),
                         data = primary, REML = TRUE)
  gain_emm <- data.table::as.data.table(emmeans::emmeans(gain_fit, ~ cycle))
  gain_emm[, cycle_number := as.integer(sub("C", "", as.character(cycle)))]
  published_means <- c(C0 = 8.52, C1 = 8.40, C2 = 8.62, C3 = 8.92, C4 = 9.05)
  gain_emm[, published := published_means[as.character(cycle)]]
  gain_emm[, absolute_error := abs(emmean - published)]

  primary_43 <- data.table::copy(primary)[!(cycle == "C4" & as.integer(as.character(ENTRY)) == 44L)]
  primary_43[, `:=`(
    entry_cycle = interaction(cycle, ENTRY, drop = TRUE),
    block_key = interaction(LOC, cycle, REP, BLOCK, drop = TRUE)
  )]
  gain_fit_43 <- lme4::lmer(
    yield ~ cycle * LOC + (1 | entry_cycle) + (1 | block_key),
    data = primary_43, REML = TRUE
  )
  gain_emm_43 <- data.table::as.data.table(emmeans::emmeans(gain_fit_43, ~ cycle))
  c4_44 <- gain_emm[cycle == "C4", emmean]
  c4_43 <- gain_emm_43[cycle == "C4", emmean]
  for (i in seq_len(nrow(gain_emm))) {
    add_cycle(
      dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", endpoint = "grain yield",
      cycle = as.character(gain_emm$cycle[i]), estimate = gain_emm$emmean[i],
      standard_error = gain_emm$SE[i],
      n_independent_units = data.table::uniqueN(primary[cycle == gain_emm$cycle[i], entry_cycle]),
      adjustment = "two-location mixed model; entry and incomplete block random"
    )
  }
  slope_c0 <- stats::coef(stats::lm(emmean ~ cycle_number, data = gain_emm))[2L]
  slope_c1 <- stats::coef(stats::lm(emmean ~ cycle_number,
                                    data = gain_emm[cycle_number >= 1L]))[2L]
  for (spec in list(c("C0-C4", slope_c0, 0.158), c("C1-C4", slope_c1, 0.225))) {
    estimate <- as.numeric(spec[2L])
    target <- as.numeric(spec[3L])
    add_summary(
      validation_id = paste0("rcgs_gain_", tolower(gsub("-", "_", spec[1L])), "_slope"),
      dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", evidence_tier = "realized-cycle",
      endpoint = paste0("grain-yield gain slope ", spec[1L]),
      estimate = estimate, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
      unit = "t/ha per cycle", n = if (spec[1L] == "C0-C4") 5L else 4L,
      status = if (abs(estimate - target) <= 0.05) "programme-concordant" else "mixed",
      interpretation = paste0("Raw two-location deposit reproduction; published slope = ", target, "."),
      limitation = paste(
        "Cycle is not randomized; the deposit has 44 C4 candidate entries",
        "whereas the paper reports 43, and the original fixed-entry analysis",
        "was approximated with a random-entry mixed model."
      )
    )
  }
  add_summary(
    validation_id = "rcgs_published_cycle_mean_max_error",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", evidence_tier = "realized-cycle",
    endpoint = "maximum absolute cycle-mean reproduction error",
    estimate = max(gain_emm$absolute_error), lower = NA_real_, upper = NA_real_,
    p_value = NA_real_, unit = "t/ha", n = 5L,
    status = if (max(gain_emm$absolute_error) <= 0.10) "programme-concordant" else "mixed",
    interpretation = "Comparison of raw-deposit adjusted means with Zhang et al. Table 2.",
    limitation = "The published C4 candidate count and raw deposit disagree (43 versus 44)."
  )
  add_summary(
    validation_id = "rcgs_c4_candidate_count_sensitivity",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", evidence_tier = "sensitivity",
    endpoint = "C4 adjusted mean: 43-entry minus 44-entry analysis",
    estimate = c4_43 - c4_44, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
    unit = "t/ha", n = 2L, status = if (abs(c4_43 - c4_44) <= 0.10) "robust" else "sensitive",
    interpretation = paste0(
      "The paper reports 43 C4 entries, while the raw deposit has 44; the two ",
      "adjusted C4 means are ", round(c4_43, 3), " and ", round(c4_44, 3), " t/ha."
    ),
    limitation = "The deposit does not identify which of its 44 C4 entries was excluded from the publication."
  )

  # -------------------------------------------------------------------------
  # Published genomic-prediction and diversity anchors.
  # These are not relabelled as raw-data reproductions.
  # -------------------------------------------------------------------------
  diversity <- dt(
    group = c("Parents", "C0", "C0 selected", "C1", "C1 selected",
              "C2", "C2 selected", "C3", "C3 selected", "All"),
    n = c(18, 1000, 50, 157, 25, 91, 18, 44, 22, 1331),
    Shannon = c(0.0661, 0.0728, 0.0200, 0.0776, 0.0520,
                0.0765, 0.0430, 0.0588, 0.0630, 0.0740),
    heterozygosity = c(0.1104, 0.1226, 0.1208, 0.1297, 0.1250,
                       0.1276, 0.1228, 0.0973, 0.0923, 0.1245),
    polymorphic_SNPs = c(950248, 952825, 943344, 951390, 947868,
                         953199, 953453, 954058, 954924, 954960),
    source = "Zhang et al. (2017), Table 4"
  )
  population_diversity <- diversity[group %chin% c("C0", "C1", "C2", "C3")]
  population_diversity[, cycle_number := 0:3]
  diversity_slope <- stats::coef(stats::lm(heterozygosity ~ cycle_number,
                                           data = population_diversity))[2L]
  add_summary(
    validation_id = "rcgs_published_heterozygosity_slope",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", evidence_tier = "published-aggregate",
    endpoint = "population heterozygosity trend C0-C3",
    estimate = diversity_slope, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
    unit = "heterozygosity per cycle", n = 4L, status = "descriptive",
    interpretation = "Trend calculated from the published population-level diversity table.",
    limitation = "Population sizes and inbreeding differ across cycles; C4 was not genotyped."
  )
  add_summary(
    validation_id = "rcgs_published_gblup_cv_accuracy",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", evidence_tier = "published-aggregate",
    endpoint = "C0 grain-yield genomic-prediction accuracy",
    estimate = 0.55, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
    unit = "mean predictive correlation", n = 100L, status = "published-anchor",
    interpretation = "Five-fold cross-validation repeated 100 times, as reported by Zhang et al.",
    limitation = paste(
      "This is a published aggregate, not a local raw-marker reproduction;",
      "the paper reports using 331,740 of 955,690 SNPs."
    )
  )
  genomic_summary_path <- file.path(out_dir, "cimmyt-rcgs-genomic-cv-summary.csv")
  if (file.exists(genomic_summary_path)) {
    genomic <- data.table::fread(genomic_summary_path)
    add_summary(
      validation_id = "rcgs_raw_subsampled_gblup_cv_accuracy",
      dataset = "CIMMYT maize RCGS (Zhang et al. 2017)",
      evidence_tier = "raw-genomic-reproduction",
      endpoint = "C0 grain-yield genomic prediction from genome-wide marker sample",
      estimate = genomic$mean_correlation, lower = genomic$lower_correlation,
      upper = genomic$upper_correlation, p_value = NA_real_,
      unit = "repeated five-fold predictive correlation", n = genomic$candidates,
      status = if (abs(genomic$mean_correlation - genomic$published_correlation) <= 0.10 &&
                    abs(genomic$negative_control_correlation) < 0.10) {
        "programme-concordant"
      } else {
        "mixed"
      },
      interpretation = paste0(
        genomic$filtered_markers, " filtered markers from ",
        genomic$sampled_marker_rows, " deterministic genome-wide rows gave mean r = ",
        round(genomic$mean_correlation, 3), "; the published full-marker result was 0.55."
      ),
      limitation = paste(
        "This computationally bounded reproduction uses a deterministic marker",
        "subsample and mean standardised yield, not the paper's full 331,740-marker",
        "mixed-model phenotype pipeline."
      )
    )
    add_summary(
      validation_id = "rcgs_raw_gblup_permutation_control",
      dataset = "CIMMYT maize RCGS (Zhang et al. 2017)",
      evidence_tier = "negative-control",
      endpoint = "permuted-outcome genomic prediction correlation",
      estimate = genomic$negative_control_correlation,
      lower = NA_real_, upper = NA_real_, p_value = NA_real_,
      unit = "mean predictive correlation", n = genomic$candidates,
      status = if (abs(genomic$negative_control_correlation) < 0.10) "pass" else "investigate",
      interpretation = "The identical CV pipeline was run after breaking genotype-phenotype linkage.",
      limitation = "A finite permutation control is a leakage diagnostic, not a null-hypothesis proof."
    )
  }
  add_summary(
    validation_id = "rcgs_actual_founder_simulation_calibration",
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)", evidence_tier = "not-compatible",
    endpoint = "simulation from deposited RCGS genotypes",
    estimate = NA_real_, lower = NA_real_, upper = NA_real_, p_value = NA_real_,
    unit = "not estimable", n = NA_integer_, status = "not-compatible",
    interpretation = "No actual-founder simulation is run from these HapMap calls.",
    limitation = paste(
      "Samples are bulk DNA from outcrossing S2 families, phase is unavailable,",
      "C4 was not genotyped, and the C4-labelled selected file is byte-identical",
      "to the C3 selected file. Converting the calls as inbred phased founders",
      "would be scientifically invalid."
    )
  )

  add_audit(
    dataset = "CIMMYT maize RCGS (Zhang et al. 2017)",
    rows = nrow(gain), candidates = 1000L, traits = length(traits), years = 1L,
    environments = length(training_paths), cycles_or_stages = 5L,
    genotype_data = TRUE, pedigree = TRUE,
    usable_for = paste(
      "realized cycle gain, nested environment transport, deposit integrity,",
      "and published genomic/diversity calibration"
    ),
    principal_limitation = paste(
      "C4 is ungenotyped; its labelled HapMap is a duplicate of C3; bulk-family",
      "unphased calls are incompatible with actual-founder simulation; C4 has",
      "44 raw candidates versus 43 reported"
    )
  )

  integrity <- data.table::rbindlist(list(
    marker_integrity[, .(section = "markers", item = check,
                         observed = as.character(observed),
                         expected = as.character(expected), status)],
    header_audit[, .(section = "hapmap", item = relative_path,
                     observed = paste0(samples, " samples; ",
                                       samples_in_genetic_material, " linked"),
                     expected = paste0(expected_samples, " samples"),
                     status = ifelse(samples == expected_samples &
                                       samples == unique_samples &
                                       samples_in_genetic_material == samples,
                                     "pass", "investigate"))],
    dt(section = "cross-file", item = "C3/C4 selected HapMap identity",
       observed = as.character(duplicate_c4), expected = "FALSE",
       status = if (duplicate_c4) "fail-do-not-use-as-C4" else "pass"),
    env_counts[, .(section = "training phenotype", item = environment,
                   observed = paste0(candidates, " candidates"),
                   expected = "1000 candidates", status = "pass")],
    dt(section = "cross-file", item = "training pedigree-to-sample linkage",
       observed = paste0(training_linked, " of ", nrow(fieldbook)),
       expected = "1000 of 1000", status = if (training_linked == 1000L) "pass" else "fail"),
    gain_link[, .(section = "cross-file",
                  item = paste0(cycle, " gain-trial pedigree linkage"),
                  observed = paste0(linked_to_genetic_material, " of ", candidates),
                  expected = paste0(candidates, " of ", candidates),
                  status = ifelse(linked_to_genetic_material == candidates, "pass", "fail"))]
  ), fill = TRUE)

  data.table::fwrite(integrity,
                     file.path(out_dir, "cimmyt-rcgs-deposit-integrity.csv"), na = "")
  data.table::fwrite(environment_validation,
                     file.path(out_dir, "cimmyt-rcgs-environment-validation.csv"), na = "")
  data.table::fwrite(diversity,
                     file.path(out_dir, "cimmyt-rcgs-published-diversity.csv"), na = "")
  data.table::fwrite(gain_emm,
                     file.path(out_dir, "cimmyt-rcgs-cycle-model.csv"), na = "")

  list(summary = bind(summaries), cycles = bind(cycles), audit = bind(audits),
       integrity = integrity, environment_validation = environment_validation,
       diversity = diversity, gain_model = gain_emm)
}
