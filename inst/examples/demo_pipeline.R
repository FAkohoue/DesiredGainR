library(DesiredGainR)
library(data.table)

ext <- system.file("extdata", package = "DesiredGainR")
values <- fread(file.path(ext, "example_pheno.csv"))
gebv <- fread(file.path(ext, "example_gebv.csv"))
traits <- c("YLD", "MY", "MI", "BL", "NBL", "VHB")
dg <- c(YLD = 1.5, MY = 0.5, MI = 0.5, BL = 1, NBL = 1, VHB = 1)
G <- stats::cov(as.matrix(values[, ..traits]))

dgsi <- run_dgsi(
  init_data = values[, .(GenoID, Family)],
  cand_data = values,
  trait_cols = traits,
  dg = dg,
  G = G,
  lower_is_better = c("BL", "NBL", "VHB"),
  n_select = 10,
  n_iter = 200,
  n_rep = 5,
  seed = 42
)

print(dgsi$replicate_diagnostics)
print(head(dgsi$ranked_geno))

w <- c(YLD = 1.5, MY = 0.5, MI = 0.5, BL = 1, NBL = 1, VHB = 1)
W <- diag(c(
  YLD = 0.05, MY = 0.02, MI = 0.02,
  BL = -0.03, NBL = -0.03, VHB = -0.03
))
dimnames(W) <- list(traits, traits)

qgsi <- run_qgsi(
  init_data = values[, .(GenoID, Family)],
  gebv_data = gebv,
  trait_cols = traits,
  linear_weights = w,
  W = W,
  lower_is_better = c("BL", "NBL", "VHB"),
  n_select = 10
)

print(qgsi)
print(qgsi$theoretical_parameters)
print(qgsi$observed_selection_differential)
