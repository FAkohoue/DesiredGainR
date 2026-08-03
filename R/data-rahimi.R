# Published reference values from Rahimi and Debnath (2023), used for the
# frozen external reproduction.
#
# These numbers were produced by the authors' own SAS PROC IML and R code on
# their maize data, published in Scientific Reports and independent of this
# package. They are transcribed from the article's Tables 2 and 3 and are
# frozen: they must not be recomputed, adjusted or "corrected" here, because
# their whole value is that they came from somewhere else.

#' Published selection-index results from Rahimi and Debnath (2023)
#'
#' Reference values transcribed from Tables 2 and 3 of Rahimi and Debnath
#' (2023), for seven traits measured on 28 maize inbred lines evaluated in a
#' randomised complete block design with three replications.
#'
#' @format A list with components:
#' \describe{
#'   \item{`traits`}{The seven trait names, in the article's order.}
#'   \item{`desired_gains`}{The \eqn{d} vector used for the Pesek-Baker index,
#'     which the article states is \eqn{\sqrt{\mathrm{diag}(\mathbf{G})}}, the
#'     genetic standard deviations (Table 2, final column).}
#'   \item{`genetic_variances`}{The implied diagonal of \eqn{\mathbf{G}},
#'     obtained as `desired_gains^2`.}
#'   \item{`economic_weights_method1`}{Unit economic weights (Table 2).}
#'   \item{`pesek_baker_gain`}{Expected genetic advance per trait for the
#'     Pesek-Baker index (Table 3).}
#'   \item{`pesek_baker_criteria`}{The article's reported \eqn{R_{HI}},
#'     \eqn{\Delta H}, RE and \eqn{CV_I} for the Pesek-Baker index.}
#'   \item{`optimum_method1_gain`}{Expected genetic advance for the optimum
#'     (Smith-Hazel) index under Method 1 weights (Table 3).}
#'   \item{`optimum_method1_criteria`}{Its reported criteria.}
#' }
#'
#' @details
#' This object freezes the article's printed table values. It is a transcription
#' check, not a full data reanalysis; see the reproduction vignette for that
#' distinction and for the archived-data accession.
#'
#' See `vignette("DesiredGainR-reproduction")` for the reproduction and its
#' stated tolerances.
#'
#' @source
#' Rahimi, M. and Debnath, S. (2023) Estimating optimum and base selection
#' indices in plant and animal breeding programs by development new and simple
#' SAS and R codes.
#' *Scientific Reports* 13, 18977. \doi{10.1038/s41598-023-46368-6}
#'
#' @examples
#' reference <- rahimi_debnath_2023()
#' # The Pesek-Baker property: expected response is exactly proportional to the
#' # desired gains. The published numbers satisfy it to four significant
#' # figures.
#' round(reference$pesek_baker_gain / reference$desired_gains, 5)
#'
#' @export
rahimi_debnath_2023 <- function() {
  traits <- c(
    "plant_height", "number_of_grain", "number_of_row", "row_length",
    "leaf_length", "hundred_grain_weight", "yield"
  )

  # Table 2, final column: "The d for Pesek and Baker index". The article
  # states d = sqrt(vecdiag(G)) in SAS and d = sqrt(diag(G1)) in R.
  desired_gains <- stats::setNames(
    c(27.502, 4.554, 36.096, 0.224, 2.484, 45.046, 2.846), traits
  )

  # Table 3, Pesek-Baker row: the expected genetic advance for each trait.
  pesek_baker_gain <- stats::setNames(
    c(9.2366, 1.5295, 12.1228, 0.0751, 0.8342, 15.1287, 0.9559), traits
  )

  # Table 3, optimum index under Method 1 weights.
  optimum_method1_gain <- stats::setNames(
    c(0.2115, -0.3692, 54.4801, 0.0238, -1.0277, 69.8089, 2.6488), traits
  )

  list(
    traits = traits,
    desired_gains = desired_gains,
    genetic_variances = desired_gains^2,
    economic_weights_method1 = stats::setNames(rep(1, length(traits)), traits),
    pesek_baker_gain = pesek_baker_gain,
    pesek_baker_criteria = c(
      R_HI = 0.0018, delta_H = 9.2231, RE = 0.1986, CV_I = 3.1535
    ),
    optimum_method1_gain = optimum_method1_gain,
    optimum_method1_criteria = c(
      R_HI = 0.9887, delta_H = 125.7762, RE = 0.5504, CV_I = 13.6215
    ),
    citation = paste(
      "Rahimi, M. and Debnath, S. (2023) Estimating optimum and base selection",
      "indices in plant and animal breeding programs by development new and",
      "simple SAS and R codes.",
      "Scientific Reports 13, 18977. doi:10.1038/s41598-023-46368-6"
    ),
    note = paste(
      "Frozen transcription from the published tables. Do not recompute or",
      "adjust these values: their purpose is to be an external check, which",
      "they stop being the moment they are derived from this package."
    )
  )
}
