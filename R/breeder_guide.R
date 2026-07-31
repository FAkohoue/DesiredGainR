# ==============================================================================
# breeder_guide.R
#
# The PDF and HTML editions of the DesiredGainR Breeder's Guide are standalone
# companions to the *Defining a breeding objective* vignette: the six
# objective-setting tools (implied_economic_weights(), implied_desired_gains(),
# gain_feasibility(), retrospective_weights(), weight_sensitivity() and
# effective_weights()) together with the choice among index families,
# explained in plain language for readers who will approve a selection
# decision but do not themselves run R.
#
# Both installed editions are generated from the single R Markdown source
# inst/guide/DesiredGainR_Breeder_Guide.Rmd; open_desiredgain_guide() locates
# either edition without requiring users to know the system.file() path.
# ==============================================================================

#' Locate or Open the DesiredGainR Breeder's Guide
#'
#' The PDF and HTML editions of the DesiredGainR Breeder's Guide are standalone
#' companions to the \emph{Defining a breeding objective} vignette. They cover
#' the same material in plain language: what a breeding objective is, why
#' economic weights are difficult to state, how to read a feasibility result,
#' how to choose among the index families, and how to judge whether a
#' recommendation is robust. No R code is required to read the guide.
#'
#' It is intended for breeders, programme managers and reviewers who will use
#' or approve a selection decision but do not necessarily run R themselves.
#'
#' The tagged PDF provides a fixed reference edition and the HTML edition
#' renders the same content for on-screen reading. These standalone files have
#' no vignette engine, so they are not indexed by
#' \code{vignette()}/\code{browseVignettes()} like the package's \code{.Rmd}
#' vignettes. Their version-controlled source is
#' \code{inst/guide/DesiredGainR_Breeder_Guide.Rmd}; both generated files ship
#' in \code{inst/extdata/} and are located via \code{\link[base]{system.file}}.
#'
#' @param open Logical, default \code{TRUE}. When \code{TRUE} and the session
#'   is interactive, the file is opened with the operating system's default
#'   application through \code{\link[utils]{browseURL}}. When \code{FALSE}, or
#'   in a non-interactive session, only the path is returned.
#' @param format Character, one of \code{"pdf"} (default) or \code{"html"}.
#'
#' @return The file path, invisibly, to the selected edition in the installed
#'   package. An explicit error is raised when the requested edition is
#'   unavailable, rather than a silent fallback to the other format.
#'
#' @seealso The \emph{Defining a breeding objective} vignette
#'   (\code{vignette("DesiredGainR-objective", package = "DesiredGainR")}) for
#'   the same material with full statistical detail and runnable code;
#'   \code{\link{gain_feasibility}}, \code{\link{implied_economic_weights}},
#'   \code{\link{weight_sensitivity}}.
#'
#' @examples
#' # Locate without opening, which is what a script should do.
#' path <- try(open_desiredgain_guide(open = FALSE), silent = TRUE)
#' if (!inherits(path, "try-error")) file.exists(path)
#'
#' @export
open_desiredgain_guide <- function(open = TRUE, format = c("pdf", "html")) {
  format <- match.arg(format)
  guide_name <- paste0("DesiredGainR_Breeder_Guide.", format)
  guide_path <- system.file("extdata", guide_name, package = "DesiredGainR")

  if (!nzchar(guide_path) || !file.exists(guide_path)) {
    stop(
      guide_name, " was not found in the installed package, where it is ",
      "expected under inst/extdata/. Reinstall DesiredGainR, or obtain the ",
      "guide directly from the package repository.",
      call. = FALSE
    )
  }

  if (isTRUE(open) && interactive()) {
    message("[open_desiredgain_guide] Opening ", guide_name, " ...")
    utils::browseURL(guide_path)
  } else {
    message("[open_desiredgain_guide] Guide located at: ", guide_path)
  }

  invisible(guide_path)
}
