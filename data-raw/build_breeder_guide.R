# Render the Breeder's Guide into inst/extdata, where
# open_desiredgain_guide() looks for it.
#
# Run this before every release. The two generated editions are deliberately
# not committed as build artefacts of the vignette engine, because they are
# standalone documents rather than package vignettes, but they MUST exist in
# the installed package or open_desiredgain_guide() fails.
#
#   source("data-raw/build_breeder_guide.R")
#
# The PDF edition needs a LaTeX installation. Where none is available, build
# the HTML edition alone and let open_desiredgain_guide(format = "html") be
# the documented route.

source_file <- "inst/guide/DesiredGainR_Breeder_Guide.Rmd"
target_dir <- "inst/extdata"

stopifnot(file.exists(source_file))
dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)

render_edition <- function(format, output_file) {
  message("Rendering ", output_file, " ...")
  result <- tryCatch(
    {
      rmarkdown::render(
        input = source_file,
        output_format = format,
        output_file = output_file,
        output_dir = target_dir,
        quiet = TRUE
      )
      TRUE
    },
    error = function(e) {
      warning("Could not render ", output_file, ": ", conditionMessage(e),
        call. = FALSE
      )
      FALSE
    }
  )
  invisible(result)
}

html_ok <- render_edition("html_document", "DesiredGainR_Breeder_Guide.html")
pdf_ok <- render_edition("pdf_document", "DesiredGainR_Breeder_Guide.pdf")

if (!html_ok && !pdf_ok) {
  stop(
    "Neither edition rendered. open_desiredgain_guide() will fail in the ",
    "installed package. Fix this before release, or remove the export.",
    call. = FALSE
  )
}
if (!pdf_ok) {
  # Do not guess at the cause. The common failures are distinguishable from
  # the log, and naming the wrong one sends the reader to the wrong fix: an
  # L3 mismatch looks nothing like a missing installation but the previous
  # version of this message claimed it was one.
  message(
    "\nThe PDF edition did not render. Read the .log file named above before ",
    "acting; the usual causes are:\n\n",
    "  1. 'Mismatched LaTeX support files' / 'Loading expl3.sty aborted'\n",
    "     The LaTeX format and the package tree are out of step. TinyTeX IS\n",
    "     installed; it is stale. Fix with:\n",
    "         tinytex::reinstall_tinytex()\n",
    "     or, to keep the existing tree:\n",
    "         tinytex::tlmgr(c('update', '--self', '--all'))\n",
    "         tinytex::tlmgr(c('install', 'latex-bin', 'l3kernel',\n",
    "                          'l3backend', 'l3packages'))\n",
    "         tinytex::tlmgr(c('--verify-repo=none', 'update', '--all'))\n\n",
    "  2. 'Fatal error occurred, no output PDF file produced' with a missing\n",
    "     .sty name. A package is absent; tinytex normally installs it\n",
    "     automatically, but a repository that fails verification blocks that.\n",
    "     Set a verified mirror:\n",
    "         tinytex::tlmgr(c('option', 'repository',\n",
    "                          'https://mirror.ctan.org/systems/texlive/tlnet'))\n\n",
    "  3. No LaTeX at all: 'pdflatex not found' or similar. Then, and only\n",
    "     then:\n",
    "         tinytex::install_tinytex()\n\n",
    "  The HTML edition is the documented default for ",
    "open_desiredgain_guide(),\n",
    "  so the guide is usable meanwhile."
  )
}
message("\nDone.")
