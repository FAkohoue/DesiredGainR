# Release checks that require a network connection or a full toolchain, and so
# cannot run inside R CMD check.
#
#   source("data-raw/release_checklist.R")
#
# Every step reports pass or fail rather than stopping at the first problem, so
# one run gives the whole list of what still needs doing.

release_state <- new.env(parent = emptyenv())
release_state$results <- list()
record <- function(name, passed, detail = "") {
  release_state$results[[length(release_state$results) + 1L]] <- list(
    name = name, passed = isTRUE(passed), detail = detail
  )
  invisible(NULL)
}

message("Running release checks ...\n")

# --- Version agreement ------------------------------------------------------

description_version <- as.character(read.dcf("DESCRIPTION")[1L, "Version"])
cff <- readLines("CITATION.cff", warn = FALSE)
cff_version <- sub(
  '^version:\\s*"?([^"]*)"?\\s*$', "\\1",
  grep("^version:", cff, value = TRUE)[1L]
)
record(
  "DESCRIPTION and CITATION.cff versions agree",
  identical(description_version, cff_version),
  paste("DESCRIPTION", description_version, "vs CITATION.cff", cff_version)
)
# inst/CITATION derives its version from package metadata, so it cannot drift;
# this confirms the file still does that rather than hard-coding a version.
citation_text <- readLines("inst/CITATION", warn = FALSE)
record(
  "inst/CITATION derives its version rather than hard-coding it",
  any(grepl("packageVersion", citation_text)),
  "A hard-coded version silently goes stale, as 0.3.1 did."
)

# --- The Breeder Guide ------------------------------------------------------

guide_html <- file.path("inst", "extdata", "DesiredGainR_Breeder_Guide.html")
record(
  "The Breeder Guide is rendered",
  file.exists(guide_html),
  "Run data-raw/build_breeder_guide.R"
)

# --- Legacy example data ----------------------------------------------------

legacy <- file.path(
  "inst", "extdata", c("example_gebv.csv", "example_pheno.csv")
)
record(
  "Legacy example CSVs are removed",
  !any(file.exists(legacy)),
  paste(
    "They are .Rbuildignore'd so they will not ship, but delete them from the",
    "repository once nothing local depends on them."
  )
)

# --- URLs and identifiers ---------------------------------------------------

if (requireNamespace("urlchecker", quietly = TRUE)) {
  urls <- tryCatch(urlchecker::url_check(), error = function(e) e)
  record(
    "All URLs resolve",
    !inherits(urls, "error") && !nrow(as.data.frame(urls)),
    if (inherits(urls, "error")) conditionMessage(urls) else "See output above."
  )
} else {
  record("All URLs resolve", FALSE, "install.packages('urlchecker')")
}

orcid <- "0000-0002-2160-0182"
record(
  "ORCID is stated",
  any(grepl(orcid, readLines("DESCRIPTION", warn = FALSE))),
  paste0("Verify manually at https://orcid.org/", orcid)
)

# --- Documentation is current ----------------------------------------------

if (requireNamespace("roxygen2", quietly = TRUE)) {
  # The package documents LazyData objects, so pkgload must make those objects
  # available while roxygen resolves the data documentation blocks.
  roxygen2::roxygenise(load_code = "pkgload")
  changed <- system("git status --porcelain man NAMESPACE", intern = TRUE)
  record(
    "Generated documentation is committed",
    !length(changed),
    paste(changed, collapse = "; ")
  )
} else {
  record(
    "Generated documentation is committed", FALSE,
    "install.packages('roxygen2')"
  )
}

# --- Full check -------------------------------------------------------------

if (requireNamespace("rcmdcheck", quietly = TRUE) &&
  requireNamespace("pkgbuild", quietly = TRUE)) {
  tarball <- pkgbuild::build(dest_path = tempdir())
  check <- rcmdcheck::rcmdcheck(tarball, args = "--as-cran", error_on = "never")
  expected_new_submission_note <-
    length(check$notes) == 1L &&
      grepl("New submission", paste(check$notes, collapse = "\n"), fixed = TRUE)
  record(
    "R CMD check --as-cran has no unexpected findings",
    !length(check$errors) && !length(check$warnings) &&
      (!length(check$notes) || expected_new_submission_note),
    sprintf(
      "%d error(s), %d warning(s), %d note(s)",
      length(check$errors), length(check$warnings), length(check$notes)
    )
  )
} else {
  record(
    "R CMD check --as-cran is clean", FALSE,
    "install.packages(c('rcmdcheck', 'pkgbuild'))"
  )
}

# --- Manual steps -----------------------------------------------------------

manual <- c(
  "Tag the release and push the tag, so Zenodo mints a DOI from .zenodo.json",
  "Add the DOI badge to README.md and the DOI to CITATION.cff",
  "Submit to CRAN with devtools::release()",
  "Submit the method paper (see paper/paper.md)"
)

# --- Report -----------------------------------------------------------------

cat("\n", strrep("-", 72), "\n", sep = "")
for (entry in release_state$results) {
  cat(
    if (entry$passed) "PASS  " else "TODO  ", entry$name, "\n",
    sep = ""
  )
  if (!entry$passed && nzchar(entry$detail)) {
    cat("        ", entry$detail, "\n", sep = "")
  }
}
cat("\nManual steps, in order:\n")
for (step in seq_along(manual)) {
  cat("  ", step, ". ", manual[step], "\n", sep = "")
}
cat(strrep("-", 72), "\n", sep = "")

outstanding <- sum(!vapply(
  release_state$results, `[[`, logical(1L), "passed"
))
if (outstanding) {
  stop("\n", outstanding, " automated check(s) still outstanding.",
    call. = FALSE
  )
} else {
  message("\nAll automated checks pass. The manual steps remain.")
}
