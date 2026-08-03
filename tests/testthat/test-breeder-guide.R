# Item 11 of the external review: the Breeder's Guide must actually ship.
#
# open_desiredgain_guide() is exported and promises a file. Until this test
# existed nothing checked that the file was ever built, and the function was
# broken in the installed package while looking fine in the source tree.

test_that("the guide source ships in the installed package", {
  source_path <- system.file(
    "guide", "DesiredGainR_Breeder_Guide.Rmd",
    package = "DesiredGainR"
  )
  expect_true(nzchar(source_path))
  expect_true(file.exists(source_path))
})

test_that("at least one rendered edition resolves through system.file()", {
  # Rendering happens in data-raw/build_breeder_guide.R and needs pandoc, and
  # for the PDF a LaTeX installation, so it is a release step rather than a
  # check-time one. The test is skipped where neither edition is present, but
  # the skip message names the fix so the omission cannot pass silently.
  editions <- vapply(
    c("html", "pdf"),
    function(format) {
      system.file(
        "extdata", paste0("DesiredGainR_Breeder_Guide.", format),
        package = "DesiredGainR"
      )
    },
    character(1L)
  )
  present <- editions[nzchar(editions) & file.exists(editions)]
  skip_if(
    !length(present),
    paste(
      "No rendered Breeder's Guide is installed. Run",
      "data-raw/build_breeder_guide.R before release, or remove",
      "open_desiredgain_guide() and its website links."
    )
  )
  expect_true(all(file.exists(present)))
  expect_gt(min(file.size(present)), 1000)
})

test_that("a missing edition errors rather than returning a bad path", {
  available <- vapply(
    c("pdf", "html"),
    function(format) {
      path <- system.file(
        "extdata", paste0("DesiredGainR_Breeder_Guide.", format),
        package = "DesiredGainR"
      )
      nzchar(path) && file.exists(path)
    },
    logical(1L)
  )
  missing_format <- names(available)[!available]
  skip_if(!length(missing_format), "Both editions are installed.")
  expect_error(
    open_desiredgain_guide(open = FALSE, format = missing_format[1L]),
    "was not found in the installed package"
  )
})

test_that("a located guide is returned invisibly and exists", {
  path <- tryCatch(
    open_desiredgain_guide(open = FALSE),
    error = function(e) NULL
  )
  skip_if(is.null(path), "The default (HTML) edition is not installed.")
  expect_true(file.exists(path))
})

test_that("the default format is the one that renders most easily", {
  # The PDF edition needs LaTeX as well as pandoc, so defaulting to it made the
  # exported function fail on installations where the guide had been built.
  expect_identical(eval(formals(open_desiredgain_guide)$format)[1L], "html")
})

test_that("a missing edition points at the one that is present", {
  installed <- vapply(
    c("html", "pdf"),
    function(format) {
      path <- system.file(
        "extdata", paste0("DesiredGainR_Breeder_Guide.", format),
        package = "DesiredGainR"
      )
      nzchar(path) && file.exists(path)
    },
    logical(1L)
  )
  skip_if(all(installed) || !any(installed), "Needs exactly one edition.")
  absent <- names(installed)[!installed]
  expect_error(
    open_desiredgain_guide(open = FALSE, format = absent),
    "edition is installed"
  )
})
