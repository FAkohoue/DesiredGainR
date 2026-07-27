# DesiredGainR GitHub setup script
# Run this from an interactive R session after editing placeholders.

# Install needed packages if missing:
# install.packages(c("usethis", "devtools", "roxygen2", "pkgdown", "testthat"))

library(usethis)
library(devtools)

# ------------------------------------------------------------------
# 1) Set your working directory to the package root before sourcing
# ------------------------------------------------------------------
# setwd("/path/to/DesiredGainR")

# ------------------------------------------------------------------
# 2) Edit package metadata placeholders before proceeding
# ------------------------------------------------------------------
# Replace <YOUR_GITHUB_USERNAME> in:
# - DESCRIPTION
# - README.md
# - _pkgdown.yml

# Optional: set your Git identity if needed
# use_git_config(user.name = "Your Name", user.email = "you@example.com")

# ------------------------------------------------------------------
# 3) Regenerate docs and check package
# ------------------------------------------------------------------
#unlink("NAMESPACE")
devtools::document()
devtools::document()
devtools::check()

# ------------------------------------------------------------------
# 4) Initialize git and make first commit
# ------------------------------------------------------------------
use_git()

# ------------------------------------------------------------------
# 5) Create GitHub repo and push
# ------------------------------------------------------------------
# This requires a GitHub PAT configured for usethis/gh.
# If needed, run: usethis::create_github_token()
# Then store it: gitcreds::gitcreds_set()
use_github(private = FALSE, protocol = "https")

# ------------------------------------------------------------------
# 6) Optional: pkgdown site
# ------------------------------------------------------------------
# build_site()

# ------------------------------------------------------------------
# 7) Optional: install locally
# ------------------------------------------------------------------
install(build_vignettes = TRUE)

message("Setup complete. Review GitHub Actions and pkgdown deployment in your repository.")
