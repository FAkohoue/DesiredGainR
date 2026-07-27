# DesiredGainR GitHub publish checklist

1. Replace `<YOUR_GITHUB_USERNAME>` in `DESCRIPTION`, `README.md`, and `_pkgdown.yml`.
2. Open R in the package root.
3. Run `source("scripts/setup_github_repo.R")` line by line.
4. If GitHub authentication is missing, create a PAT with `usethis::create_github_token()`.
5. Save the PAT with `gitcreds::gitcreds_set()`.
6. Confirm the first push completed successfully.
7. Check the **Actions** tab in GitHub for `R-CMD-check` and `pkgdown`.
8. If pkgdown does not deploy automatically, verify repository **Settings > Pages** and workflow permissions.
9. After the repository is live, install from GitHub with:

```r
remotes::install_github("<YOUR_GITHUB_USERNAME>/DesiredGainR", build_vignettes = TRUE)
```
