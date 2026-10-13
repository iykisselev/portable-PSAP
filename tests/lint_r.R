#!/usr/bin/env Rscript

# Check readability without treating data.table column expressions as ordinary
# variable references. Keep the three selected checks here with their runner.
paths <- c("workflow/lib/r", "workflow/scripts/apply", "tests/r")
files <- c("tests/test_genotype_models.R", unlist(lapply(paths, function(path) {
  list.files(path, pattern = "\\.R$", full.names = TRUE)
})))
linters <- list(
  line_length_linter = lintr::line_length_linter(120L),
  trailing_whitespace_linter = lintr::trailing_whitespace_linter(),
  semicolon_linter = lintr::semicolon_linter()
)
lints <- unlist(lapply(files, lintr::lint, linters = linters), recursive = FALSE)
if (length(lints)) {
  print(lints)
  quit(status = 1L)
}
