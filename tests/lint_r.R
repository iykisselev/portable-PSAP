#!/usr/bin/env Rscript

# Check readability without treating data.table column expressions as ordinary
# variable references. The selected linters are configured in .lintr.
paths <- c("workflow/lib/r", "workflow/scripts/apply", "tests/r")
files <- c("tests/test_genotype_models.R", unlist(lapply(paths, function(path) {
  list.files(path, pattern = "\\.R$", full.names = TRUE)
})))
lints <- unlist(lapply(files, lintr::lint), recursive = FALSE)
if (length(lints)) {
  print(lints)
  quit(status = 1L)
}
