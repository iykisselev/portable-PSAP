# Canonical genomic identifiers shared with workflow/lib/python/genomic_contracts.py.

canonical_chromosome <- function(chromosome) {
  value <- trimws(as.character(chromosome))
  value <- sub("^chr", "", value, ignore.case = TRUE)
  value <- toupper(value)
  value[value == "M"] <- "MT"
  value
}


canonical_allele_key <- function(chromosome, position, reference, alternate) {
  paste(
    canonical_chromosome(chromosome),
    as.integer(position),
    toupper(as.character(reference)),
    toupper(as.character(alternate)),
    sep = ":"
  )
}


clinical_variant_key <- function(chromosome, start, end, reference, alternate) {
  paste0(
    canonical_chromosome(chromosome), ":", as.integer(start), "-", as.integer(end),
    "_", toupper(as.character(reference)), ">", toupper(as.character(alternate))
  )
}
