standardize_genotypes <- function(genotypes, hemizygote_model) {
  genotypes[
    ,
    GT := chartr("|", "/", trimws(as.character(GT)))
  ]

  missing <- is.na(genotypes$GT) |
    !nzchar(genotypes$GT) |
    grepl(".", genotypes$GT, fixed = TRUE)
  diploid <- !missing & grepl("^[0-9]+/[0-9]+$", genotypes$GT)
  haploid <- !missing & grepl("^[0-9]+$", genotypes$GT)
  invalid <- !missing & !diploid & !haploid

  if (any(invalid)) {
    examples <- head(unique(genotypes$GT[invalid]), 5L)
    abort(
      "Input contains invalid genotype value(s): %s",
      paste(examples, collapse = ", ")
    )
  }

  if (any(haploid) && !hemizygote_model) {
    abort(
      "Input contains haploid genotype calls, but the hemizygote model is disabled"
    )
  }

  if (any(haploid)) {
    chromosome <- canonical_chromosome(genotypes$CHROM)
    unsupported <- haploid & !(chromosome %chin% c("X", "Y"))
    if (any(unsupported)) {
      examples <- head(unique(chromosome[unsupported]), 5L)
      abort(
        "Haploid genotype calls are supported only on chromosomes X and Y; found: %s",
        paste(examples, collapse = ", ")
      )
    }
  }

  diploid_rows <- which(diploid)
  if (length(diploid_rows) > 0L) {
    allele_parts <- tstrsplit(
      genotypes$GT[diploid_rows],
      "/",
      fixed = TRUE,
      type.convert = TRUE
    )
    first_allele <- pmin(allele_parts[[1]], allele_parts[[2]])
    second_allele <- pmax(allele_parts[[1]], allele_parts[[2]])
    genotypes$GT[diploid_rows] <- paste0(
      first_allele,
      "/",
      second_allele
    )
  }

  keep <- !missing & genotypes$GT != "0" & genotypes$GT != "0/0"
  genotypes <- genotypes[keep]

  if (nrow(genotypes) == 0L) {
    genotypes[, gt_class := character()]
    return(genotypes)
  }

  is_haploid <- grepl("^[0-9]+$", genotypes$GT)
  genotypes[, gt_class := "haploid"]

  diploid_rows <- which(!is_haploid)
  if (length(diploid_rows) > 0L) {
    allele_parts <- tstrsplit(
      genotypes$GT[diploid_rows],
      "/",
      fixed = TRUE
    )
    genotypes$gt_class[diploid_rows] <- fifelse(
      allele_parts[[1]] == allele_parts[[2]],
      "hom",
      "het"
    )
  }

  genotypes
}


apply_sex_classification <- function(genotypes, sex, hemizygote_model) {
  chromosome <- canonical_chromosome(genotypes$CHROM)

  gt_class <- tolower(trimws(as.character(genotypes$gt_class)))
  is_hom <- gt_class %chin% c("hom", "homozygous")
  is_het <- gt_class == "het"
  is_haploid <- gt_class == "haploid"
  is_x <- chromosome == "X"
  is_y <- chromosome == "Y"
  is_xy <- is_x | is_y

  if (any(is_haploid & !is_xy)) {
    abort("Haploid genotype calls are supported only on chromosomes X and Y")
  }
  if (any(is_haploid) && !hemizygote_model) {
    abort(
      "Haploid genotype calls require the hemizygote model to be enabled"
    )
  }

  remove <- rep(FALSE, nrow(genotypes))
  removed_male_xy_het <- 0
  removed_unknown_sex_xy <- 0

  if (!is.na(sex) && sex == 1) {
    if (!hemizygote_model) {
      remove <- is_xy
    } else {
      # The generated HEM contract is chromosome-X only.  Do not reinterpret
      # Y calls as X-linked hemizygous evidence.
      remove <- is_y | (is_x & is_het)
      removed_male_xy_het <- sum(is_x & is_het)
      genotypes[
        is_x & (is_hom | is_haploid),
        gt_class := "hem"
      ]
    }
  } else if (!is.na(sex) && sex == 2) {
    if (any(is_haploid)) {
      abort("Haploid X/Y genotype calls cannot be classified for a female sample")
    }
    genotypes[
      is_x & is_het,
      gt_class := "het_female"
    ]
    genotypes[
      is_x & is_hom,
      gt_class := "hom_female"
    ]
  } else {
    if (any(is_haploid)) {
      abort("Haploid X/Y genotype calls require a known male sample sex")
    }
    # Sex is required to distinguish female-X diploid models from male-X HEM.
    # Unknown-sex X/Y calls are retained nowhere rather than leaking into the
    # autosomal DOM/REC paths.
    remove <- is_xy
    removed_unknown_sex_xy <- sum(remove)
  }

  list(
    genotypes = genotypes[!remove],
    removed_male_xy_het = removed_male_xy_het,
    removed_unknown_sex_xy = removed_unknown_sex_xy
  )
}
