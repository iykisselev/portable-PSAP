variant_key <- function(chr, start, end, ref, alt) {
    paste0(chr, ":", start, "-", end, "_", ref, ">", alt)
}


add_vcf_variant_identity <- function(genotypes) {
    genotypes[, `:=`(
        Chr = trimws(as.character(CHROM)),
        Start = suppressWarnings(as.integer(POS)),
        End = suppressWarnings(as.integer(POS)),
        Ref = trimws(as.character(REF)),
        Alt = trimws(as.character(ALT))
    )]

    invalid <- is.na(genotypes$Chr) | !nzchar(genotypes$Chr) |
        is.na(genotypes$Start) | is.na(genotypes$End) |
        is.na(genotypes$Ref) | !nzchar(genotypes$Ref) |
        is.na(genotypes$Alt) | !nzchar(genotypes$Alt)

    if (any(invalid)) {
        abort(
            "Input contains %d rows with invalid variant coordinates or alleles",
            sum(invalid)
        )
    }

    genotypes[, Key := variant_key(Chr, Start, End, Ref, Alt)]
    genotypes
}

