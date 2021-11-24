

bioc_avail_primates <- function() {
    avl_pkg_names <- BSgenome::available.genomes()
    ucsc_names <-
        c("panTro", "panPan", "gorGor", "ponAbe", "nomLeu", "chlSab",
          "macFas", "rheMac", "papAnu", "papHam", "nasLar", "rhiRox",
          "calJac", "saiBol", "tarSyr", "micMur", "otoGar")
    lgl.match <- lapply(ucsc_names, function(ucsc_name)
        grepl(ucsc_name, avl_pkg_names, ignore.case = TRUE))
    lgl.match <- Reduce(lgl.match, init = FALSE, f = `|`)
    avl_pkg_names[lgl.match]
}

# #' @export
# #' @example
# #' vcf <- TOPMedSampleVcf()
# #' vcf
# TOPMedSampleVcf <- function() {
#     ans <- system.file("extdata/topmed_freeze5_sample.vcf.bgz", package = "AceSNP", mustWork = TRUE)
#     ans <- VariantAnnotation::VcfFile(ans)
#     ans <- VariantAnnotation::readVcf(ans)
#     ans
# }

# Extract information from vcf file and fix strand
extractVcf <- function(vcf) {
    if (is.character(vcf))
        vcf <- VariantAnnotation::VcfFile(vcf)
    if (is(vcf, "VcfFile"))
        vcf <- VariantAnnotation::readVcf(vcf, row.names = FALSE)
    stopifnot(is(vcf, "VCF"))

    #vcf <- expand(vcf)
    gr <- rowRanges(vcf)

    ## Remove row names
    names(gr) <- NULL

    ## Add extra columns
    mcols(gr) <- cbind(mcols(gr), info(vcf))

    ## Remove some columns that I do not care about
    ## - AF can be calculated from AN and AC
    for (colname in c("paramRangeID", "QUAL", "FILTER", "AF", "VRT"))
        mcols(gr)[[colname]] <- NULL

    ## Unlist some columns
    for (colname in c("ALT", "AC", "Het", "Hom")) {
        stopifnot(is(mcols(gr)[[colname]], "List"))
        stopifnot(unique(lengths(mcols(gr)[[colname]])) == 1)
        mcols(gr)[[colname]] <- unlist(mcols(gr)[[colname]])
    }

    ## Fix strand
    stopifnot(unique(strand(gr)) == "*")
    strand(gr) <- "+"

    # Filter only SNPs
    # FIXME: Note that when we calculate the frequency of the reference allele,
    # we need to taking in consideration of the deletino... It's complex...
    gr <- gr[width(gr) == 1 & width(gr$ALT) == 1]

    gr
}

collapseRanges <- function(gr) {
    ans <- gr
    mcols(ans) <- NULL
    rownames(ans) <- NULL

    ans <- unique(ans)
    mcols(ans) <- split()
    ans <- ans[unique(ans)]
    ans
}

if (FALSE) {
    vcf <- readVcf(
        VcfFile("../TreeFriends/raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz"),
        param = ScanVcfParam(which = GRanges("chr8", IRanges(1, 1000000)))
    )
    gr <- extractVcf(vcf)
    gr <- collapseRanges(gr)
    gr
}


