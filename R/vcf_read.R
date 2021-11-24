
#' Sample records from vcf file
#'
#' \code{sampleVcf} randomly sample records from a large vcf file (i.e. a vcf file that is hard to read
#' into memory) for analysis. It wraps the \code{filterVcf} function from \code{VariantAnnotation}
#' package.
#'
#' @param file,destination,verbose,index,
#' @param prob
#' @param seed
#' @param yieldSize
#'
#' @return
#' @export
#'
#' @examples
#' vcf <- VcfFile(system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation"))
#' outvcf <- sampleVcf(vcf, tempfile(fileext = ".vcf"), index = FALSE)
#' print(outvcf)
#' system2("ls", c("-lh", path(outvcf)))
sampleVcf <- function(file, destination, prob = 0.001, seed = 42,
                      verbose = FALSE, index = FALSE, yieldSize = 10000L) {

    ## FIXME: Note that the records are randomly selected. Thus we miss
    ## records on the same location.

    set.seed(seed)
    file <- VariantAnnotation::VcfFile(file)
    # Set yield size
    original_yieldsize <- Rsamtools::yieldSize(file)
    Rsamtools::yieldSize(file) <- yieldSize
    # Restore yield size when exit
    on.exit(Rsamtools::yieldSize(file) <- original_yieldsize)

    ans <- VariantAnnotation::filterVcf(
        file,
        destination = destination,
        verbose = verbose,
        index = index,
        prefilters = S4Vectors::FilterRules(
            function(x) {
                runif(length(x)) < prob
            }
        )
    )
    VariantAnnotation::VcfFile(ans)
}

if (FALSE) {
    vcf_file <- VariantAnnotation::VcfFile("../TreeFriends/raw_data/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz")
    sampleVcf(vcf_file, destination = "topmed_sample.vcf", index = FALSE, prob = 0.001,
              verbose = TRUE, yieldSize = 1000000L)
    bgzip("topmed_sample.vcf")
    indexTabix("topmed_sample.vcf.bgz", format = "vcf")
}
