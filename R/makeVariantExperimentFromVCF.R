#' The function to convert VCF files directly into VariantExperiment
#' object.
#' @name makeVariantExperimentFromVCF
#' @rdname makeVariantExperimentFromVCF
#' @description \code{makeVariantExperimentFromVCF} is the function
#'     to convert a vcf file into \code{VariantExperiment} object. The
#'     genotype data will be written as \code{GDSArray} format, which
#'     is saved in the \code{assays} slot. The annotation info for
#'     variants or samples will be written as \code{DelayedDataFrame}
#'     object, and saved in the \code{rowData} or \code{colData} slot.
#' @param vcf.fn the file name(s) of (compressed) VCF format; or a
#'     ‘connection’ object.
#' @param out.dir The directory to save the gds format of the vcf
#'     data, and the newly generated VariantExperiment object with
#'     array data in \code{GDSArray} format and annotation data in
#'     \code{DelayedDataFrame} format. The default is a temporary
#'     folder.
#' @param replace Whether to replace the directory if it already
#'     exists. The default is FALSE.
#' @param header if NULL, ‘header’ is set to be
#'     ‘seqVCF_Header(vcf.fn)’, which is a list (with a class name
#'     "SeqVCFHeaderClass", S3 object).
#' @param info.import characters, the variable name(s) in the INFO
#'     field for import; default is ‘NULL’ for all variables.
#' @param fmt.import characters, the variable name(s) in the FORMAT
#'     field for import; default is ‘NULL’ for all variables.
#' @param sample.info characters (with) file path for the sample info
#'     data. The data must have colnames (for phenotypes), rownames
#'     (sample ID's). No blank line allowed. The default is ‘NULL’ for
#'     no sample info.
#' @param ignore.chr.prefix a vector of character, indicating the
#'     prefix of chromosome which should be ignored, like "chr"; it is
#'     not case-sensitive.
#' @param reference genome reference, like "hg19", "GRCh37"; if the
#'     genome reference is not available in VCF files, users could
#'     specify the reference here.
#' @param start the starting variant if importing part of VCF files.
#' @param count the maximum count of variant if importing part of VCF
#'     files, -1 indicates importing to the end.
#' @param parallel ‘FALSE’ (serial processing), ‘TRUE’ (parallel
#'     processing), a numeric value indicating the number of cores, or
#'     a cluster object for parallel processing; ‘parallel’ is passed
#'     to the argument ‘cl’ in ‘seqParallel’, see
#'     ‘?SeqArray::seqParallel’ for more details. The default is
#'     "FALSE".
#' @param verbose whether to print the process messages. The default
#'     is FALSE.
#' @return An \code{VariantExperiment} object.
#' @importFrom utils read.table
#' @export
#' @examples
#' ## the vcf file
#' vcf <- SeqArray::seqExampleFileName("vcf")
#' ## conversion
#' ve <- makeVariantExperimentFromVCF(vcf)
#' ve
#' ## the filepath to the gds file.
#' gdsfile(ve)
#' 
#' ## only read in specific info columns
#' ve <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(),
#'                                    info.import=c("OR", "GP"))
#' ve
#' ## convert without the INFO and FORMAT fields
#' ve <- makeVariantExperimentFromVCF(vcf, out.dir = tempfile(),
#'                                    info.import=character(0),
#'                                    fmt.import=character(0))
#' ve
#' ## now the assay data does not include the
#' #"annotation/format/DP/data", and the rowData(ve) does not include
#' #any info columns.
 
makeVariantExperimentFromVCF <- function(vcf.fn,
                                         out.dir = tempfile(),
                                         replace = FALSE,
                                         header = NULL,
                                         info.import = NULL,
                                         fmt.import = NULL,
                                         sample.info = NULL,
                                         ignore.chr.prefix = "chr",
                                         reference = NULL,
                                         start = 1L, count = -1L,
                                         parallel = FALSE,
                                         verbose = FALSE)
{
    stopifnot(is.character(vcf.fn), length(vcf.fn)==1L)
    if (!isSingleString(out.dir))
        stop(wmsg("'dir' must be a single string specifying the path ",
                  "to the directory where to save the \"VariantExperiment\" ",
                  "object (the directory will be created)"))
    
    ## stopifnot(is.character(out.dir), length(out.dir)==1L)
    
    ## run VCF to GDS
    .create_dir(out.dir, replace)
    out.gds.fn <- file.path(out.dir, "se.gds")
    seqVCF2GDS(vcf.fn, out.gds.fn, header = header,
               storage.option = "LZMA_RA", 
               info.import = info.import, fmt.import = fmt.import,
               genotype.var.name = "GT",
               ignore.chr.prefix = ignore.chr.prefix,
               reference = reference, start = start, count = count,
               optimize=TRUE, raise.error=TRUE, digest = TRUE,
               parallel = parallel,
               verbose=verbose)
    
    ## add sample info into gds file.
    if (!is.null(sample.info)) {
        sample.info <- read.table(file = sample.info, header = TRUE,
                                  blank.lines.skip = TRUE,
                                  stringsAsFactors = FALSE)
        .addsmpanno(out.gds.fn, value = sample.info)
    }
    
    ## run GDS to VE
    makeVariantExperimentFromGDS(file=out.gds.fn,
                                 assayNames = NULL,
                                 rowDataOnDisk = TRUE,
                                 colDataOnDisk = TRUE, 
                                 infoColumns = info.import)
}

.addsmpanno <- function(gdsfile, value) {
    gfile <- acquireGDS(gdsfile, readonly = FALSE)
    sample.id <- read.gdsn(index.gdsn(gfile, "sample.id"))
    idx <- match(sample.id, rownames(value))
    if (any(is.na(idx)))
        warning(paste0(sum(is.na(idx)), " out of ",
                       length(sample.id), " samples ",
                       "does not have sample info available. ",
                       "(", paste(head(sample.id[is.na(idx)]), collapse=", "),
                       ifelse(sum(is.na(idx)) > 6, "...", ""), ")"))
    value <- value[idx, , drop=FALSE]
    sampAnnot <- index.gdsn(gfile, "sample.annotation")
    for (i in seq_along(value)) {
        SeqArray:::.AddVar(storage.option = SeqArray::seqStorageOption("LZMA_RA"),
                           node = sampAnnot, varname = names(value)[i],
                           val = value[,i],
                           closezip = TRUE)
        ## FIXME: returns warning: Missing characters are converted to "".
    }
}
    
