test_that("initiate seq gds file works", {
    file <- SeqArray::seqExampleFileName("gds")
    .initiate_seqgds <- VariantExperiment:::.initiate_seqgds

    ve <- makeVariantExperimentFromGDS(file)
    gds_path <- tempfile(fileext=".gds")
    .initiate_seqgds(ve, gds_path, compress="LZMA_RA")
    ff <- acquireGDS(gds_path)
    ## expect_true(validObject(f))
    ## Error in `validObject(.Object)`: invalid class “SeqVarGDSClass” object: sample.id, variant.id, position, chromosome, allele, and genotype are required variables.
    expect_true(is(ff, "gds.class"))
    expect_true(all(c("sample.id", "variant.id", "chromosome", "position") %in%
                    gdsfmt::ls.gdsn(ff)))
})

test_that("initiate snp gds file works", {
    .initiate_snpgds <- VariantExperiment:::.initiate_snpgds
    file <- SNPRelate::snpgdsExampleFileName()

    ve <- makeVariantExperimentFromGDS(file)

    gds_path <- tempfile(fileext=".gds")
    .initiate_snpgds(ve, gds_path, compress="LZMA_RA")
    ff <- acquireGDS(gds_path)
    expect_true(validObject(ff))
    expect_true(is(ff, "gds.class"))
    expect_true(all(c("sample.id", "snp.id", "snp.chromosome", "snp.position") %in%
                    gdsfmt::ls.gdsn(ff)))
})

test_that("write GDSSE works", {
    file <- SeqArray::seqExampleFileName("gds")

    .initiate_seqgds <- VariantExperiment:::.initiate_seqgds
    .write_ve_as_gds <- VariantExperiment:::.write_ve_as_gds
    .write_ve_as_newve <- VariantExperiment:::.write_ve_as_newve

    ve <- makeVariantExperimentFromGDS(file, rowDataOnDisk=FALSE, colDataOnDisk=FALSE)
    gds_path <- tempfile(fileext=".gds")
    .initiate_seqgds(ve, gds_path, compress="LZMA_RA")
    suppressWarnings(
        .write_ve_as_gds(ve, "SEQ_ARRAY", gds_path, chunk_size = 10000,
                     compress = "LZMA_RA", verbose = FALSE))
    ve1 <- .write_ve_as_newve(ve, gds_path, "SEQ_ARRAY", "variant.id", "sample.id", TRUE, TRUE)
    expect_true(validObject(ve1))
    expect_s4_class(ve1, "VariantExperiment")
    expect_identical(dim(ve), dim(ve1))
})
