test_that("doubletFinder_v3 works", {
  library("Seurat")
  # Make test deterministic
  set.seed(0)
  # Err if wrong length, type, or if any missingness
  testthat::expect_error(
    doubletFinder_v3(pbmc_small, PCs = 1:4, pN = 0.25, pK = 0.01, nExp = 5, reuse.pANN = FALSE, sct=FALSE,
                     annotations = "blah")
  )
  testthat::expect_error(
    doubletFinder_v3(pbmc_small, PCs = 1:4, pN = 0.25, pK = 0.01, nExp = 5, reuse.pANN = FALSE, sct=FALSE,
                     annotations = as.factor(pbmc_small$RNA_snn_res.1))
  )
  testthat::expect_error(
    doubletFinder_v3(pbmc_small, PCs = 1:4, pN = 0.25, pK = 0.01, nExp = 5, reuse.pANN = FALSE, sct=FALSE,
                     annotations = as.character(NA*as.numeric(pbmc_small$RNA_snn_res.1)))
  )
  testthat::expect_output(
    pbmc_small_dubs <- doubletFinder_v3(pbmc_small, PCs = 1:4, pN = 0.25, pK = 0.01, nExp = 5, reuse.pANN = FALSE, sct=FALSE,
                                        annotations = as.character(pbmc_small$RNA_snn_res.1))
  )
  testthat::expect_output(
    pbmc_small_dubs <- doubletFinder_v3(pbmc_small, PCs = 1:4, pN = 0.5, pK = 0.01, nExp = 5, reuse.pANN = FALSE, sct=FALSE,
                                        annotations = as.character(pbmc_small$RNA_snn_res.1))
  )
  testthat::expect_length( pbmc_small_dubs$DF.doublet.contributors_0.5_0.01_5_0, 80 )
  testthat::expect_length( pbmc_small_dubs$DF.doublet.contributors_0.5_0.01_5_1, 80 )
  testthat::expect_length( pbmc_small_dubs$DF.doublet.contributors_0.5_0.01_5_2, 80 )
})
