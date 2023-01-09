test_that("Maelstrom_Motif2TF works", {
  # TPMfile = paste0(tempdir(),'/TPM.tsv')
  # Countfile = paste0(tempdir(),'/RNA_Counts.tsv')
  # on.exit(unlink(c(Countfile,TPMfile)))
  #
  # expect_false(file.exists(TPMfile))
  # expect_false(file.exists(Countfile))

  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_Maelstrom_obj.Rds"))

  for(method in c('max_cor','means')){
  sce_returned <- Maelstrom_Motif2TF(sce,
                                 cluster_id = 'seurat_clusters',
                                 maelstrom_dir = testthat::test_path(),
                                 combine_motifs = method,
                                 expr_tresh = 2,
                                 cor_tresh = 0.01,
                                 cor_method = "pearson",
                                 output_dir = tempdir())

  expect_true(length(sce_returned@assays) == 6)
  expect_true("MotifTFcor" %in% names(sce_returned@assays))
  expect_true("MotifTFanticor" %in% names(sce_returned@assays))}
})

