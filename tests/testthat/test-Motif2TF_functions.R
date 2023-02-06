test_that("Maelstrom_Motif2TF works", {
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_Maelstrom_obj.Rds"))
  
  for (method in c('max_cor', 'means', 'max_var')) {
    sce_returned <- Maelstrom_Motif2TF(
      sce,
      cluster_id = 'seurat_clusters',
      maelstrom_dir = testthat::test_path(),
      combine_motifs = method,
      expr_tresh = 2,
      cor_tresh = 0.01,
      cor_method = "pearson"
    )
    expect_true(length(sce_returned@assays) == 6)
    expect_true("MotifTFcor" %in% names(sce_returned@assays))
    expect_true("MotifTFanticor" %in% names(sce_returned@assays))
  }
})

test_that("Factor_Motif_Plot works", {
  sce <- readRDS(testthat::test_path("sce_MotifTF_obj.Rds"))
  plot <- Factor_Motif_Plot(
    sce,
    c('gene1', 'gene2'),
    assay_maelstrom = 'MotifTFanticor',
    dim_reduction = 'pca',
    logo_dir = testthat::test_path('logos/'),
    
  )
  expect_true(ggplot2::is.ggplot(plot))
})
