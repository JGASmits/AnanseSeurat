#test_functions

test_that("DEGS are called from the single cell object", {
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))
  DEG_dir <- paste0(tempdir(), '/deseq2')
  print(DEG_dir)
  on.exit(unlink(DEG_dir, recursive = T))
  
  DEGS_scANANSE(
    sce,
    min_cells = 2,
    output_dir = tempdir(),
    cluster_id = 'seurat_clusters',
    RNA_count_assay = "RNA",
    additional_contrasts = list('cluster1_cluster2')
  )
  expect_true(file.exists(
    paste0(
      DEG_dir,
      '/',
      'hg38-anansesnake_cluster1_average.diffexp.tsv'
    )
  ))
  expect_true(file.exists(
    paste0(
      DEG_dir,
      '/',
      'hg38-anansesnake_cluster2_average.diffexp.tsv'
    )
  ))
  expect_true(file.exists(
    paste0(
      DEG_dir,
      '/',
      'hg38-anansesnake_cluster3_average.diffexp.tsv'
    )
  ))
  expect_true(file.exists(
    paste0(
      DEG_dir,
      '/',
      'hg38-anansesnake_cluster1_cluster2.diffexp.tsv'
    )
  ))
})
