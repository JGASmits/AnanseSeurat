#test_functions

test_that("import seurat succesfully imports an influence file, wich contains an influence assay",
          {
            #load dummy single cell object
            sce <- readRDS(testthat::test_path("sce_obj.Rds"))
            sce <- import_seurat_scANANSE(
              sce,
              cluster_id = 'seurat_clusters',
              anansnake_inf_dir = testthat::test_path('influence/')
            )
            expect_true(length(sce@assays) == 3)
            expect_true('influence' %in% names(sce@assays))
            inf_df <- as.data.frame(sce@assays$influence@data)
            inf_correct <-
              read.table(testthat::test_path("influence_assay.tsv"))
            expect_true(all.equal(inf_df, inf_correct, tolerance = 1.0e-7))
            
          })

test_that(
  "import maelstrom succesfully imports a maelstrom file, wich contains an motif enrichment scores",
  {
    #load dummy single cell object
    sce <- readRDS(testthat::test_path("sce_obj.Rds"))
    sce <- import_seurat_maelstrom(
      sce,
      cluster_id = 'seurat_clusters',
      maelstrom_file = testthat::test_path('final.out.txt')
    )
    expect_true(length(sce@assays) == 3)
    expect_true('maelstrom' %in% names(sce@assays))
    maelstrom_df <- as.data.frame(sce@assays$maelstrom@data)
    maelstrom_correct <-
      read.table(testthat::test_path("maelstrom_table.tsv"))
    expect_true(all.equal(maelstrom_df, maelstrom_correct))
  }
)

test_that("per cell df works", {
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_inf_obj.Rds"))
  TF_influence <-
    per_cluster_df(sce, assay = 'influence', cluster_id = 'seurat_clusters')
  expect_true(all.equal(dim(TF_influence), c(4, 3)))
  
  TF_influence_correct <-
    read.table(testthat::test_path("influence_table.tsv"))
  expect_true(all.equal(TF_influence, TF_influence_correct))
})
