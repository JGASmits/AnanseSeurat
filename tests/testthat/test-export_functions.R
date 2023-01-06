#test_functions

test_that("CPM & Count files are sucessfully generated", {
  TPMfile = paste0(tempdir(),'/TPM.tsv')
  Countfile = paste0(tempdir(),'/RNA_Counts.tsv')
  on.exit(unlink(c(Countfile,TPMfile)))

  expect_false(file.exists(TPMfile))
  expect_false(file.exists(Countfile))

  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))

  export_CPM_scANANSE(sce,
                      min_cells <- 1,
                      output_dir =tempdir(),
                      cluster_id = 'seurat_clusters',
                      RNA_count_assay = 'RNA')
  expect_true(file.exists(TPMfile))
  expect_true(file.exists(Countfile))

  #check if the files are equal to the expected output
  TPM_correct <- read.table(testthat::test_path("TPM.tsv"))
  COUNT_correct <- read.table(testthat::test_path("RNA_Counts.tsv"))

  TPM_df <- read.table(TPMfile)
  COUNT_df <- read.table(Countfile)

  expect_true(all.equal(TPM_df, TPM_correct))
  expect_true(all.equal(COUNT_df, COUNT_correct))
})

test_that("Peak count matrix file is sucessfully generated", {
  Peakfile = paste0(tempdir(),'/Peak_Counts.tsv')
  on.exit(unlink(c(Peakfile)))

  expect_false(file.exists(Peakfile))

  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))

  export_ATAC_scANANSE(sce,
                       min_cells <- 1,
                       output_dir =tempdir(),
                       cluster_id = 'seurat_clusters',
                       ATAC_peak_assay = 'peaks')

  expect_true(file.exists(Peakfile))

  #check if the files are equal to the expected output
  Peak_correct <- read.table(testthat::test_path("Peak_Counts.tsv"))
  Peak_df <- read.table(Peakfile)
  expect_true(all.equal(Peak_df, Peak_correct))
})

test_that("Config and sample file are sucessfully generated", {
  Configfile = paste0(tempdir(),'/config.yaml')
  samplefile = paste0(tempdir(),'/samplefile.tsv')

  on.exit(unlink(c(Configfile,samplefile)))

  expect_false(file.exists(Configfile))
  expect_false(file.exists(samplefile))


  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))

  contrasts = list('cluster1_cluster2',
                   'cluster2_cluster1',
                   'cluster3_cluster2')
  config_scANANSE(sce,
                  min_cells <- 1,
                  output_dir = tempdir(),
                  cluster_id = 'seurat_clusters',
                  genome = paste0(tempdir(),'/hg38'),
                  additional_contrasts = contrasts)

  expect_true(file.exists(Configfile))
  expect_true(file.exists(samplefile))

  #check if the sample fileis  equal to the expected output
  samples_correct <- read.table(testthat::test_path("samplefile.tsv"))
  samples_df <- read.table(samplefile)
  expect_true(all.equal(samples_correct, samples_df))
})

testthat::test_path("config.yaml")

scan(testthat::test_path("config.yaml"), what = '', sep='\n')
