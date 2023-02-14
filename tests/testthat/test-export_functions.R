#test_functions

test_that("CPM & Count files are sucessfully generated", {
  outdir <- paste0(tempdir(),'/CPM')
  TPMfile <- paste0(outdir, '/TPM.tsv')
  Countfile <- paste0(outdir, '/RNA_Counts.tsv')
  on.exit(unlink(c(Countfile, TPMfile, outdir)))
  
  expect_false(file.exists(TPMfile))
  expect_false(file.exists(Countfile))
  
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))
  
  export_CPM_scANANSE(
    sce,
    min_cells <- 1,
    output_dir = outdir,
    cluster_id = 'seurat_clusters',
    RNA_count_assay = 'RNA'
  )
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
  outdir <- paste0(tempdir(),'/peaks')
  Peakfile <- paste0(outdir, '/Peak_Counts.tsv')
  on.exit(unlink(c(Peakfile, outdir)))
  
  expect_false(file.exists(Peakfile))
  
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))
  
  export_ATAC_scANANSE(
    sce,
    min_cells <- 1,
    output_dir = outdir,
    cluster_id = 'seurat_clusters',
    ATAC_peak_assay = 'peaks'
  )
  
  expect_true(file.exists(Peakfile))
  
  #check if the files are equal to the expected output
  Peak_correct <- read.table(testthat::test_path("Peak_Counts.tsv"))
  Peak_df <- read.table(Peakfile)
  expect_true(all.equal(Peak_df, Peak_correct))
})

test_that("Config and sample file are sucessfully generated", {
  outdir <- paste0(tempdir(),'/config')
  Configfile <- paste0(outdir, '/config.yaml')
  samplefile <- paste0(outdir, '/samplefile.tsv')
  
  on.exit(unlink(c(Configfile, samplefile, outdir)))
  
  expect_false(file.exists(Configfile))
  expect_false(file.exists(samplefile))
  
  
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))
  
  contrasts <- list('cluster1_cluster2',
                   'cluster2_cluster1',
                   'cluster3_cluster2')
  config_scANANSE(
    sce,
    min_cells <- 1,
    output_dir = outdir,
    cluster_id = 'seurat_clusters',
    genome = paste0(outdir, '/hg38'),
    additional_contrasts = contrasts
  )
  
  expect_true(file.exists(Configfile))
  expect_true(file.exists(samplefile))
  
  #check if the sample fileis  equal to the expected output
  samples_correct <-
    read.table(testthat::test_path("samplefile.tsv"))
  samples_df <- read.table(samplefile)
  expect_true(all.equal(samples_correct, samples_df))
})

test_that("Maelstrom count matrix is sucessfully generated", {
  outdir <- paste0(tempdir(),'/maelstrom')
  peakfile <- paste0(outdir, '/Peaks_scaled.tsv')
  on.exit(unlink(c(peakfile,outdir)))
  expect_false(file.exists(peakfile))
  
  #load dummy single cell object
  sce <- readRDS(testthat::test_path("sce_obj.Rds"))
  
  export_ATAC_maelstrom(
    sce,
    min_cells <- 1,
    output_dir = outdir,
    cluster_id = 'seurat_clusters',
    ATAC_peak_assay = 'peaks'
  )
  
  expect_true(file.exists(peakfile))
  
  #check if the sample fileis  equal to the expected output
  peaks_correct <-
    read.table(testthat::test_path("Peaks_scaled.tsv"))
  peaks_df <- read.table(peakfile)
  expect_true(all.equal(peaks_correct, peaks_df))
})
