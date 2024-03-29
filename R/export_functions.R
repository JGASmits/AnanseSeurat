#' export_CPM_scANANSE
#'
#' This functions exports CPM values from a seurat object
#' @param seurat_object the seurat object used to export the CPM values from
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param RNA_count_assay assay of the seurat object containing the RNA count data
#' @param cluster_id ID used for finding clusters of cells
#' @return None, outputs CPM and counts files in the output directory
#' @examples
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' export_CPM_scANANSE(sce_small, min_cells = 2, output_dir = tempdir())
#' @export
export_CPM_scANANSE <- function(seurat_object,
                                output_dir,
                                min_cells = 50,
                                RNA_count_assay = "RNA",
                                cluster_id = 'seurat_clusters') {
  if (missing(output_dir)) {
    stop('no output_dir specified')
  }
  dir.create(file.path(output_dir), showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  message('calculating CPM')
  seurat_object <- Seurat::NormalizeData(
    seurat_object,
    assay  = RNA_count_assay,
    normalization.method = 'RC',
    scale.factor = 1e6
  )
  
  rna_count_lists <- list()
  FPKM_count_lists <- list()
  cluster_names <- list()
  i <- 1
  for (cluster in levels(Seurat::Idents(seurat_object))) {
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells) {
      message(paste0('gather data from ', cluster, ' with ', n_cells, ' cells'))
      cluster_names[i] <- cluster
      
      #lets grab the scRNAseq counts
      mat_counts <-
        seurat_object_cluster@assays[[RNA_count_assay]]@counts
      mat_counts <- as.matrix(rowSums(as.matrix(mat_counts)))
      colnames(mat_counts) <- cluster
      rna_count_lists[i] <- as.data.frame(mat_counts)
      
      #lets also grab the normalized intensities, wich we will
      #use as a FPKM like matrix, since the reads are UMI normalized
      #we do not expect a gene-length bias
      mat_norm <-
        seurat_object_cluster@assays[[RNA_count_assay]]@data
      mat_norm <- as.matrix(rowSums(as.matrix(mat_norm)))
      colnames(mat_norm) <- cluster
      FPKM_count_lists[i] <- as.data.frame(mat_norm)
      i <- i + 1
    } else{
      warning(paste0(
        cluster,
        ' less than ',
        min_cells,
        ' not including this cluster'
      ))
    }
  }
  
  #RNA count matrix
  count_matrix <- do.call(rbind, rna_count_lists)
  count_matrix <- as.data.frame(count_matrix)
  count_matrix <- t(count_matrix)
  colnames(count_matrix) <- cluster_names
  rownames(count_matrix) <- rownames(mat_counts)
  count_matrix <- as.data.frame(count_matrix)
  count_matrix$average <- round(rowMeans(count_matrix))
  
  #FPKM matrix
  FPKM_matrix <- do.call(rbind, FPKM_count_lists)
  FPKM_matrix <- as.data.frame(FPKM_matrix)
  FPKM_matrix <- t(FPKM_matrix)
  colnames(FPKM_matrix) <- cluster_names
  rownames(FPKM_matrix) <- rownames(mat_norm)
  FPKM_matrix <- as.data.frame(FPKM_matrix)
  FPKM_matrix$average <- rowMeans(FPKM_matrix)
  
  count_file <- paste(output_dir, "RNA_Counts.tsv", sep = '/')
  CPM_file <- paste(output_dir, "TPM.tsv", sep = '/')
  
  utils::write.table(as.matrix(count_matrix), count_file, sep = '\t', quote = FALSE)
  utils::write.table(as.matrix(FPKM_matrix), CPM_file, sep = '\t', quote = FALSE)
  return('done exporting cluster data files')
}
#' export_ATAC_scANANSE
#'
#' This functions exports ATAC values from a seurat object
#' @param seurat_object object
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param ATAC_peak_assay assay of the seurat object containing the peaks and peakcounts
#' @param cluster_id ID used for finding clusters of cells
#' @return None, outputs ATAC peak count file in the output directory
#' @examples
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' export_ATAC_scANANSE(sce_small, min_cells = 2, output_dir = tempdir())
#' @export
export_ATAC_scANANSE <- function(seurat_object,
                                 output_dir,
                                 min_cells = 50,
                                 ATAC_peak_assay = "peaks",
                                 cluster_id = 'seurat_clusters') {
  if (missing(output_dir)) {
    stop('no output_dir specified')
  }
  dir.create(file.path(output_dir), showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  
  peak_count_lists <- list()
  cluster_names <- list()
  i <- 1
  for (cluster in levels(Seurat::Idents(seurat_object))) {
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells) {
      message(paste0('gather data from ', cluster, ' with ', n_cells, ' cells'))
      cluster_names[i] <- cluster
      
      #lets grab the ATAC signal
      mat <- seurat_object_cluster@assays[[ATAC_peak_assay]]@counts
      mat <- as.matrix(rowSums(as.matrix(mat)))
      colnames(mat) <- cluster
      peak_count_lists[i] <- as.data.frame(mat)
      i <- i + 1
    } else{
      warning(paste0(
        cluster,
        ' less than ',
        min_cells,
        ' not including this cluster'
      ))
    }
  }
  
  #ATAC peak matrix
  activity_matrix <- do.call(rbind, peak_count_lists)
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix <- t(activity_matrix)
  colnames(activity_matrix) <- cluster_names
  #output the peaks as a bed file:
  peaks <- rownames(mat)
  peaks <- stringr::str_split_fixed(peaks, '-', 3)
  peaknames <- paste0(peaks[, 1], ":", peaks[, 2], '-', peaks[, 3])
  rownames(activity_matrix) <- peaknames
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix$average <- round(rowMeans(activity_matrix))
  Peak_file <- paste(output_dir, "Peak_Counts.tsv", sep = '/')
  utils::write.table(as.matrix(activity_matrix),
                     Peak_file,
                     sep = '\t',
                     quote = FALSE)
}


#' config_scANANSE
#'
#' This functions generates a sample file and config file for running Anansnake based on the seurat object
#' @param seurat_object seurat object
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param genome genomepy name or location of the genome fastq file
#' @param cluster_id ID used for finding clusters of cells
#' @param additional_contrasts additional contrasts to add between clusters within cluster_ID
#' @return None, outputs snakemake config file in the output directory
#' @examples
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' config_scANANSE(sce_small, min_cells = 2, output_dir = tempdir())
#' @export
config_scANANSE <- function(seurat_object,
                            output_dir,
                            min_cells = 50,
                            cluster_id = 'seurat_clusters',
                            genome = './scANANSE/data/hg38',
                            additional_contrasts = c()) {
  if (missing(output_dir)) {
    stop('no output_dir specified')
  }
  dir.create(file.path(output_dir), showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  Peak_file <- paste(output_dir, "Peak_Counts.tsv", sep = '/')
  count_file <- paste(output_dir, "RNA_Counts.tsv", sep = '/')
  CPM_file <- paste(output_dir, "TPM.tsv", sep = '/')
  cluster_names <- list()
  i <- 1
  for (cluster in levels(Seurat::Idents(seurat_object))) {
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells) {
      cluster_names[i] <- cluster
      i <- i + 1
    }
  }
  
  #lets generate the snakemake sample file
  sample_names <- c(cluster_names, 'average')
  sample_file_df <- as.data.frame(sample_names)
  sample_file_df <- as.data.frame(t(sample_file_df))
  colnames(sample_file_df) <- 'sample'
  sample_file_df$assembly <- basename(genome)
  sample_file_df$anansesnake <- sample_file_df$sample
  sample_file_location <- paste0(output_dir, "/samplefile.tsv")
  
  utils::write.table(
    sample_file_df,
    sample_file_location,
    sep = '\t',
    row.names = FALSE,
    quote = FALSE
  )
  #lets generate the snakemake config file
  contrast_list <-
    paste0('anansesnake_', cluster_names, '_average')
  if (length(additional_contrasts)>0){
    message('adding additional contrasts')
    additional_contrasts <- 
      paste0('anansesnake_', additional_contrasts)
    contrast_list <- c(contrast_list, additional_contrasts)
  }
  
  file <- paste0(output_dir, "/config.yaml")
  
  lines <- c(
    paste0("result_dir: ", output_dir, '\n'),
    paste0("rna_samples: ", sample_file_location, '\n'),
    paste0("rna_tpms: ", CPM_file, '\n'),
    paste0("rna_counts: ", count_file, '\n'),
    paste0("atac_samples: ", sample_file_location, '\n'),
    paste0("atac_counts: ", Peak_file, '\n'),
    paste0("genome: ", genome, "\n"),
    "database: gimme.vertebrate.v5.0 \n",
    "jaccard: 0.1 \n",
    "edges: 500_000 \n",
    "padj: 0.05 \n",
    'plot_type: "png" \n',
    "get_orthologs: false \n",
    "contrasts: \n"
  )
  
  string <- paste(lines, collapse = '')
  cat(string,
      file = paste0(output_dir, "/config.yaml"),
      append = FALSE)
  
  for (contr in contrast_list) {
    cat(
      paste0('  - "', contr, '"'),
      "\n",
      file = paste0(output_dir, "/config.yaml"),
      append = TRUE
    )
  }
}

#' export_seurat_Maelstrom
#'
#' normalize and export the peak table of a seurat object based on clusters
#' @param seurat_object object
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param ATAC_peak_assay assay of the seurat object containing the peaks and peakcounts
#' @param cluster_id ID used for finding clusters of cells
#' @param select_top_rows only output the top variable rows, or all rows if false
#' @param n_top_rows amount of variable rows to export
#' @return None, outputs maelstrom peak counts table in the output directory
#' @examples
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' config_scANANSE(sce_small, min_cells = 2, output_dir = tempdir())
#' @export
export_ATAC_maelstrom <- function(seurat_object,
                                  output_dir,
                                  min_cells = 50,
                                  ATAC_peak_assay = "peaks",
                                  cluster_id = 'seurat_clusters',
                                  select_top_rows = TRUE,
                                  n_top_rows = 100000) {
  if (missing(output_dir)) {
    stop('no output_dir specified')
  }
  dir.create(file.path(output_dir), showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  peak_count_lists <- list()
  cluster_names <- list()
  i <- 1
  for (cluster in levels(Seurat::Idents(seurat_object))) {
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells) {
      message(paste0('gather data from ', cluster, ' with ', n_cells, ' cells'))
      cluster_names[i] <- cluster
      #lets grab the ATAC signal
      mat <- seurat_object_cluster@assays[[ATAC_peak_assay]]@counts
      mat <- as.matrix(rowSums(as.matrix(mat)))
      colnames(mat) <- cluster
      peak_count_lists[i] <- as.data.frame(mat)
      i <- i + 1
    } else{
      message(paste0(cluster, ' less than ', min_cells, ' skipping'))
    }
  }
  
  #ATAC peak matrix
  activity_matrix <- do.call(rbind, peak_count_lists)
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix <- t(activity_matrix)
  colnames(activity_matrix) <- cluster_names
  #output the peaks as a bed file:
  peaks <- rownames(mat)
  peaks <- stringr::str_split_fixed(peaks, '-', 3)
  peaknames <- paste0(peaks[, 1], ":", peaks[, 2], '-', peaks[, 3])
  rownames(activity_matrix) <- peaknames
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix <- log2(activity_matrix + 1)
  activity_matrix <- scale(activity_matrix)
  activity_matrix <- as.data.frame(activity_matrix)
  if (select_top_rows) {
    if (nrow(activity_matrix) > n_top_rows) {
      message(paste0(
        'large dataframe detected, selecting top variable rows n = ',
        n_top_rows
      ))
      message('if entire dataframe is required, add select_top_rows = False as a parameter')
      message('or change ammount of rows via the n_top_rows parameter')
      row_variance <- apply(activity_matrix, 1, stats::var)
      activity_matrix$RowVar <- row_variance
      activity_matrix <-
        activity_matrix[order(activity_matrix$RowVar, decreasing = TRUE),]
      activity_matrix <- utils::head(activity_matrix, n_top_rows)
      activity_matrix$RowVar <- NULL
    }
  }
  Peak_file <- paste(output_dir, "Peaks_scaled.tsv", sep = '/')
  utils::write.table(activity_matrix,
                     Peak_file,
                     sep = '\t',
                     quote = FALSE)
}