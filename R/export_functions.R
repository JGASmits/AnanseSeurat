#' export_CPM_scANANSE
#'
#' This functions exports CPM values from a seurat object
#' @param seurat_object the seurat object used to export the CPM values from
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param RNA_count_assay assay of the seurat object containing the RNA count data
#' @param cluster_id ID used for finding clusters of cells
#' @return CPM and count files
#' @export
export_CPM_scANANSE <- function(seurat_object,
                                min_cells = 50,
                                output_dir = '~/',
                                RNA_count_assay = "RNA",
                                cluster_id = 'seurat_clusters'
) {
  dir.create(file.path(output_dir),showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  print('calculating CPM')
  seurat_object <- Seurat::NormalizeData(seurat_object,
                                         assay  = RNA_count_assay,
                                         normalization.method = 'RC',
                                         scale.factor = 1e6
  )
  rna_count_lists = list()
  FPKM_count_lists = list()
  cluster_names = list()
  i = 1
  for (cluster in levels(Seurat::Idents(seurat_object))){
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells){
      print(paste0('gather data from ', cluster, ' with ', n_cells, ' cells'))
      cluster_names[i] <- cluster

      #lets grab the scRNAseq counts
      mat_counts <- seurat_object_cluster@assays[[RNA_count_assay]]@counts
      mat_counts <- as.matrix(rowSums(as.matrix(mat_counts)))
      colnames(mat_counts) <- cluster
      rna_count_lists[i] <- as.data.frame(mat_counts)

      #lets also grab the normalized intensities, wich we will
      #use as a FPKM like matrix, since the reads are UMI normalized
      #we do not expect a gene-length bias
      mat_norm <- seurat_object_cluster@assays[[RNA_count_assay]]@data
      mat_norm <- as.matrix(rowSums(as.matrix(mat_norm)))
      colnames(mat_norm) <- cluster
      FPKM_count_lists[i] <- as.data.frame(mat_norm)
      i = i + 1
    }else{print(paste0(cluster,' less than ',min_cells,' skipping'))}
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
  CPM_file <-paste(output_dir, "TPM.tsv", sep = '/')

  utils::write.table(count_matrix, count_file, sep = '\t', quote = F)
  utils::write.table(FPKM_matrix, CPM_file, sep = '\t', quote = F)
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
#' @export
export_ATAC_scANANSE <- function(seurat_object,
                                 min_cells = 50,
                                 output_dir = '~/',
                                 ATAC_peak_assay = "peaks",
                                 cluster_id = 'seurat_clusters'
) {
  dir.create(file.path(output_dir),showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id

  peak_count_lists = list()
  cluster_names = list()
  i = 1
  for (cluster in levels(Seurat::Idents(seurat_object))){
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells){
      print(paste0('gather data from ', cluster, ' with ', n_cells, ' cells'))
      cluster_names[i] <- cluster

      #lets grab the ATAC signal
      mat <- seurat_object_cluster@assays[[ATAC_peak_assay]]@counts
      mat <- as.matrix(rowSums(as.matrix(mat)))
      colnames(mat) <- cluster
      peak_count_lists[i] <- as.data.frame(mat)
      i = i + 1
    }else{print(paste0(cluster,' less than ',min_cells,' skipping'))}
  }

  #ATAC peak matrix
  activity_matrix <- do.call(rbind, peak_count_lists)
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix <- t(activity_matrix)
  colnames(activity_matrix) <- cluster_names
  #output the peaks as a bed file:
  peaks <- rownames(mat)
  peaks <- stringr::str_split_fixed(peaks, '-', 3)
  peaknames <- paste0(peaks[,1],":",peaks[,2],'-',peaks[,3])
  rownames(activity_matrix) <- peaknames
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix$average <- round(rowMeans(activity_matrix))

  Peak_file <-paste(output_dir, "Peak_Counts.tsv", sep = '/')
  utils::write.table(activity_matrix, Peak_file, sep = '\t', quote = F)
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
#' @export
config_scANANSE <- function(seurat_object,
                            min_cells = 50,
                            output_dir = '~/',
                            cluster_id = 'seurat_clusters',
                            genome = './scANANSE/data/hg38',
                            additional_contrasts = 'None'
) {
  dir.create(file.path(output_dir), showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  Peak_file <-paste(output_dir, "Peak_Counts.tsv", sep = '/')
  count_file <- paste(output_dir, "RNA_Counts.tsv", sep = '/')
  CPM_file <-paste(output_dir, "TPM.tsv", sep = '/')
  cluster_names = list()
  i = 1
  for (cluster in levels(Seurat::Idents(seurat_object))){
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells){
      cluster_names[i] <- cluster
      i = i +1
    }}

  #lets generate the snakemake sample file
  sample_names <- c(cluster_names, 'average')
  sample_file_df <- as.data.frame(sample_names)
  sample_file_df <- as.data.frame(t(sample_file_df))
  colnames(sample_file_df) <- 'sample'
  sample_file_df$assembly <- basename(genome)
  sample_file_df$anansesnake <- sample_file_df$sample
  sample_file_location <- paste0(output_dir,"/samplefile.tsv")

  utils::write.table(sample_file_df,sample_file_location,
                     sep='\t',
                     row.names = FALSE,
                     quote = F)
  #lets generate the snakemake config file
  contrast_list <- as.list(paste0('anansesnake_',cluster_names,'_average'))
  if (typeof(additional_contrasts) == 'list'){print('adding additional contrasts')
    additional_contrasts <- paste0('anansesnake_',additional_contrasts)
    contrast_list = c(contrast_list,additional_contrasts)}

  file=paste0(output_dir,"/config.yaml")

  lines= c(
    paste0("result_dir: ",output_dir,'\n'),
    paste0("rna_samples: ",sample_file_location,'\n'),
    paste0("rna_tpms: ",CPM_file,'\n'),
    paste0("rna_counts: ",count_file,'\n'),
    paste0("atac_samples: ",sample_file_location,'\n'),
    paste0("atac_counts: ",Peak_file,'\n'),
    paste0("genome: ",genome, "\n"),
    "database: gimme.vertebrate.v5.0 \n",
    "jaccard: 0.1 \n",
    "edges: 500_000 \n",
    "padj: 0.05 \n",
    'plot_type: "png" \n',
    "get_orthologs: false \n",
    "contrasts: \n")

  string= paste(lines, collapse = '')
  cat(string,file=paste0(output_dir,"/config.yaml"),append = F)

  for (contr in contrast_list){
    cat(paste0('  - "',contr, '"'),"\n",file=paste0(output_dir,"/config.yaml"),append = T)
    print(contr)
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
#' @export
export_ATAC_maelstrom <- function(seurat_object,
                                  min_cells = 50,
                                  output_dir = '~/',
                                  ATAC_peak_assay = "peaks",
                                  cluster_id = 'seurat_clusters',
                                  select_top_rows = T,
                                  n_top_rows = 100000

) {
  dir.create(file.path(output_dir),showWarnings = FALSE)
  Seurat::Idents(seurat_object) <- cluster_id
  peak_count_lists = list()
  cluster_names = list()
  i = 1
  for (cluster in levels(Seurat::Idents(seurat_object))){
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells){
      print(paste0('gather data from ', cluster, ' with ', n_cells, ' cells'))
      cluster_names[i] <- cluster
      #lets grab the ATAC signal
      mat <- seurat_object_cluster@assays[[ATAC_peak_assay]]@counts
      mat <- as.matrix(rowSums(as.matrix(mat)))
      colnames(mat) <- cluster
      peak_count_lists[i] <- as.data.frame(mat)
      i = i + 1
    }else{print(paste0(cluster,' less than ',min_cells,' skipping'))}
  }

  #ATAC peak matrix
  activity_matrix <- do.call(rbind, peak_count_lists)
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix <- t(activity_matrix)
  colnames(activity_matrix) <- cluster_names
  #output the peaks as a bed file:
  peaks <- rownames(mat)
  peaks <- stringr::str_split_fixed(peaks, '-', 3)
  peaknames <- paste0(peaks[,1],":",peaks[,2],'-',peaks[,3])
  rownames(activity_matrix) <- peaknames
  activity_matrix <- as.data.frame(activity_matrix)
  activity_matrix <- log2(activity_matrix + 1)
  activity_matrix <- scale(activity_matrix)
  activity_matrix <- as.data.frame(activity_matrix)
  print(nrow(activity_matrix))
  if (select_top_rows){
    if (nrow(activity_matrix)>n_top_rows){
      print(paste0('large dataframe detected, selecting top variable rows n = ',n_top_rows))
      print('if entire dataframe is required, add select_top_rows = False as a parameter')
      print('or change ammount of rows via the n_top_rows parameter')
      row_variance <- apply(activity_matrix, 1, stats::var)
      activity_matrix$RowVar <-row_variance
      activity_matrix <-activity_matrix[order(activity_matrix$RowVar, decreasing = T),]
      activity_matrix <- utils::head(activity_matrix,n_top_rows)
      activity_matrix$RowVar <- NULL
    }}
  Peak_file <-paste(output_dir, "Peaks_scaled.tsv", sep = '/')
  utils::write.table(activity_matrix, Peak_file, sep = '\t', quote = F)
}

#' import_seurat_Maelstrom
#'
#' load Maelstrom enriched motifs
#' @param seurat_object object
#' @param cluster_id ID used for finding clusters of cells
#' @param maelstrom_dir directory where the maelstrom results are stored
#' @param return_df return both the seurat object and a dataframe with maelstrom scores as a list
#' @export
import_seurat_maelstrom <- function(seurat_object,
                                    cluster_id = 'seurat_clusters',
                                    maelstrom_dir = '~/',
                                    return_df = F

){
  maelstrom_file <- paste0(maelstrom_dir,"final.out.txt")
  maelstrom_df <- utils::read.table(maelstrom_file, header = T,row.names=1, sep = '\t',check.names=FALSE)
  rownames(maelstrom_df) <- gsub('_', '-',rownames(maelstrom_df))
  maelstrom_Zscore <- maelstrom_df[grep("z-score ", colnames(maelstrom_df))]
  maelstrom_corr <- maelstrom_df[grep("corr ", colnames(maelstrom_df))]
  #add the motif intensities to the seurat object
  motif_df <- t(maelstrom_Zscore)
  cell_signal_list <- list()
  rownames(motif_df) <- mapply(gsub, pattern="z-score ", x=rownames(motif_df), replacement='')
  i = 1

  #get cell IDs from the single cell object
  for (cell_type in rownames(motif_df)){
    #lets check if any cells with this annotation are present in the seurat object
    print(cell_type)
    cells_in_seurat <- seurat_object[[cluster_id]] ==cell_type
    if (TRUE %in% cells_in_seurat){
      relevant_IDS <- seurat_object[,seurat_object[[cluster_id]] ==cell_type][[cluster_id]]
      TF_signal <- motif_df[rownames(motif_df) == cell_type,]
      TF_signal <- as.data.frame(TF_signal)
      TF_signal <- t(TF_signal)
      TF_signal <- TF_signal[rep(seq_len(nrow(TF_signal)), each = nrow(relevant_IDS)), ]
      rownames(TF_signal) <- rownames(as.data.frame(relevant_IDS))
      #print(TF_signal)
      cell_signal_list[[i]] <- TF_signal}else{print(paste0(paste0('no cells of id ',cell_type),' found in seurat object'))}
    i = i + 1
  }
  TF_array <- do.call("rbind", cell_signal_list)
  #print(head(TF_array))
  #add cell IDs that do not have Maelstrom signal and give them a NA vallue
  cell_barcodes_all <- rownames(as.data.frame(Seurat::Idents(object = seurat_object)))
  cell_barcodes_missing <- cell_barcodes_all[!cell_barcodes_all %in% rownames(TF_array)]
  if (length(cell_barcodes_missing) > 0){
    cell_barcodes_missing <- as.data.frame(cell_barcodes_missing)
    for (TF in colnames(TF_signal)){
      cell_barcodes_missing[[TF]] <- NA
    }
    rownames(cell_barcodes_missing) <- cell_barcodes_missing$cell_barcodes_missing
    cell_barcodes_missing$cell_barcodes_missing <- NULL
    TF_array <- rbind(TF_array,cell_barcodes_missing)
  }
  seurat_object[['maelstrom']] <- Seurat::CreateAssayObject(t(TF_array))
  motif_df <- t(motif_df)
  if (return_df){return(list(seurat_object,motif_df))}
  else{return(seurat_object)}
}
