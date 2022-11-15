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

#' DEGS_scANANSE
#'
#' Calculate the differential genes needed for ananse influence
#' @param seurat_object seurat object
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param cluster_id ID used for finding clusters of cells
#' @param additional_contrasts additional contrasts to add between clusters within cluster_ID
#' @export
DEGS_scANANSE <- function(seurat_object,
                                      min_cells = 50,
                                      output_dir = '~/',
                                      cluster_id = 'seurat_clusters',
                                      additional_contrasts = 'None'
) {
  Seurat::Idents(seurat_object) <- cluster_id
  cluster_names = list()
  i = 1
  for (cluster in levels(Seurat::Idents(seurat_object))){
    seurat_object_cluster <- subset(x = seurat_object, idents = cluster)
    #only use ANANSE on clusters with more than 50 cells
    n_cells <- dim(seurat_object_cluster)[2]
    if (n_cells > min_cells){
      cluster_names[i] <- cluster
      i = i+1 # Increase i to add clusters iteratively to the cluster_names list
    }}

  #lets generate the snakemake config file
  contrast_list <- as.list(paste0('anansesnake_',cluster_names,'_average'))
  if (typeof(additional_contrasts) == 'list'){print('adding additional contrasts')
    additional_contrasts <- paste0('anansesnake_',additional_contrasts)
    contrast_list = c(contrast_list,additional_contrasts)}

  for (contr in contrast_list){
    print(paste0('calculating DEGS for contrast ',contr))
    comparison1 <- stringr::str_split(contr, "_")[[1]][2]
    comparison2 <- stringr::str_split(contr, "_")[[1]][3]

    DEG_file <- paste0(output_dir,'/deseq2/','hg38-anansesnake_',comparison1,'_', comparison2,'.diffexp.tsv')

    if (file.exists(DEG_file)){print('skip')} else {

      if (comparison2 == 'average'){
        DEGS <- Seurat::FindMarkers(seurat_object,
                            ident.1 = comparison1,
                            logfc.threshold = 0.1,
                            min.pct = 0.05,
                            return.thresh = 0.1)
      }else{DEGS <- Seurat::FindMarkers(seurat_object,
                                ident.1 = comparison1,
                                ident.2 = comparison2,
                                logfc.threshold = 0.1,
                                min.pct = 0.05,
                                return.thresh = 0.1)}
      DEGS <- DEGS[c('avg_log2FC','p_val_adj')]
      colnames(DEGS) <- c('log2FoldChange','padj')
      utils::write.table(DEGS,DEG_file, sep = '\t', quote = F)
    }
  }
}

#' import_seurat_scANANSE
#'
#' import the influences from a anansnake directory into a seurat object
#' @param seurat_object seurat object
#' @param cluster_id ID used for finding clusters of cells
#' @param anansnake_inf_dir influence directory generated by anansnake
#' @export
import_seurat_scANANSE <- function(seurat_object,
                                   cluster_id = 'seurat_clusters',
                                   anansnake_inf_dir = 'None'
){
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
      row_variance <- apply(activity_matrix, 1, var)
      activity_matrix$RowVar <-row_variance
      activity_matrix <-activity_matrix[order(activity_matrix$RowVar, decreasing = T),]
      activity_matrix <- head(activity_matrix,n_top_rows)
      activity_matrix$RowVar <- NULL
    }}
  Peak_file <-paste(output_dir, "Peaks_scaled.tsv", sep = '/')
  utils::write.table(activity_matrix, col.names = NA, row.names = TRUE, Peak_file, sep = '\t', quote = F)
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
export_ATAC_Maelstrom <- function(seurat_object,
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
#' @param Maelstrom_dir directory where the maelstrom results are stored
#' @export
import_seurat_Maelstrom <- function(seurat_object,
                                    cluster_id = 'seurat_clusters',
                                    Maelstrom_dir = ''
){
  maelstrom_file <- paste0(Maelstrom_dir,"final.out.txt")
  maelstrom_df <- utils::read.table(maelstrom_file, header = T,row.names=1, sep = '\t',check.names=FALSE)
  maelstrom_Zscore <- maelstrom_df[grep("z-score ", colnames(maelstrom_df))]
  maelstrom_corr <- maelstrom_df[grep("corr ", colnames(maelstrom_df))]
  #add the motif intensities to the seurat object
  motif_df <- t(maelstrom_Zscore)
  cell_signal_list <- list()
  rownames(motif_df) <- mapply(gsub, pattern="z-score ", x=rownames(motif_df), replacement='')
  i = 1
  #get cell IDs from the pbmcs
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
  #add cell IDs that do not have ANANSE signal and give them a NA vallue
  cell_barcodes_all <- rownames(as.data.frame(Seurat::Idents(object = seurat_object)))
  cell_barcodes_missing <- cell_barcodes_all[!cell_barcodes_all %in% rownames(TF_array)]
  if (length(rownames(cell_barcodes_missing)) > 0){
    cell_barcodes_missing <- as.data.frame(cell_barcodes_missing)
    for (TF in colnames(TF_signal)){
      print(TF)
      cell_barcodes_missing[[TF]] <- NA
    }
    rownames(cell_barcodes_missing) <- cell_barcodes_missing$cell_barcodes_missing
    cell_barcodes_missing$cell_barcodes_missing <- NULL
    TF_array <- rbind(TF_array,cell_barcodes_missing)
  }
  seurat_object[['Maelstrom']] <- Seurat::CreateAssayObject(t(TF_array))
  return(list(seurat_object,motif_df))
}

#' generate_TF_linked_motif_tables
#'
#' create motif-factor links & export tables for printing motif score alongside its binding factor
#' @param seurat_object object
#' @param cluster_id ID used for finding clusters of cells
#' @param Maelstrom_dir directory where the GimmeMotifs m2f table is stored
#' @param RNA_expression_assay Seurat assay containing factor expression info
#' @param RNA_expression_slot slot within assay used for calculating mean factor expression per cluster
#' @param repressors include repressing (TRUE; negative correlation) or only activating (FALSE; positive correlation) factors 
#' @param expr_tresh minimum sum of gene counts over all cells in RNA_expression_assay to filter genes by
#' @param cor_tresh minimum value of to filter the cor() output by
#' @param cor_method specify one of the cor() methods 
#' @param output_dir where to store the results
#' @param select_most_var choose highest variable motif per factor (TRUE) otherwise take the mean of all linked motifs (FALSE)
#' @export

generate_TF_Motif_tables <- function(maelstrom_object,
                                     m2f_df = NULL, 
                                     cluster_id = 'seurat_clusters',
                                     maelstrom_dir = '/ceph/rimlsfnwi/data/moldevbio/veenstra/snabel/scrna/CMdiff_10XMO_EB_d8/R_analysis/20221107_TFMotifVisualization/runmaelstrom/',
                                     RNA_expression_assay = "RNA",
                                     RNA_expression_slot = "data",
                                     repressors = FALSE,
                                     expr_tresh = 10,
                                     cor_tresh = 0.30,
                                     cor_method = "pearson",
                                     output_dir = "./",
                                     select_most_var = FALSE){
  ## Check if m2f_df object provided or path to Maelstrom output.
  ## Check if m2f_df object provided contains the right columns.
  
  m2f_df_cols <- c("Motif","Factor")
  if (is.null(m2f_df)){
    if (is.null(maelstrom_dir)){
      print("Provide path for maelstrom_dir or provide table m2f_df with Motif and Factor column.")
    } else {m2f_file <- paste0(maelstrom_dir, "/nonredundant.motifs.motif2factors.txt")
    m2f_df <- utils::read.table(m2f_file, header = T, sep = '\t',check.names=FALSE)
    }
  } else if (!unique(m2f_df_cols %in% colnames(m2f_df))){
    print("Provide m2f_df with at least 2 columns with names Motif and Factor.")
  }
  
  ## Load needed objects
  seurat_object <- maelstrom_object[[1]]
  mot_mat <- t(maelstrom_object[[2]])
  
  ## Set up Seurat object
  Seurat::DefaultAssay(seurat_object) <- RNA_expression_assay
  genes_expressed <- rownames(seurat_object)[rowSums(seurat_object[[RNA_expression_assay]]@counts) >= expr_tresh]
  seurat_object <- seurat_object[genes_expressed,]
  
  ## Select motifs with binding TFs present in object
  m2f_df <- m2f_df[m2f_df$Factor %in% rownames(seurat_object),]
  
  ## Check if data is normalized
  if (identical(seurat_object[[RNA_expression_assay]]@data, seurat_object[[RNA_expression_assay]]@counts) & RNA_expression_slot == "data"){
    print("Your data slot was not yet normalized.")
    print(paste0("Seurat NormalizeData with default settings will be run on all the genes in the ", RNA_expression_assay, " assay."))
    seurat_object <- Seurat::NormalizeData(seurat_object, assay = RNA_expression_assay)
  }
  
  ## Obtain df with mean expression
  exp_mat <- Seurat::AverageExpression(seurat_object, assays = RNA_expression_assay,
                                       slot = RNA_expression_slot, 
                                       features = m2f_df$Factor,
                                       group.by = cluster_id)[[1]]
  
  ## make sure that all genes in matrix have mean expression > 0
  exp_mat <- exp_mat[!rowSums(exp_mat) <= 0,]
  ## Select the same exp_mat columns as in mot_mat columns (if the grouping var is the same)
  exp_mat <- exp_mat[,colnames(mot_mat)]
  
  ## limit table to motifs and TFs present in dataset
  mot_mat <- mot_mat[rownames(mot_mat) %in% m2f_df$Motif,]
  TF_mat <- exp_mat[rownames(exp_mat) %in% m2f_df$Factor,]
  
  m2f_df_match <- m2f_df[m2f_df$Motif %in% rownames(mot_mat) & m2f_df$Factor %in% rownames(TF_mat),]
  
  m2f_df_match$cor <- NA
  ## perform correlations between cluster expression and cluster motif enrichment
  for (i in 1:nrow(m2f_df_match)){
    m2f_df_match$cor[i] <- stats::cor(mot_mat[m2f_df_match$Motif[i],],exp_mat[m2f_df_match$Factor[i],], method = cor_method)
  }
  
  ## Only keep motif-TF combinations with an absolute R higher than treshold
  print(paste0("Only keep motif-TF combinations with an R > ", cor_tresh))
  m2f_df_match <- m2f_df_match[base::abs(m2f_df_match$cor) > cor_tresh,]
  
  # ## Select highest correlating TF per motif
  # m2f_df_unique <- m2f_df_match %>% dplyr::group_by(m2f_df_match$Motif) %>%
  #   dplyr::arrange(desc(cor)) %>% dplyr::filter(dplyr::row_number() == 1)
  
  ## Select highest absolute correlation of TF and motif
  m2f_df_unique <- m2f_df_match %>% dplyr::group_by(m2f_df_match$Motif) %>%
    dplyr::arrange(dplyr::desc(base::abs(m2f_df_match$cor))) %>% dplyr::filter(dplyr::row_number() == 1)
  
  m2f <- m2f_df_unique
  # Select only positive correlations or only negative correlations (repressors)
  if (repressors == TRUE){
    print("Selecting repressors")
    m2f <- m2f_df_unique[m2f_df_unique$cor < 0,]
    assay_name <- "repressor"
  } else {
    print("Selecting activators")
    m2f <- m2f_df_unique[m2f_df_unique$cor > 0,]
    assay_name <- "activator" 
  }
  
  ## Order motifs according to m2f
  mot_plot <- mot_mat[match(m2f$Motif,rownames(mot_mat)),]
  ## Replace motif name by TF name
  rownames(mot_plot) <- m2f$Factor
  
  ## Make motif score per TF (selecting most variable motif per TF or make mean of all motifs associated). 
  if (select_most_var == TRUE){
    print("Most variable binding motif is selected per TF")
    mot_plot <- cbind(mot_plot, matrix(apply(mot_plot, 1, stats::var)))
    colnames(mot_plot)[length(colnames(mot_plot))]<-"var"
    mot_plot <- mot_plot[order(mot_plot[,"var"], decreasing = T),]
    # Remove duplicated occurences of rowname (ordered on variation)
    mot_plot <- mot_plot[!duplicated(rownames(mot_plot)),]
    mot_plot <- mot_plot[,!colnames(mot_plot) %in% c("var")]
  } else {
    print("Take mean motif score of all binding motifs per TF")
    ## Take mean of motifs linked to the same TF
    mot_plot <- stats::aggregate(mot_plot, list(row.names(mot_plot)), mean)
    mot_plot <- as.data.frame(mot_plot, row.names = mot_plot[,1])[,-1]
  }
  
  ## order expression matrix and motif matrix the same way
  exp_plot <- exp_mat[match(rownames(mot_plot),rownames(exp_mat)),]
  
  exp_plot_scale <- t(scale(t(exp_plot)))
  mot_plot_scale <- t(scale(t(mot_plot)))
  
  expression_file <- paste(output_dir, "expression_means_scaled.tsv", sep = '/')
  utils::write.table(exp_plot_scale, expression_file, sep = '\t', quote = F)
  motif_file <- paste(output_dir, "motif_intensities_scaled.tsv", sep = '/')
  utils::write.table(mot_plot_scale, motif_file, sep = '\t', quote = F)
  matrix_list <- list()
  matrix_list[["expression"]] <- exp_plot
  matrix_list[["motif_score"]] <- mot_plot
  matrix_list[["scaled_expression"]] <- exp_plot_scale
  matrix_list[["scaled_motif_score"]] <- mot_plot_scale
  matrix_list[["motif_tf_correlations"]] <- m2f_df_match
  
  ## Create seurat assay with binding factor assay
  new_assay <- as.data.frame(matrix(data = NA, ncol = length(colnames(seurat_object)), nrow = length(rownames(mot_plot))))
  colnames(new_assay) <- colnames(seurat_object)
  rownames(new_assay) <- rownames(mot_plot)
  for (cluster in colnames(mot_plot)){
    cluster_cells <- colnames(seurat_object[,seurat_object@meta.data[[cluster_id]] == cluster])
    for (TF in rownames(new_assay)){
      new_assay[TF,cluster_cells] <- mot_plot[TF,cluster]
    }
  }
  assay_name <- paste0(assay_name, "Score")
  seurat_object[[assay_name]] <- Seurat::CreateAssayObject(new_assay)
  matrix_list[["seurat_object"]] <- seurat_object
  
  return(matrix_list)
}
