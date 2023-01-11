#' DEGS_scANANSE
#'
#' Calculate the differential genes needed for ananse influence
#' @param seurat_object seurat object
#' @param min_cells minimum of cells a cluster needs to be exported
#' @param output_dir directory where the files are outputted
#' @param cluster_id ID used for finding clusters of cells
#' @param RNA_count_assay assay containing the RNA data
#' @param additional_contrasts additional contrasts to add between clusters within cluster_ID
#' @return None, outputs DEG files in the output directory
#' @examples 
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' DEGS_scANANSE(sce_small, min_cells = 2, output_dir = tempdir())
#' @export
DEGS_scANANSE <- function(seurat_object,
                          min_cells = 50,
                          output_dir = '~/',
                          cluster_id = 'seurat_clusters',
                          RNA_count_assay = "RNA",
                          additional_contrasts = 'None'
) {
  dir.create(file.path(paste0(output_dir,'/deseq2/')),showWarnings = FALSE)
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
                                    assay = RNA_count_assay,
                                    logfc.threshold = 0.1,
                                    min.pct = 0.05,
                                    return.thresh = 0.1)
      }else{DEGS <- Seurat::FindMarkers(seurat_object,
                                        ident.1 = comparison1,
                                        ident.2 = comparison2,
                                        assay = RNA_count_assay,
                                        logfc.threshold = 0.1,
                                        min.pct = 0.05,
                                        return.thresh = 0.1)}
      DEGS <- DEGS[c('avg_log2FC','p_val_adj')]
      colnames(DEGS) <- c('log2FoldChange','padj')
      utils::write.table(DEGS,DEG_file, sep = '\t', quote = F)
    }
  }
}

