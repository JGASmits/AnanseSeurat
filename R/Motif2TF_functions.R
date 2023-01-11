#' Maelstrom_Motif2TF
#'
#' create motif-factor links & export tables for printing motif score alongside its binding factor
#' @param seurat_object object
#' @param mot_mat motif_matrix, if not provided extracts one from the single cell object from the maelstrom assay
#' @param m2f_df motif to factor dataframe, if not provided extracts from the maelstrom directory
#' @param cluster_id ID used for finding clusters of cells
#' @param maelstrom_dir directory where the GimmeMotifs m2f table is stored
#' @param RNA_expression_assay Seurat assay containing factor expression info
#' @param RNA_expression_slot slot within assay used for calculating mean factor expression per cluster
#' @param expr_tresh minimum sum of gene counts over all cells in RNA_expression_assay to filter genes by
#' @param cor_tresh minimum value of to filter the cor() output by
#' @param cor_method specify one of the cor() methods
#' @param curated_motifs use only curated motifs (T), or all motifs in the database (F)
#' @param output_dir where to store the results
#' @param combine_motifs means (take mean multiple motifscores), max_var (take motif with highest variance), or max_cor (take motif with best correlation to gene expression)
#' @param return_df return both the seurat object and two dataframes with maelstrom scores and expression values as a list
#' @return seurat object with two assays added, MotifTFcor for TFs with positive correlation to the linked motif, and MotifTFanticor for TFs with positive correlation to the linked motif
#' @examples 
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' maelstrom_dir_path <- system.file("extdata","maelstrom",package = 'AnanseSeurat') 
#' sce_small <- Maelstrom_Motif2TF(sce_small, maelstrom_dir = maelstrom_dir_path)
#' @export
#' @importFrom rlang .data
Maelstrom_Motif2TF <- function(seurat_object,
                               mot_mat = NULL,
                               m2f_df = NULL,
                               cluster_id = 'seurat_clusters',
                               maelstrom_dir = './maelstrom/',
                               combine_motifs = 'means',
                               RNA_expression_assay = "RNA",
                               RNA_expression_slot = "data",
                               expr_tresh = 10,
                               cor_tresh = 0.30,
                               curated_motifs = F,
                               cor_method = "pearson",
                               output_dir = "./",
                               return_df = F
){
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

  if((combine_motifs != 'means')&(combine_motifs != 'max_cor')&(combine_motifs != 'max_var')){
      stop("use either 'max_var', 'max_cor' or 'mav_var' as a selection method to select the most relevant motifs per TF")
  }

  if (curated_motifs){
    print('using only curated motifs from database')
    m2f_df <- m2f_df[m2f_df$Curated == 'Y',]}

  ## Load needed objects
  if(is.null(mot_mat)){
    print(paste0('loading maelstrom values from maelstrom assay using the cluster identifier ',cluster_id))
    mot_mat <- per_cluster_df(seurat_object,
                               assay = 'maelstrom',
                               cluster_id = cluster_id)
    mot_mat <- as.matrix(mot_mat)
  }

  ## Set up Seurat object
  print("non-expressed genes are removed")
  Seurat::DefaultAssay(seurat_object) <- RNA_expression_assay
  genes_expressed <- rownames(seurat_object)[rowSums(as.matrix(seurat_object[[RNA_expression_assay]]@counts)) >= expr_tresh]
  seurat_object[[RNA_expression_assay]] <- Seurat::CreateAssayObject(counts = seurat_object[[RNA_expression_assay]]@data[genes_expressed,])

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
  exp_mat <- exp_mat[!rowSums(as.matrix(exp_mat)) <= 0,]
  ## Select the same exp_mat columns as in mot_mat columns (if the grouping var is the same)
  exp_mat <- exp_mat[,colnames(mot_mat)]

  ## limit table to motifs and TFs present in dataset
  mot_mat <- mot_mat[rownames(mot_mat) %in% m2f_df$Motif,]
  TF_mat <- exp_mat[rownames(exp_mat) %in% m2f_df$Factor,]

  m2f_df_match <- m2f_df[m2f_df$Motif %in% rownames(mot_mat) & m2f_df$Factor %in% rownames(TF_mat),]
  ## Creat unique motif-factor entries, losing curation information
  m2f_df_match <- unique(m2f_df_match[,!colnames(m2f_df_match) %in% c("Evidence", "Curated")])

  ## perform correlations between cluster expression and cluster motif enrichment
  ## calculate motif score variation over clusters
  m2f_df_match$cor <- NA
  m2f_df_match$var <- NA
  for (i in 1:nrow(m2f_df_match)){
    m2f_df_match$cor[i] <- stats::cor(mot_mat[m2f_df_match$Motif[i],],exp_mat[m2f_df_match$Factor[i],], method = cor_method)
    m2f_df_match$var[i] <- stats::var(mot_mat[m2f_df_match$Motif[i],])
  }

  ## Only keep motif-TF combinations with an absolute R higher than treshold
  print(paste0("Only keep motif-TF combinations with an R > ", cor_tresh))
  m2f_df_match <- m2f_df_match[base::abs(m2f_df_match$cor) > cor_tresh,]

  ## Select highest absolute correlation of TF and motif
  m2f_df_unique <- as.data.frame(m2f_df_match %>% dplyr::group_by(m2f_df_match$Motif) %>%
    dplyr::arrange(dplyr::desc(base::abs(m2f_df_match$cor))) %>% dplyr::filter(dplyr::row_number() == 1))

  #print(paste0('total length m2f_df_unique ', length(m2f_df_unique$cor)))
  #Select only positive correlations or only negative correlations (repressors)
  for (typeTF in c('MotifTFcor','MotifTFanticor')){
    m2f <- m2f_df_unique
    if (typeTF == 'MotifTFanticor'){
      print("Selecting anticorrelating TFs")
      #print(paste0('total m2f', length(m2f$cor)))
      m2f <- m2f_df_unique[m2f_df_unique$cor < 0,]
      #print(paste0('total m2f', length(m2f$cor)))
    } else {
      print("Selecting correlating TFs")
      #print(paste0('total m2f', length(m2f$cor)))
      m2f <- m2f_df_unique[m2f_df_unique$cor > 0,]
      #print(paste0('total m2f', length(m2f$cor)))
    }

    m2f$associated_motifs <- NA
    for (tf in unique(m2f$Factor)){
      ## Generate a string with all associated motifs and their correlation to the tf
      motif_vector <- c()
      for (motif in unique(m2f[m2f$Factor == tf, c("Motif")])){
        motif_cor <- paste0(c(motif, unique(m2f[m2f$Factor== tf & m2f$Motif == motif,"cor"])), collapse = ":")
        motif_vector <- paste0(c(motif_vector, motif_cor), collapse = "_")
      }
      m2f[m2f$Factor == tf,]$associated_motifs <- motif_vector
    }
    ## Order motifs according to m2f & replace with TF name
    mot_plot <- as.matrix(mot_mat[match(m2f$Motif,rownames(mot_mat)),])

    ## Make motif score per TF (selecting most variable motif per TF or make mean of all motifs associated).
    if (combine_motifs == 'means'){
      rownames(mot_plot) <- m2f$Factor
      print("Take mean motif score of all binding motifs per TF")
      ## Take mean of motifs linked to the same TF
      mot_plot <- stats::aggregate(mot_plot, list(row.names(mot_plot)), mean)
      print(mot_plot)
      mot_plot <- as.data.frame(mot_plot, row.names = mot_plot[,1])[,-1]
      m2f <- as.data.frame(m2f[!duplicated(m2f$Factor), c("Factor","associated_motifs"), drop = FALSE])
    }
    if(combine_motifs == 'max_cor') {
      print("Motif best (absolute) correlated to expression is selected per TF")
      ## Using m2f file for selecting highest correlating motif to factor:
      m2f <- m2f[order(base::abs(m2f[,"cor"]), decreasing = T)]
      m2f <- m2f[!duplicated(m2f$Factor),c("Factor","Motif","cor"), drop = FALSE]
      mot_plot <- mot_plot[match(m2f$Motif, rownames(mot_plot)),]
      #mot_plot <- as.data.frame(mot_plot)
      rownames(mot_plot) <- m2f$Factor
    }
    if(combine_motifs == 'max_var'){
      print("Most variable binding motif is selected per TF")
      ## Using m2f file for selecting highest variable motif to factor:
      m2f <- m2f[order(base::abs(m2f[,"var"]), decreasing = T)]
      m2f <- m2f[!duplicated(m2f$Factor),]
      mot_plot <- mot_plot[match(m2f$Motif, rownames(mot_plot)),]
      rownames(mot_plot) <- m2f$Factor
    }

    ## order expression matrix and motif matrix the same way
    exp_plot <- exp_mat[match(rownames(mot_plot),rownames(exp_mat)),]

    exp_plot_scale <- t(scale(t(exp_plot)))
    mot_plot_scale <- t(scale(t(mot_plot)))

    #expression_file <- paste(output_dir, "expression_means_scaled.tsv", sep = '/')
    #utils::write.table(exp_plot_scale, expression_file, sep = '\t', quote = F)
    #motif_file <- paste(output_dir, "motif_intensities_scaled.tsv", sep = '/')
    #utils::write.table(mot_plot_scale, motif_file, sep = '\t', quote = F)
    matrix_list <- list()
    matrix_list[["expression"]] <- exp_plot
    matrix_list[["motif_score"]] <- mot_plot
    matrix_list[["scaled_expression"]] <- exp_plot_scale
    matrix_list[["scaled_motif_score"]] <- mot_plot_scale
    matrix_list[[paste0("tf2motif_selected_", combine_motifs)]] <- m2f

    ## Create seurat assay with binding factor assay
    new_assay <- as.data.frame(matrix(data = NA, ncol = length(colnames(seurat_object)), nrow = length(rownames(mot_plot))))
    colnames(new_assay) <- colnames(seurat_object)
    rownames(new_assay) <- rownames(mot_plot)
    for (cluster in colnames(mot_plot)){
      cluster_cells <- colnames(seurat_object[,seurat_object@meta.data[[cluster_id]] == cluster])
      for (TF in rownames(new_assay)){
        new_assay[TF,cluster_cells] <- as.matrix(mot_plot[TF,cluster])
      }
    }

    seurat_object[[typeTF]] <- Seurat::CreateAssayObject(new_assay)
    ## Adding meta.features with information about the motifs used in the matrix
    m2f <- as.data.frame(m2f[match(rownames(new_assay), m2f$Factor),])
    rownames(m2f) <- m2f$Factor
    seurat_object[[typeTF]]@meta.features <- m2f[,!colnames(m2f) == "Factor", drop = FALSE]
    }
  matrix_list[["seurat_object"]] <- seurat_object

  if (return_df){return(list(seurat_object,matrix_list))}
  else{return(seurat_object)}
}

#' Factor_Motif_Plot
#'
#' plot both expression of a TF, and the motif accessibility of the associated motif. Finally, fetch the motif logo from the Maelstrom directory.
#' @param seurat_object seurat object
#' @param TF_list list of TFs to plot the expression and linked motif Z-score for
#' @param assay_RNA RNA_count_assay assay containing the RNA data
#' @param assay_maelstrom maelstrom assay used for zscore vizualization, often either TFcor or TFanticor
#' @param logo_dir directory containing motif logos generated by gimme maelstrom
#' @param col colours used for zscore vizualization
#' @param dim_reduction dimensionality reduction method to use
#' @return patchwork plot containing a expression dimreduction plot, a maelstrom motif score dimreduction plot, and a png image of the motif
#' @examples 
#' sce_small <- readRDS(system.file("extdata","sce_small.Rds",package = 'AnanseSeurat'))
#' logos_dir_path <- system.file("extdata","maelstrom","logos",package = 'AnanseSeurat') 
#' sce_small <- Factor_Motif_Plot(sce_small, 
#'   c('gene1', 'gene2'), 
#'   dim_reduction = 'pca',
#'   logo_dir = logos_dir_path)
#' @export
Factor_Motif_Plot <- function(seurat_object,
                                TF_list,
                                assay_RNA = 'RNA',
                                assay_maelstrom = 'MotifTFanticor',
                                logo_dir = '~/maelstrom/logos',
                                col = c('darkred','white','darkgrey'),
                                dim_reduction = 'umap'){
  Seurat::DefaultAssay(object = seurat_object) <- assay_RNA
  plot_expression1 <- Seurat::FeaturePlot(seurat_object, features = TF_list, ncol = 1, reduction = dim_reduction)
  Seurat::DefaultAssay(object = seurat_object) <- assay_maelstrom
  plot_Maelstrom_raw <- Seurat::FeaturePlot(seurat_object,ncol = 1, features = TF_list, combine = F)
  TF_motif_table <- seurat_object@assays[[assay_maelstrom]][[]]

  #replace the TF name with the motif name for the maelstrom enrichment score
  plot_Maelstrom <- lapply(plot_Maelstrom_raw, function(x){
    TF_name <- names(x$data)[4][[1]]
    motif_name <- TF_motif_table[TF_name,]$Motif
    x = x + ggplot2::labs(title = motif_name)
    x + ggplot2::scale_colour_gradient2(low = col[1], mid = col[2], high = col[3], midpoint = 0)
    })
  plot_Maelstrom <- patchwork::wrap_plots(plot_Maelstrom , ncol = 1)
  plot_logo <- lapply(plot_Maelstrom_raw, function(x){
    TF_name <- names(x$data)[4][[1]]
    motif_name <- TF_motif_table[TF_name,]$Motif
    motif_name <- gsub('\\.','_',motif_name)
    motif_path <- paste0(logo_dir,'/',motif_name,'.png')
    logo_image <- png::readPNG(motif_path)
    ggplot2::ggplot() + ggpubr::background_image(logo_image) #+ ggplot2::coord_fixed()
    })
  plot_logo <- patchwork::wrap_plots(plot_logo, ncol = 1)
  return(plot_expression1|plot_Maelstrom|plot_logo)
}
