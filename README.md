## `AnanseSeurat` package

### Export single cell RNA&ATACseq data from a seurat object for ANANSE GRN analysis


The `AnanseSeurat` package takes pre-processed clustered single cell objects of scRNAseq and scATACseq or a multiome combination, and generates files for gene regulatory network (GRN) analysis.

### Installation

```{r eval=FALSE}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
Sys.unsetenv("GITHUB_PAT")
remotes::install_github("JGASmits/AnanseSeurat@main")
```


### Usage

First, load the package, and load your pre-processed single cell object.

```{r eval=FALSE}
library("AnanseSeurat")
rds_file <- './scANANSE/preprocessed_PDMC.Rds'
pbmc <- readRDS(rds_file)
```

Next you can output the data from your single cell object, the file format, config file and sample file are all ready to automate anansnake GRN analysis.
https://github.com/vanheeringen-lab/anansnake

```{r eval=FALSE}
export_CPM_scANANSE(pbmc,
                    min_cells <- 25,
                    output_dir ='./scANANSE/analysis',
                    cluster_id = 'predicted.id',
                    RNA_count_assay = 'RNA')

export_ATAC_scANANSE(pbmc,
                     min_cells <- 25,
                     output_dir ='./scANANSE/analysis',
                    cluster_id = 'predicted.id',
                    ATAC_peak_assay= 'peaks')

# Specify additional contrasts:
contrasts <-  c('B-naive_B-memory',
                   'B-memory_B-naive',
                   'B-naive_CD14-Mono',
                   'CD14-Mono_B-naive')

config_scANANSE(pbmc,
                min_cells <- 25,
                output_dir ='./scANANSE/analysis',
                cluster_id = 'predicted.id',
                additional_contrasts = contrasts)

DEGS_scANANSE(pbmc,
              min_cells <- 25,
              output_dir ='./scANANSE/analysis',
              cluster_id = 'predicted.id',
              additional_contrasts = contrasts)

```

### install and run anansnake 

Follow the instructions its respective github page, https://github.com/vanheeringen-lab/anansnake
Next automatically use the generated files to run GRN analysis using your single cell cluster data:


```{bash eval=FALSE}
snakemake --use-conda --conda-frontend mamba \
--configfile scANANSE/analysis/config.yaml \
--snakefile scANANSE/anansnake/Snakefile \
--resources mem_mb=48_000 --cores 12
```

### import ANANSE results back to your single cell object


```{r eval=FALSE}
ananse_result_list <- import_seurat_scANANSE(pbmc,
                        cluster_id = 'predicted.id',
                        anansnake_inf_dir= "./scANANSE/analysis/influence")

pbmc <- ananse_result_list[1]
TF_influence_scores <- ananse_result_list[2]
```


### Thanks to:

* Julian Arts and his Pycharm equivalent of this package https://github.com/Arts-of-coding/AnanseScanpy 
* Siebren Frohlich and his anansnake implementation https://github.com/vanheeringen-lab/anansnake
* Branco Heuts for testing
