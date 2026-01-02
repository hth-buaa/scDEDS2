#' @title Extract Gene Expression and Activity Profiles Along Pseudotime
#'
#' @description
#' This function processes single-cell multi-omics data to extract temporal expression and activity profiles for transcription factors (TFs) and target genes (TGs) along developmental trajectories.
#' It performs rigorous quality filtering, maps cells to pseudotime indices, normalizes data, and returns organized temporal profiles for downstream GRN inference.
#'
#' @param interest_cell_type_Branches The output of function order_pseudotime_and_divide_branches.
#' @param interest_cell_type_data The output of function get_interest_cell_type_data
#' @param alpha_gene
#' Numeric (0-1). Threshold for gene filtering based on zero-expression rate.
#' Genes with zero-expression rate higher than this value will be removed. Default is 0.7.
#' @param alpha_cell
#' Numeric (0-1). Threshold for cell filtering based on zero-expression rate.
#' Cells with zero-expression rate higher than this value will be removed. Default is 0.9.
#' @param ncores See ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list structure organized by cell type, where for each cell type contains:
#' \itemize{
#' \item \code{Branches_TFE}: List of data frames containing temporal expression vectors for TFs in each branch
#' \item \code{Branches_TGA}: List of data frames containing temporal activity vectors for TGs in each branch
#' \item \code{Branches_TGE}: List of data frames containing temporal expression vectors for TGs in each branch
#' \item \code{Branches_n_cell}: List of cell counts for each branch
#' \item \code{Branches}: Original branch information with pseudotime mapping
#' \item \code{Branches_TFs}: List of TF names present in each branch
#' \item \code{Branches_TGs}: List of TG names present in each branch
#' }
#' All vectors are ordered by pseudotime indices within each branch.
#'
#' @details
#' The function performs the following processing steps for each cell type and branch:
#' \enumerate{
#' \item Applies iterative quality filtering based on zero-expression rates
#' \item Maps cells to pseudotime indices (k~ values) for temporal ordering
#' \item Normalizes both expression and activity matrices using Seurat's methods
#' \item Extracts and orders temporal profiles for TFs and TGs
#' \item Returns organized data structures ready for dynamic GRN inference
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Alpha parameters are hyperparameters that significantly affect results
#' \item Output vectors are ordered by pseudotime, enabling time-series analysis
#' \item Uses parallel processing in Linux for efficient handling of multiple cell types
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_Branches = base::readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.rds")
#' interest_cell_type_data = base::readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.rds")
#' alpha_gene = 0.7
#' alpha_cell = 0.9
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_genes_pseudotime_info = get_genes_pseudotime_info(
#'   interest_cell_type_Branches = interest_cell_type_Branches,
#'   interest_cell_type_data = interest_cell_type_data,
#'   alpha_gene = alpha_gene,
#'   alpha_cell = alpha_cell,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_genes_pseudotime_info, file = "./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
#' }
get_genes_pseudotime_info = function(
    interest_cell_type_Branches = interest_cell_type_Branches,
    interest_cell_type_data = interest_cell_type_data,
    alpha_gene = 0.7,
    alpha_cell = 0.9,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Obtaining Gene Expression/Activity Vectors along Pseudotime ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.2 Data Processing - Cell Grouping"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining a function to process each cell type individually.
  get_genes_pseudotime_info_process = function(Branches,
                                               interest_scRNAseq_for_GRN, interest_scATACseq_for_GRN,
                                               interest_gene_for_GRN, interest_TGs, interest_TFs,
                                               alpha_gene, alpha_cell)
  {
    # Obtaining data for each branch.
    message("Initializing lists for storing scRNAseq, scATACseq, gene_for_GRN, TGs, TFs, n_cell of each branch.")
    Branches_RNA = base::list()
    Branches_ATAC = base::list()
    Branches_gene_for_GRN = base::list()
    Branches_TGs = base::list()
    Branches_TFs = base::list()
    Branches_n_cell = base::list()

    for (i in base::seq_along(Branches)) {
      message("Obtaining gene expression matrices for branch ", i, " (retaining genes/cells with non-zero expression values and keeping genes having zero-value frequencies not exceeding ", alpha_gene, " and cells having zero-value frequencies not exceeding ", alpha_cell, ").")
      Branches_RNA[[i]] = interest_scRNAseq_for_GRN[, base::rownames(Branches[[i]])]
      for (ii in 1:5) {
        Branches_RNA[[i]] = Branches_RNA[[i]][base::rowSums(base::as.matrix(Branches_RNA[[i]]@assays$RNA$counts)) > 0, ]
        Branches_RNA[[i]] = Branches_RNA[[i]][, base::colSums(base::as.matrix(Branches_RNA[[i]]@assays$RNA$counts)) > 0]
        zero_frequency = base::rowMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")) == 0))
        Branches_RNA[[i]] = subset(Branches_RNA[[i]], features = base::names(zero_frequency[zero_frequency <= alpha_gene]))
        zero_frequency = base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "counts")) == 0))
        Branches_RNA[[i]] = subset(Branches_RNA[[i]], cells = base::names(zero_frequency[zero_frequency <= alpha_cell]))
      }

      message("Extracting gene activity matrix for branch ", i," (retaining genes/cells with non-zero activity values).")
      Branches_ATAC[[i]] = interest_scATACseq_for_GRN[base::rownames(Branches_RNA[[i]]), base::colnames(Branches_RNA[[i]])]
      for (ii in 1:5) {
        Branches_ATAC[[i]] = Branches_ATAC[[i]][base::rowSums(base::as.matrix(Branches_ATAC[[i]]@assays$ACTIVITY$counts)) > 0, ]
        Branches_ATAC[[i]] = Branches_ATAC[[i]][, base::colSums(base::as.matrix(Branches_ATAC[[i]]@assays$ACTIVITY$counts)) > 0]
        zero_frequency = base::rowMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "counts")) == 0))
        Branches_ATAC[[i]] = subset(Branches_ATAC[[i]], features = base::names(zero_frequency[zero_frequency <= alpha_gene]))
        zero_frequency = base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "counts")) == 0))
        Branches_ATAC[[i]] = subset(Branches_ATAC[[i]], cells = base::names(zero_frequency[zero_frequency <= alpha_cell]))
      }
      Branches_RNA[[i]] = interest_scRNAseq_for_GRN[base::rownames(Branches_ATAC[[i]]), base::colnames(Branches_ATAC[[i]])]

      message("Retrieving TGs and TFs for branch ", i, ".")
      Branches_gene_for_GRN[[i]] = base::intersect(interest_gene_for_GRN, base::rownames(Branches_RNA[[i]]))
      Branches_TGs[[i]] = base::intersect(interest_TGs, Branches_gene_for_GRN[[i]])
      Branches_TFs[[i]] = base::intersect(interest_TFs, Branches_gene_for_GRN[[i]])

      message("Retrieving the cell number for branch ", i, ".")
      Branches_gene_for_GRN[[i]] = base::union(Branches_TGs[[i]], Branches_TFs[[i]])
      Branches_RNA[[i]] = Branches_RNA[[i]][base::intersect(base::rownames(Branches_RNA[[i]]), interest_gene_for_GRN), ]
      Branches_ATAC[[i]] = Branches_ATAC[[i]][base::rownames(Branches_RNA[[i]]), ]
      Branches_n_cell[[i]] = base::ncol(Branches_RNA[[i]])
    }

    # Mapping cells in expression and activity matrices to pseudotime point indices k~ for each branch.
    message("Mapping cells in expression and activity matrices to pseudotime point indices k~ for each branch.")
    for (i in base::seq_along(Branches)) {
      message("Creating a cell-to-pseudotime mapping for branch ", i, ".")
      Branches[[i]]["cell"] = base::rownames(Branches[[i]])
      Branches[[i]]["k~"] = 1:base::nrow(Branches[[i]])
      name_mapping = stats::setNames(Branches[[i]]$`k~`, Branches[[i]]$cell)

      message("Renaming cells in branch ", i, " to pseudotime indices for each branch.")
      base::colnames(Branches_RNA[[i]]) = name_mapping[base::colnames(Branches_RNA[[i]])]
      base::colnames(Branches_ATAC[[i]]) = name_mapping[base::colnames(Branches_ATAC[[i]])]
    }

    # Normalizing expression and activity matrices for each branch.
    message("Normalizing expression and activity matrices for each branch.")
    for (i in base::seq_along(Branches)) {
      message("Normalizing scRNA-seq matrix for branch ", i, ".")
      Branches_RNA[[i]] = Seurat::NormalizeData(Branches_RNA[[i]])

      message("Normalizing scATAC-seq matrix for branch ", i, ".")
      Branches_ATAC[[i]] = Seurat::NormalizeData(Branches_ATAC[[i]])
    }

    # For each branch, obtaining temporal expression vectors for each TF, temporal expression vectors for each TG, and temporal activity vectors for each TG.
    message("Initializing lists for storing temporal expression vectors of TFs, temporal expression vectors of TGs, temporal activity vectors of TGs across all branches.")
    Branches_TFE = base::list()
    Branches_TGA = base::list()
    Branches_TGE = base::list()
    for (i in base::seq_along(Branches)) {
      message("Obtaining temporal expression vectors for each TF in branch ", i, ".")
      Branches_TFE[[i]] =
        base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data")))[Branches_TFs[[i]], ]
      Branches_TFE[[i]] = Branches_TFE[[i]][, base::order(base::as.integer(base::colnames(Branches_TFE[[i]])))]

      message("Obtaining temporal activity vectors for each TG in branch ", i, ".")
      Branches_TGA[[i]] =
        base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(Branches_ATAC[[i]], layer = "data")))[Branches_TGs[[i]], ]
      Branches_TGA[[i]] = Branches_TGA[[i]][, base::order(base::as.integer(base::colnames(Branches_TGA[[i]])))]

      message("Obtaining temporal expression vectors for each TG in branch ", i, ".")
      Branches_TGE[[i]] =
        base::as.data.frame(SeuratObject::GetAssayData(Branches_RNA[[i]], layer = "data"))[Branches_TGs[[i]], ]
      Branches_TGE[[i]] = Branches_TGE[[i]][, base::order(base::as.integer(base::colnames(Branches_TGE[[i]])))]
    }

    return(base::list("Branches_TFE" = Branches_TFE,
                      "Branches_TGA" = Branches_TGA,
                      "Branches_TGE" = Branches_TGE,
                      "Branches_n_cell" = Branches_n_cell,
                      "Branches" = Branches,
                      "Branches_TFs" = Branches_TFs,
                      "Branches_TGs" = Branches_TGs))
  }

  ### Invoking the newly defined function in parallel.
  interest_cell_type_genes_pseudotime_info = parallel::mclapply(
    X = base::names(interest_cell_type_Branches),
    FUN = function(cell_type) {
      message("Preparing to extract pseudotime vectors of TFE, TGA, and TGE across all branches in cell type ", cell_type, ".")
      get_genes_pseudotime_info_process(
        Branches = interest_cell_type_Branches[[cell_type]],
        interest_scRNAseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scRNAseq_for_GRN"]],
        interest_scATACseq_for_GRN = interest_cell_type_data[[cell_type]][["interest_scATACseq_for_GRN"]],
        interest_gene_for_GRN = interest_cell_type_data[[cell_type]][["interest_gene_for_GRN"]],
        interest_TGs = interest_cell_type_data[[cell_type]][["interest_TGs"]],
        interest_TFs = interest_cell_type_data[[cell_type]][["interest_TFs"]],
        alpha_gene = alpha_gene,
        alpha_cell = alpha_cell
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_genes_pseudotime_info) = base::names(interest_cell_type_Branches)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Obtaining Gene Expression/Activity Vectors along Pseudotime ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_genes_pseudotime_info)
}
