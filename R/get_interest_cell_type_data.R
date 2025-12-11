#' @title Extract and Process Data for Specific Cell Types
#'
#' @description
#' This function processes single-cell multi-omics data for specified cell types by performing rigorous quality control,
#' filtering genes and cells based on expression thresholds,
#' and ensuring consistency between RNA expression and chromatin accessibility data.
#'
#' @param alpha_gene
#' Numeric (0-1). Threshold for gene filtering preliminary based on zero-expression rate.
#' Genes with zero-expression rate higher than this value will be removed. Default is 0.99.
#' It's not recommended to reset. Using default settings is advised.
#' @param alpha_cell Threshold for cell filtering preliminary based on zero-expression rate.
#' Cells with zero-expression rate higher than this value will be removed. Default is 0.99.
#' It's not recommended to reset. Using default settings is advised.
#' @param interest_cell_type Character vector. Names of cell types to extract and process (in meta.data of Seurat object).
#' @param Basic_Info The output of function get_expression_and_activity_matrix
#' @param ncores
#' Integer. Number of cores for parallel processing. Default is 1.
#' Recommended to run this function in Linux for parallel computing capability.
#' If running on Windows, ncores must be set to 1.
#'
#' @returns
#' A named list (names correspond to cell types) where each element contains:
#' \itemize{
#' \item \code{interest_scRNAseq_for_GRN}: Filtered scRNA-seq data for the cell type
#' \item \code{interest_scATACseq_for_GRN}: Filtered scATAC-seq data for the cell type
#' \item \code{interest_gene_for_GRN}: Filtered gene vector for GRN construction
#' \item \code{interest_TGs}: Filtered target genes present in the cell type
#' \item \code{interest_TFs}: Filtered transcription factors present in the cell type
#' }
#'
#' @details
#' The function performs the following operations for each specified cell type:
#' \enumerate{
#' \item Extracts expression matrix and removes genes/cells with all zero counts
#' \item Applies stringent filtering based on zero-expression rates using alpha thresholds
#' \item Extracts and aligns activity matrix with the filtered RNA-seq data
#' \item Ensures consistency between RNA and ATAC datasets
#' \item Returns processed data ready for downstream pseudotime and GRN analysis
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Requires that scATACseq_for_GRN contains a "seurat_annotations" metadata column (displays cell labels)
#' \item Temporarily changes working directory but restores upon completion
#' \item Uses parallel processing (Linux) for efficiency with multiple cell types
#' \item Alpha parameters are critical hyperparameters that may need tuning
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type = c("CD4 TEM", "cDc")
#' # interest_cell_type = "Treg"
#' Basic_Info = base::readRDS("./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.rds")
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_data = get_interest_cell_type_data(
#'   interest_cell_type = interest_cell_type,
#'   Basic_Info = Basic_Info,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_data, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.rds")
#' }
get_interest_cell_type_data = function(
    alpha_gene = 0.99,
    alpha_cell = 0.99,
    interest_cell_type = interest_cell_type,
    Basic_Info = Basic_Info,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting Data for Cell Types of Interest ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Inner function.
  get_interest_cell_type_data_process = function(alpha_gene = 0.99,
                                                 alpha_cell = 0.99,
                                                 cell_type, scRNAseq_for_GRN, scATACseq_for_GRN, gene_for_GRN, TGs, TFs)
  {
    message("The cell label is ", cell_type, ".")

    # Extracting expression matrix for cell types of interest (retaining genes/cells with non-zero expression values).
    message("Extracting expression matrix for cell types of interest (retaining genes/cells with non-zero expression values).")
    interest_scRNAseq_for_GRN = base::subset(scRNAseq_for_GRN, subset = seurat_annotations == cell_type)
    for (t in 1:3) {
      interest_scRNAseq_for_GRN =
        interest_scRNAseq_for_GRN[base::rowSums(base::as.matrix(interest_scRNAseq_for_GRN@assays$RNA$counts)) > 0, ]
      interest_scRNAseq_for_GRN =
        interest_scRNAseq_for_GRN[, base::colSums(base::as.matrix(interest_scRNAseq_for_GRN@assays$RNA$counts)) > 0]
    }

    # Retaining genes with missing value rate less than alpha_gene and cells with missing value rate less than alpha_cell.
    message("Retaining genes with missing value rate less than ", alpha_gene, " and cells with missing value rate less than ", alpha_cell, ".")
    # library(Matrix)
    for (t in 1:3) {
      zero_frequency =
        base::rowMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts"))) == 0)
      interest_scRNAseq_for_GRN =
        base::subset(interest_scRNAseq_for_GRN, features = base::names(zero_frequency[zero_frequency <= alpha_gene]))

      zero_frequency =
        base::colMeans(base::as.data.frame(base::as.matrix(SeuratObject::GetAssayData(interest_scRNAseq_for_GRN, layer = "counts"))) == 0)
      interest_scRNAseq_for_GRN =
        base::subset(interest_scRNAseq_for_GRN, cells = base::names(zero_frequency[zero_frequency <= alpha_cell]))
    }

    # Extracting activity matrix for cell types of interest (retaining genes/cells with non-zero activity values).
    message("Extracting activity matrix for Cell types of interest (retaining genes/cells with non-zero activity values).")
    interest_scATACseq_for_GRN = base::subset(scATACseq_for_GRN, subset = seurat_annotations == cell_type)
    interest_scATACseq_for_GRN =
      interest_scATACseq_for_GRN[base::rownames(interest_scRNAseq_for_GRN), base::colnames(interest_scRNAseq_for_GRN)]

    # Update data.
    message("Update data.")
    interest_gene_for_GRN = base::intersect(gene_for_GRN, base::rownames(interest_scRNAseq_for_GRN))
    interest_TGs = base::intersect(TGs, interest_gene_for_GRN)
    interest_TFs = base::intersect(TFs, interest_gene_for_GRN)
    interest_gene_for_GRN = base::union(interest_TGs, interest_TFs)
    interest_scRNAseq_for_GRN =
      interest_scRNAseq_for_GRN[base::intersect(base::rownames(interest_scRNAseq_for_GRN), interest_gene_for_GRN), ]
    interest_scATACseq_for_GRN = interest_scATACseq_for_GRN[base::rownames(interest_scRNAseq_for_GRN), ]

    base::stopifnot(base::identical(
      base::rownames(interest_scRNAseq_for_GRN),
      base::rownames(interest_scATACseq_for_GRN)
    ))

    # End.
    return(base::list("interest_scRNAseq_for_GRN" = interest_scRNAseq_for_GRN,
                      "interest_scATACseq_for_GRN" = interest_scATACseq_for_GRN,
                      "interest_gene_for_GRN" = interest_gene_for_GRN,
                      "interest_TGs" = interest_TGs,
                      "interest_TFs" = interest_TFs))
  }

  ### Iterating inner function.
  interest_cell_type_data = parallel::mclapply(
    X = interest_cell_type,
    FUN = function(cell_type) {get_interest_cell_type_data_process(
      alpha_gene = alpha_gene,
      alpha_cell = alpha_cell,
      cell_type = cell_type,
      scRNAseq_for_GRN = Basic_Info[["scRNAseq_for_GRN"]],
      scATACseq_for_GRN = Basic_Info[["scATACseq_for_GRN"]],
      gene_for_GRN = Basic_Info[["gene_for_GRN"]],
      TGs = Basic_Info[["TGs"]],
      TFs = Basic_Info[["TFs"]]
    )},
    mc.cores = ncores)
  base::names(interest_cell_type_data) = interest_cell_type

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd())
  t_end = base::Sys.time()
  message("Finish: Extracting Data for Cell Types of Interest ", t_end)
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_data)
}
