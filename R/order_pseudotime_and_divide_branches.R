#' @title Perform Pseudotime Ordering and Branch Partitioning
#'
#' @description
#' This function performs pseudotime analysis using Monocle2 and enables interactive branch partitioning based on trajectory visualization results.
#' It processes scRNA-seq data for multiple cell types to reconstruct developmental trajectories and identify branching points.
#'
#' @param interest_cell_type_data The output of function get_interest_cell_type_data.
#'
#' @returns
#' A nested list structure where:
#' \itemize{
#' \item First level: Named by cell types
#' \item Second level: For each cell type, a list of data frames containing:
#' \itemize{
#' \item \code{Pseudotime}: Pseudotime values for cells in the branch
#' \item \code{State}: Cell state assignments for the branch
#' }
#' \item Cells are ordered by pseudotime within each branch
#' }
#'
#' @details
#' The function performs the following steps for each cell type:
#' \enumerate{
#' \item Converts Seurat object to Monocle CellDataSet
#' \item Estimates size factors and dispersions
#' \item Performs dimensionality reduction using DDRTree
#' \item Orders cells by pseudotime
#' \item Generates visualization plots (cell states and pseudotime)
#' \item Interactively queries user for branch definitions based on visualizations
#' \item Extracts and orders cells for each defined branch
#' }
#'
#' @section Interactive Input:
#' The function requires user interaction to:
#' \itemize{
#' \item Specify the number of branches based on visualization results
#' \item Define cell state sequences for each branch
#' }
#' Users should examine the generated plots in the output directory before providing input.
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Requires \code{igraph} package version 2.0.3 for compatibility with Monocle
#' \item Generates multiple visualization files in the output directory
#' \item Creates subdirectories for each cell type's analysis results
#' \item The function will pause for user input during execution
#' \item Save the returned object as it contains valuable trajectory information
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Perform pseudotime analysis and branch partitioning
#' library(DDRTree)
#' packageVersion("igraph")
#' remove.packages("igraph")
#' install_version("igraph", version = "2.0.3")
#' interest_cell_type_data = base::readRDS("./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_data.rds")
#' interest_cell_type_Branches = order_pseudotime_and_divide_branches(interest_cell_type_data)
#' # Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): 6
#' # Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).
#' # Please enter the integer vector of cell state order for branch 1 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#' #   1 11
#' # Please enter the integer vector of cell state order for branch 2 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#' #   1 2 10
#' # Please enter the integer vector of cell state order for branch 3 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#' #   1 2 3 9
#' # Please enter the integer vector of cell state order for branch 4 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#' #   1 2 3 4 8
#' # Please enter the integer vector of cell state order for branch 5 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#' #   1 2 3 4 5 6
#' # Please enter the integer vector of cell state order for branch 6 (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5):
#' #   1 2 3 4 5 7
#' base::saveRDS(interest_cell_type_Branches, file = "./2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/interest_cell_type_Branches.rds")
#' }
order_pseudotime_and_divide_branches = function(interest_cell_type_data = interest_cell_type_data)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Performing pseudotime ordering and branch partitioning for cell types of interest ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining a function for pseudotime ordering and branch partitioning of individual cell types.
  order_pseudotime_and_divide_branches_process = function (interest_scRNAseq_for_GRN, interest_gene_for_GRN) {
    # Converting counts expression matrix to sparse matrix format and creating a CellDataSet object.
    message("Converting counts expression matrix to sparse matrix format and creating a CellDataSet object.")
    interest_RNAsparseMatrix_for_GRN =
      methods::as(base::as.matrix(SeuratObject::JoinLayers(interest_scRNAseq_for_GRN)@assays$RNA$counts), 'sparseMatrix')
    interest_cds = monocle::newCellDataSet(
      cellData = interest_RNAsparseMatrix_for_GRN,
      phenoData =
        methods::new('AnnotatedDataFrame',
                     data = SeuratObject::JoinLayers(interest_scRNAseq_for_GRN)@meta.data),
      featureData =
        methods::new('AnnotatedDataFrame',
                     data = base::data.frame(gene_short_name = base::row.names(interest_RNAsparseMatrix_for_GRN),
                                             row.names = base::row.names(interest_RNAsparseMatrix_for_GRN))),
      lowerDetectionLimit = 0.1,
      expressionFamily = VGAM::negbinomial.size()
    )

    # Estimating normalization factors and dispersion.
    message("Estimating normalization factors and dispersion.")
    interest_cds = BiocGenerics::estimateSizeFactors(interest_cds)
    interest_cds = BiocGenerics::estimateDispersions(interest_cds)

    # Specifying henes for subsequent cell clustering or trajectory inference analysis.
    message("Specifying genes for subsequent cell clustering or trajectory inference analysis.")
    interest_cds = monocle::setOrderingFilter(interest_cds, interest_gene_for_GRN)

    # Dimensionality reduction.
    message("Dimensionality reduction.")
    # library(DDRTree)
    interest_cds = monocle::reduceDimension(interest_cds, method = 'DDRTree', max_components = 2)

    # Ordering cells by trajectory (pseudotime).
    message("Ordering cells by trajectory (pseudotime).")
    message("It is recommended to use igraph R package version 2.0.3, as later versions are not compatible with monocle.")
    interest_cds = monocle::orderCells(interest_cds)
    message("Saving the CDS variable 'interest_cds' for monocle pseudotime analysis.")
    base::saveRDS(interest_cds, file = "interest_cds.rds")

    # Generating 2D visualizations by pseudotime and cell state saved in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment.
    message("Generating 2D visualizations by pseudotime and cell states saved in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment.")
    ggplot2::ggsave("interest_AllCellState_in_trajectory.png", width = 5, height = 5,
                    plot = monocle::plot_cell_trajectory(interest_cds, color_by = "State"))
    ggplot2::ggsave(
      "interest_EveryCellState_in_trajectory.png",
      width = 3 * base::ceiling(
        base::length(base::unique(Biobase::pData(interest_cds)$State)) /
          base::ceiling(base::sqrt(base::length(base::unique(Biobase::pData(interest_cds)$State))))
      ),
      height = 3 * base::ceiling(
        base::sqrt(base::length(base::unique(Biobase::pData(interest_cds)$State)))
      ),
      plot = monocle::plot_cell_trajectory(interest_cds, color_by = "State") +
        ggplot2::facet_wrap(
          ~State,
          nrow = base::ceiling(base::sqrt(base::length(base::unique(Biobase::pData(interest_cds)$State)))),
          axes = "all"
        )
    )
    ggplot2::ggsave("interest_Pseudotime_in_trajectory.png", width = 5, height = 5,
                    plot = monocle::plot_cell_trajectory(interest_cds, color_by = "Pseudotime"))

    # Interactively input the number of branches.
    message("Interactively input the number of branches (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).")
    n_branch = base::as.integer(base::readline(prompt = "Please enter the number of branches based on the pseudotime 2D visualization (e.g., 3) (e.g., 5): "))
    if (base::is.na(n_branch) || n_branch <= 1.9999999) {
      stop("Invalid input. Please enter a positive integer for the number of branches.")
    }

    # Interactively input integer vector of cell state order for each branch.
    message("Interactively input the integer vector of cell state order for each branch (based on visualization results in: 2.1 Data Processing - Pseudotime Analysis And Cell Branching Assignment/cell_type/).")
    branches = base::list()
    for (i in 1:n_branch) {
      base::cat(base::sprintf("Please enter the integer vector of cell state order for branch %d (space-separated) (e.g. 1 2 3 4) (e.g. 1 2 5): ", i))
      integer_vector = base::as.integer(base::strsplit(base::readline(), " ")[[1]])
      if (base::any(base::is.na(integer_vector))) {
        stop(base::sprintf("Invalid input. Please ensure the input for branch %d is an integer vector.", i))
      }
      branches[[i]] = integer_vector
    }

    # Extracting pseudotime and cell states for each branch.
    message("Extracting pseudotime and cell states for each branch.")
    Branches = base::list()
    for (i in base::seq_along(branches)) {
      Branches[[i]] = Biobase::pData(interest_cds)[
        Biobase::pData(interest_cds)$State %in% branches[[i]],
      ][, c("Pseudotime", "State")]
      Branches[[i]] = Branches[[i]][base::order(Branches[[i]][, 1]), ]
    }

    return(Branches)
  }

  ### Executing the defined function for each cell type.
  interest_cell_type_Branches = stats::setNames(
    base::vector("list", base::length(interest_cell_type_data)), base::names(interest_cell_type_data)
  )
  for (name in base::names(interest_cell_type_Branches)) {
    # Start.
    message("Performing Pseudotime Ordering and Branch Partitioning for ", name, ".")
    original_dir1 = base::getwd()
    new_folder1 = name
    if (!base::dir.exists(new_folder1)) {
      base::dir.create(new_folder1, recursive = TRUE)
      message("Folder already creates: ", new_folder1, ".")
    } else {message("Folder already exists: ", new_folder, ".")}
    base::setwd(new_folder1)
    message("The current working directory has been switched to: ", base::getwd(), ".")

    # Do.
    interest_cell_type_Branches[[name]] = order_pseudotime_and_divide_branches_process(
      interest_scRNAseq_for_GRN = interest_cell_type_data[[name]]$interest_scRNAseq_for_GRN,
      interest_gene_for_GRN = interest_cell_type_data[[name]]$interest_gene_for_GRN
    )

    # End.
    base::setwd(original_dir1)
    message("The current working directory has been switched to: ", base::getwd(), ".")
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Performing pseudotime ordering and branch partitioning for cell types of interest ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_Branches)
}
