#' @title Infer Peak-TG Relationships
#'
#' @description
#' Infers regulatory relationships between chromatin accessibility peaks and target genes (TGs) by matching filtered GRN predictions with promoter-annotated peaks.
#' The function processes each cell type independently at two stringency levels (loose and strict) to identify target genes present in the predicted GRN, then retrieves promoter peaks associated with these genes.
#'
#' @param interest_cell_type_pGRN The output of function filter_GRN.
#' @param results_identify_TGs The output of function identify_TGs.
#'
#' @returns
#' A nested list named \code{interest_cell_type_peak_TG_pred} organized by cell type and filtering level, containing for each level:
#' \itemize{
#' \item \code{TGs}: Character vector of unique target gene symbols present in the filtered GRN predictions
#' \item \code{peak_TG}: Subset of promoter-annotated peaks where the associated gene symbol matches one of the identified TGs
#' }
#'
#' @details
#' The function performs the following operations for each cell type:
#' \enumerate{
#' \item Extracts target gene symbols from the loose and strict filtered GRN predictions
#' \item Identifies unique target genes at each filtering level
#' \item Subsets the promoter peak annotation data to include only peaks associated with the identified target genes
#' \item Organizes results in a structured list for downstream analysis
#' }
#' This enables integration of chromatin accessibility data with GRN predictions to identify potential regulatory elements.
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Requires promoter-annotated peaks in \code{results_identify_TGs}; ensure proper annotation format
#' \item The "loose" level provides more comprehensive target gene sets, while "strict" yields high-confidence targets
#' \item Missing promoter annotations for some target genes may result in incomplete peak-gene mappings
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
#' interest_cell_type_peak_TG_pred = infer_peak_to_TG (
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   results_identify_TGs = results_identify_TGs
#' )
#' base::saveRDS(interest_cell_type_peak_TG_pred, file = "./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")
#' }
infer_peak_to_TG = function(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    results_identify_TGs = results_identify_TGs
)
{
  t_start = base::Sys.time()
  message("Run: Inferring Peak-TG Relationships ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  interest_cell_type_peak_TG_pred = base::list()
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    message("Inferring peak-TG relationships of cell type ", cell_type, ".")
    interest_cell_type_peak_TG_pred[[cell_type]] = base::list()

    for (level in c("loose", "strict")) {
      message("Inferring peak-TG relationships of cell type ", cell_type, " in ", level, " level.")
      interest_cell_type_peak_TG_pred[[cell_type]][[level]] = base::list()
      interest_cell_type_peak_TG_pred[[cell_type]][[level]][["TGs"]] =
        base::unique(base::sapply(base::strsplit(base::rownames(
          interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]]
        ), "_to_"), function(x) x[2]))
      interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]] = base::subset(
        results_identify_TGs[["peak_anno_promoter"]][, "SYMBOL", drop = FALSE],
        SYMBOL %in% interest_cell_type_peak_TG_pred[[cell_type]][[level]][["TGs"]]
      )
    }
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring Peak-TG Relationships ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_peak_TG_pred)
}
