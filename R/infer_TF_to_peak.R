#' @title Infer TF-Peak Regulatory Relationships
#'
#' @description
#' Infers regulatory relationships between transcription factors (TFs) and chromatin accessibility peaks by integrating filtered GRN predictions with peak-TG associations.
#' The function processes each cell type at two stringency levels (loose and strict) to identify TFs present in predicted GRNs, maps them to their target genes, and identifies peaks associated with these target genes through promoter annotations.#'
#'
#' @param interest_cell_type_pGRN The output of function filter_GRN.
#' @param interest_cell_type_peak_TG_pred The output of function infer_peak_to_TG.
#'
#' @returns
#' A nested list named \code{interest_cell_type_TF_peak_pred} organized by cell type and filtering level, containing for each level:
#' \itemize{
#' \item \code{TFs}: Character vector of unique transcription factor symbols present in the filtered GRN predictions
#' \item \code{TF_TGs}: Data frame mapping TFs to their target genes, with columns "TF" and "TG"
#' \item \code{TF_peaks}: Nested list for each TF, containing:
#' \itemize{
#' \item \code{TGs}: Target genes regulated by the TF
#' \item \code{peaks}: Named character vector of promoter peaks associated with the target genes, where names are the target gene symbols
#' }
#' }
#'
#' @details
#' The function performs the following operations for each cell type and filtering level:
#' \enumerate{
#' \item Extracts unique transcription factor symbols from the filtered GRN predictions
#' \item Creates a data frame mapping each TF to its target genes based on the GRN predictions
#' \item For each TF, identifies promoter peaks associated with its target genes using the peak-TG annotations
#' \item Organizes results in a structured nested list, enabling analysis of TF regulatory networks with chromatin accessibility data
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Requires both GRN predictions and peak-TG annotations as input
#' \item The "loose" level provides more comprehensive TF-peak associations, while "strict" yields high-confidence interactions
#' \item Peaks are identified through promoter annotations of target genes, which may not capture distal regulatory elements
#' \item The output structure supports downstream analysis of TF regulatory networks
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_peak_TG_pred = base::readRDS("./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")
#' interest_cell_type_TF_peak_pred = infer_TF_to_peak(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_peak_TG_pred = interest_cell_type_peak_TG_pred
#' )
#' base::saveRDS(interest_cell_type_TF_peak_pred, file = "./5 Infer GRN/interest_cell_type_TF_peak_pred.rds")
#' }
infer_TF_to_peak = function(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    interest_cell_type_peak_TG_pred = interest_cell_type_peak_TG_pred
)
{
  t_start = base::Sys.time()
  message("Run: Inferring TF-Peak Relationships ", t_start, ".")
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

  interest_cell_type_TF_peak_pred = base::list()
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    message("Inferring TF-Peak relationships of cell type ", cell_type, ".")
    interest_cell_type_TF_peak_pred[[cell_type]] = base::list()

    for (level in c("loose", "strict")) {
      message("Inferring TF-Peak relationships of cell type ", cell_type, " in ", level, " level.")
      interest_cell_type_TF_peak_pred[[cell_type]][[level]] = base::list()
      interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TFs"]] =
        base::unique(base::sapply(base::strsplit(base::rownames(
          interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]]
        ), "_to_"), function(x) x[1]))
      interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]] =
        base::as.data.frame(stringr::str_split_fixed(base::rownames(
          interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]]
        ), "_to_", 2))
      base::colnames(interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]]) = c("TF", "TG")
      base::rownames(interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]]) = base::paste0(
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]]$TF,
        "_to_",
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]]$TG
      )
      interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]] = base::list()
      for (tf in interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TFs"]]) {
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]][[tf]] = base::list()
        TGs = data.table::setDT(interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]])[TF == tf, TG]
        peaks = base::rownames(interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]][
          interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]]$SYMBOL %in% TGs, , drop = FALSE
        ])
        for (i in 1:base::length(peaks)) {
          base::names(peaks)[i] = interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]][peaks[i], "SYMBOL"]
        }
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]][[tf]][["TGs"]] = TGs
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]][[tf]][["peaks"]] = peaks
      }
    }
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring TF-Peak Relationships ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_TF_peak_pred)
}
