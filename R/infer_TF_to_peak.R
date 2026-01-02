#' @title Infer TF-Peak Regulatory Relationships
#'
#' @description
#' Infers regulatory relationships between transcription factors (TFs) and chromatin accessibility peaks by integrating filtered GRN predictions with peak-TG associations.
#' The function processes each cell type at two stringency levels (loose and strict) to identify TFs present in predicted GRNs, maps them to their target genes, and identifies peaks associated with these target genes through promoter annotations.
#'
#' @param interest_cell_type_pGRN The output of function filter_GRN.
#' @param interest_cell_type_peak_TG_pred The output of function infer_peak_to_TG.
#' @param interest_cell_type_iGRN The output of function get_iGRN_by_TFBS_pwm_by_JASPAR2024.
#' @param ncores See ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list named \code{interest_cell_type_TF_peak_pred} organized by cell type and filtering level, containing for each level:
#' \itemize{
#' \item \code{TFs}: Character vector of unique transcription factor symbols present in the filtered GRN predictions
#' \item \code{TF_TGs}: Data frame mapping TFs to their target genes, with columns "TF" and "TG"
#' \item \code{TF_peaks}: Nested list for each TF, containing:
#' \itemize{
#' \item \code{TGs}: Target genes regulated by the TF
#' \item \code{peaks}: Named character vector of promoter peaks associated with the target genes, where names are the target gene symbols and values represent peak genomic coordinates
#' \item Regulatory scores including theta_i (TF-target interaction strength), theta_p (peak accessibility influence), fit_loss (model fitting error) with their percentile rankings, and so on (also see in ?infer_peak_to_TG)
#' }
#' }
#'
#' @details
#' The function performs the following operations for each cell type and filtering level:
#' \enumerate{
#' \item Extracts unique transcription factor symbols from the filtered GRN predictions
#' \item Creates a data frame mapping each TF to its target genes based on the GRN predictions
#' \item For each TF, identifies promoter peaks associated with its target genes using the peak-TG annotations
#' \item Calculates regulatory scores by integrating information from both initial GRN (iGRN) and predicted GRN (pGRN) to quantify TF-peak regulatory relationships
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
#' \item Parallel computation is implemented using mclapply for efficient processing of large datasets
#' }
#' Afterwards, users can filter predicted TGs based on TG_importance_store or TG_importance_score_percentile, filter predicted peaks based on peak_importance_store or peak_importance_score_percentile, and filter predicted peak TG associations based on peak_TG_importance_store or peak_TG_importance_score_percentile.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_peak_TG_pred = base::readRDS("./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_TF_peak_pred = infer_TF_to_peak(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_peak_TG_pred = interest_cell_type_peak_TG_pred,
#'   interest_cell_type_iGRN = interest_cell_type_iGRN,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_TF_peak_pred, file = "./5 Infer GRN/interest_cell_type_TF_peak_pred.rds")
#' }
infer_TF_to_peak = function(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    interest_cell_type_peak_TG_pred = interest_cell_type_peak_TG_pred,
    interest_cell_type_iGRN = interest_cell_type_iGRN,
    ncores = ncores
){
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

      tf_list = interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TFs"]]

      interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]] = base::list()
      process_tf = function(tf) {
        result = base::list()
        result[["tf"]] = tf
        result[["TGs"]] = base::as.character(interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]][
          interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_TGs"]]$TF == tf,
          "TG"
        ])

        TGs = result[["TGs"]]
        peaks = interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]][
          interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]]$SYMBOL %in% TGs, , drop = FALSE
        ]

        peaks$theta_i = 0
        peaks$theta_i_percentile = NA
        peaks$theta_p = 0
        peaks$theta_p_percentile = NA
        peaks$fit_loss = Inf
        peaks$fit_loss_percentile = NA

        iGRN = interest_cell_type_iGRN[[cell_type]][TGs, tf, drop = FALSE]
        pGRN = interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]][, c("theta_p", "fit_loss")]

        for (tg in TGs) {
          if (tg %in% base::rownames(iGRN)) {
            peaks[peaks$SYMBOL == tg, "theta_i"] = iGRN[tg, tf]
          }
          re = base::paste0(tf, "_to_", tg)
          if (re %in% base::rownames(pGRN)) {
            peaks[peaks$SYMBOL == tg, "theta_p"] = pGRN[re, "theta_p"]
            peaks[peaks$SYMBOL == tg, "fit_loss"] = pGRN[re, "fit_loss"]
          }
        }

        rank_theta_i = base::rank(peaks$theta_i)
        peaks$theta_i_percentile = rank_theta_i / base::max(rank_theta_i)

        rank_theta_p = base::rank(peaks$theta_p)
        peaks$theta_p_percentile = rank_theta_p / base::max(rank_theta_p)

        rank_fit_loss = base::rank(peaks$fit_loss)
        peaks$fit_loss_percentile = rank_fit_loss / base::max(rank_fit_loss)

        result[["peaks"]] = peaks
        return(result)
      }

      results = parallel::mclapply(tf_list, process_tf, mc.cores = ncores)

      for (result in results) {
        tf = result[["tf"]]
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]][[tf]] = base::list()
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]][[tf]][["TGs"]] = result[["TGs"]]
        interest_cell_type_TF_peak_pred[[cell_type]][[level]][["TF_peaks"]][[tf]][["peaks"]] = result[["peaks"]]
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
