#' @title Binarization Threshold of GRN for Each Cell Type
#'
#' @description
#' Computes a cell-type-specific threshold (tao) for binarizing GRN.
#' The threshold is calculated as the minimum non-zero interaction strength in the iGRN minus a small constant, with a lower bound of 0.005.
#'
#' @param interest_cell_type_iGRN The output of function get_iGRN_by_TFBS_pwm_by_JASPAR2024.
#'
#' @returns A named numeric vector where names correspond to cell types and values represent the computed binarization thresholds (tao) for each cell type.
#'
#' @details
#' The function performs the following steps:
#' 1. For each cell type, identifies all non-zero interaction strengths (theta) in the iGRN.
#' 2. Computes the minimum non-zero theta value and subtracts 0.005.
#' 3. Sets the threshold (tao) to the maximum of this computed value and 0.005, ensuring a lower bound of 0.005.
#' 4. Returns a vector of thresholds, one per cell type.
#'
#' This threshold is typically used to convert continuous interaction strengths into binary values (presence/absence of regulatory interactions) for downstream analyses.
#'
#' @note
#' Important considerations:
#' - The threshold is cell-type-specific, allowing for adaptive binarization based on the distribution of interaction strengths in each cell type.
#' - The lower bound of 0.005 prevents overly sensitive thresholds in cases where the minimum non-zero theta is very small.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_iGRN.rds")
#' interest_cell_type_tao = get_cell_type_tao(interest_cell_type_iGRN)
#' base::saveRDS(interest_cell_type_tao, file = "./3 get iGRN/interest_cell_type_tao.rds")
#' }
get_cell_type_tao = function(interest_cell_type_iGRN)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting tao for each cell type ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Process.
  message("Initializing.")
  interest_cell_type_tao = base::rep(999, base::length(interest_cell_type_iGRN))
  base::names(interest_cell_type_tao) = base::names(interest_cell_type_iGRN)
  message("Calulating.")
  for (cell_type in base::names(interest_cell_type_iGRN)) {
    interest_cell_type_tao[cell_type] = base::max(
      base::min(interest_cell_type_iGRN[[cell_type]][interest_cell_type_iGRN[[cell_type]] > 0], na.rm = TRUE) - 0.005,
      0.005
    )
  }

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Extracting tao for each cell type ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_tao)
}
