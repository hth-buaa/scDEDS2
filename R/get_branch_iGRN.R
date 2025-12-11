#' @title Extract Branch-Specific Initial Gene Regulatory Networks (iGRNs)
#'
#' @description
#' This function extracts branch-specific initial gene regulatory networks (iGRNs) from the complete cell-type-specific iGRNs.
#' It creates subnetworks for each branch within each cell type by filtering the global TF-TG association matrix to include only the TFs and TGs present in each specific branch.
#'
#' @param interest_cell_type_iGRN_all_TGTF_pairs The output of function get_iGRN_by_TFBS_pwm_by_JASPAR2024.
#' @param interest_cell_type_group The output of function cell_grouping.
#' @param ncores See in ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list structure where:
#' \itemize{
#' \item First level: Named by cell types
#' \item Second level: For each cell type, a list of branch-specific iGRNs
#' \item Each branch iGRN is a data frame with:
#' \itemize{
#' \item Rows: Target genes (TGs) present in the branch
#' \item Columns: Transcription factors (TFs) present in the branch
#' \item Values: Normalized PWM match scores from the complete iGRN
#' }
#' }
#' The branch iGRNs preserve the regulatory strength values from the complete network but are filtered to include only the TFs and TGs that are expressed/active in each branch.
#'
#' @details
#' The function performs the following operations for each cell type and branch:
#' \enumerate{
#' \item Identifies TFs and TGs present in each branch using grouping information
#' \item Creates branch-specific subnetwork matrices initialized with zeros
#' \item Populates the matrices with regulatory strength values from the complete iGRN
#' \item Returns organized branch-specific networks for downstream analysis
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item This function does not recalculate regulatory strengths but filters existing networks
#' \item Branch networks may be sparse if few TFs/TGs are shared between branches
#' \item The function uses parallel processing for efficiency with multiple cell types
#' \item Output directory structure is maintained for consistency with other functions
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_iGRN_all_TGTF_pairs = base::readRDS("./3 get iGRN/interest_cell_type_iGRN_all_TGTF_pairs.rds")
#' interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_branch_iGRN = get_branch_iGRN(
#'   interest_cell_type_iGRN_all_TGTF_pairs = interest_cell_type_iGRN_all_TGTF_pairs,
#'   interest_cell_type_group = interest_cell_type_group,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_branch_iGRN, file = "./3 get iGRN/interest_cell_type_branch_iGRN.rds")
#' }
get_branch_iGRN = function(
    interest_cell_type_iGRN_all_TGTF_pairs = interest_cell_type_iGRN_all_TGTF_pairs,
    interest_cell_type_group = interest_cell_type_group,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Extracting iGRN for each branch from every cell type's iGRN ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defining inner function.
  get_branch_iGRN_process = function(cell_type, interest_cell_type_iGRN_all_TGTF_pairs, interest_cell_type_group)
  {
    branch_iGRN = base::list()
    for (n in 1:base::length(interest_cell_type_group[[cell_type]][["n_min"]])) {
      branch_TGs = base::rownames(interest_cell_type_group[[cell_type]][["Branches_TGE_T"]][[n]])
      branch_TFs = base::rownames(interest_cell_type_group[[cell_type]][["Branches_TFE_T"]][[n]])
      branch_iGRN[[n]] = base::as.data.frame(base::matrix(0, nrow = base::length(branch_TGs), ncol = base::length(branch_TFs)))
      base::rownames(branch_iGRN[[n]]) = branch_TGs
      base::colnames(branch_iGRN[[n]]) = branch_TFs
      for (TG in branch_TGs) {
        for (TF in branch_TFs) {
          if (TG %in% base::rownames(interest_cell_type_iGRN_all_TGTF_pairs[[cell_type]]) &&
              TF %in% base::colnames(interest_cell_type_iGRN_all_TGTF_pairs[[cell_type]])) {
            branch_iGRN[[n]][TG, TF] = interest_cell_type_iGRN_all_TGTF_pairs[[cell_type]][TG, TF]
          }
        }
      }
    }
    return(branch_iGRN)
  }

  ### Invoking functions in parallel.
  interest_cell_type_branch_iGRN = parallel::mclapply(
    X = base::names(interest_cell_type_iGRN_all_TGTF_pairs),
    FUN = function(cell_type) {
      message("Preparing to extract iGRN for each branch in cell type ", cell_type, ".")
      get_branch_iGRN_process(
        cell_type = cell_type,
        interest_cell_type_iGRN_all_TGTF_pairs = interest_cell_type_iGRN_all_TGTF_pairs[cell_type],
        interest_cell_type_group = interest_cell_type_group[cell_type]
      )
    },
    mc.cores = ncores
  )
  base::names(interest_cell_type_branch_iGRN) = base::names(interest_cell_type_iGRN_all_TGTF_pairs)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Extracting iGRN for each branch from every cell type's iGRN ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_branch_iGRN)
}
