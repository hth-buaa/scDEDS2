#' @title Prepare Data for Gene Regulatory Network Inference
#'
#' @description
#' Prepares comprehensive datasets for gene regulatory network (GRN) inference by combining iGRN with cell group state features.
#' The function generates all possible TF-TG pairs, adds their corresponding regulatory strength values, and integrates cell group state measurements (TFE, TGA, TGE) for each time point.
#' Data is processed in parallel for efficiency and can be partitioned to handle large datasets.
#'
#' @param interest_cell_type_branch_iGRN The output of function get_branch_iGRN.
#' @param interest_cell_type_group The output of function cell_grouping.
#' @param n_part Integer. Number of partitions to slice the data for parallel processing. Default is 1 (no partitioning).
#' @param ncores See in ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list organized by cell type and branch, containing:
#' \itemize{
#' \item For each branch: A data frame with all potential TF-TG pairs and their features, including:
#' \itemize{
#' \item \code{TG}: Target gene identifier
#' \item \code{TF}: Transcription factor identifier
#' \item \code{theta_i}: Regulatory strength from iGRN
#' \item \code{TFE_K}, \code{TGA_K}, \code{TGE_K}: Cell group state features for each time point K
#' }
#' \item A "GRN" element containing the combined dataset across all branches
#' }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#' \item For each cell type, processes each branch independently
#' \item Generates all possible TF-TG combinations from the iGRN matrix
#' \item Extracts regulatory strength (theta_i) for each TF-TG pair
#' \item Partitions the data into \code{n_part} subsets for parallel processing
#' \item In parallel, adds cell group state features (TFE, TGA, TGE) for each time point
#' \item Combines partitioned results and organizes them by branch
#' \item Creates a combined dataset across all branches
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Parallel processing significantly speeds up the feature addition step, especially for large datasets
#' \item The number of partitions (\code{n_part}) can be adjusted based on dataset size and available memory
#' \item The resulting data frame includes all possible TF-TG pairs, which can be large for datasets with many genes
#' \item The "GRN" element provides a consolidated view across all branches for comprehensive analysis
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_iGRN = base::readRDS("./3 get iGRN/interest_cell_type_branch_iGRN.rds")
#' interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' n_part = ncores
#' interest_cell_type_branch_data_for_GRN_infer = get_data_for_inferring_GRN(
#'   interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
#'   interest_cell_type_group = interest_cell_type_group,
#'   n_part = n_part,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_branch_data_for_GRN_infer, file = "./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")
#' }
get_data_for_inferring_GRN = function(
    interest_cell_type_branch_iGRN = interest_cell_type_branch_iGRN,
    interest_cell_type_group = interest_cell_type_group,
    n_part = 1,
    ncores = 1
)
{
  t_start = base::Sys.time()
  message("Run: Getting data for inferring GRN for each branch ", t_start, ".")
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

  interest_cell_type_branch_data_for_GRN_infer = base::list()
  for (cell_type in names(interest_cell_type_branch_iGRN)) {
    message("Process cell type ", cell_type, ".")
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]] = base::list()
    branch_GRN = interest_cell_type_branch_iGRN[[cell_type]]

    message("Generating all row-column (TG-TF) combinations, and adding regulatory strength theta for each TF_TG pair.")
    n_branch = base::length(branch_GRN)
    all_combinations = base::list()
    for (n in 1:n_branch) {
      all_combinations[[base::paste0("branch_", n, "_GRN")]] = base::expand.grid(
        TG = base::rownames(branch_GRN[[n]]), TF = base::colnames(branch_GRN[[n]]), stringsAsFactors = FALSE
      )
      all_combinations[[base::paste0("branch_", n, "_GRN")]]$theta_i = base::apply(
        all_combinations[[base::paste0("branch_", n, "_GRN")]], 1, function(x) branch_GRN[[n]][x[1], x[2]]
      )
      rownames(all_combinations[[base::paste0("branch_", n, "_GRN")]]) = paste0(
        all_combinations[[base::paste0("branch_", n, "_GRN")]]$TF,
        "_to_",
        all_combinations[[base::paste0("branch_", n, "_GRN")]]$TG
      )
    }
    all_combinations[["GRN"]] = base::do.call(rbind, all_combinations)
    all_combinations[["GRN"]] = all_combinations[["GRN"]][!base::duplicated(all_combinations[["GRN"]]), ]
    rownames(all_combinations[["GRN"]]) = paste0(
      all_combinations[["GRN"]]$TF,
      "_to_",
      all_combinations[["GRN"]]$TG
    )

    message("Adding cell group states (TFE, TGA, TGE).")
    for (n in 1:n_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      GRN = all_combinations[[base::paste0("branch_", n, "_GRN")]]
      N_re_all = base::nrow(GRN)
      message("There are ", N_re_all, " potential regulatory relationship.")

      message("Partitioning data.")
      each = base::floor(N_re_all / n_part)
      each_last = N_re_all - each * (n_part - 1)
      GRN_list = base::list()
      for (s in 1:n_part) {
        if (s != n_part) {
          GRN_list[[s]] = GRN[(each * (s - 1) + 1):(each * s), ]
        } else (
          GRN_list[[s]] = GRN[(each * (s - 1) + 1):N_re_all, ]
        )
      }

      message("Adding cell group states (TFE, TGA, TGE) multithreadedly.")
      process_GRN_slice = function(s,
                                   GRN_list = GRN_list,
                                   TFs = base::colnames(branch_GRN[[n]]),
                                   TGs = base::rownames(branch_GRN[[n]]),
                                   TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
                                   TGA_T = interest_cell_type_group[[cell_type]][["Branches_TGA_T"]],
                                   TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]) {
        current_branch = GRN_list[[s]]
        N_re = base::nrow(current_branch)
        res = base::rownames(current_branch)
        for (TF_TG in res) {
          if (current_branch[TF_TG, "TF"] %in% TFs && current_branch[TF_TG, "TG"] %in% TGs) {
            for (K in 1:ncol(TFE_T[[n]])) {
              current_branch[TF_TG, c(
                base::paste0("TFE_", K), base::paste0("TGA_", K), base::paste0("TGE_", K)
              )] = c(TFE_T[[n]][current_branch[TF_TG, "TF"], K],
                     TGA_T[[n]][current_branch[TF_TG, "TG"], K],
                     TGE_T[[n]][current_branch[TF_TG, "TG"], K])
            }
          } else {
            current_branch = current_branch[!rownames(current_branch) %in% res[n_re], ]
          }
        }
        return(current_branch)
      }
      data_for_prediction_list = parallel::mclapply(
        X = 1:n_part,
        FUN = function(s) {
          process_GRN_slice(
            s,
            GRN_list = GRN_list,
            TFs = base::colnames(branch_GRN[[n]]),
            TGs = base::rownames(branch_GRN[[n]]),
            TFE_T = interest_cell_type_group[[cell_type]][["Branches_TFE_T"]],
            TGA_T = interest_cell_type_group[[cell_type]][["Branches_TGA_T"]],
            TGE_T = interest_cell_type_group[[cell_type]][["Branches_TGE_T"]]
          )
        },
        mc.cores = ncores
      )
      interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]] = base::do.call("rbind", data_for_prediction_list)
    }
    interest_cell_type_branch_data_for_GRN_infer[[cell_type]][["GRN"]] = all_combinations[["GRN"]]
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Getting data for inferring GRN for each branch ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_branch_data_for_GRN_infer)
}
