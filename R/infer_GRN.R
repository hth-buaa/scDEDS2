#' @title Infer Gene Regulatory Networks Using Trained Models
#'
#' @description
#' Infers gene regulatory networks (GRN) by applying trained models to predict regulatory strengths for all potential TF-TG interactions.
#' The function processes data in parallel, calculates predicted regulatory strengths (theta_p) and fitting losses, and organizes results by cell type and branch.
#' It supports selective processing of specific cell types and branches for targeted analysis.
#'
#' @param interest_cell_type_branch_data_for_GRN_infer The output of function get_data_for_inferring_GRN.
#' @param interest_cell_type_tao The output of function get_cell_type_tao.
#' @param interest_cell_type_branch_model_train The output of function model_train.
#' @param interest_cell_type_group The output of function cell_grouping.
#' @param interest_cell_type_genes_pseudotime_info The output of function get_genes_pseudotime_info.
#' @param select_cell_type Character vector or NULL. Specific cell types to process. If NULL, processes all cell types. Default is NULL.
#' @param select_branch Numeric vector or NULL. Specific branches to process. If NULL, processes all branches. Default is NULL.
#' @param each Integer. Maximum number of regulatory interactions to process in each parallel chunk. Default is 1,000,000.
#' @param ncores See in ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list named \code{interest_cell_type_pGRN} organized by cell type and branch, containing:
#' \itemize{
#' \item For each branch: Two data frames:
#' \itemize{
#' \item \code{theta_p}: Predicted regulatory strengths for all TF-TG pairs (rows) using each trained model (columns)
#' \item \code{loss}: Fitting losses corresponding to each prediction
#' }
#' }
#' Additionally, saves each branch's results as an RDS file in the directory "./5 Infer GRN/".
#'
#' @details
#' The inference process follows these steps:
#' \enumerate{
#' \item \strong{Data Preparation}: Selects target cell types/branches
#' \item \strong{Parallel Processing}: Partitions data into chunks of size \code{each} for efficient parallel computation
#' \item \strong{Model Application}: For each TF-TG pair, applies all trained models to predict regulatory strength:
#' \itemize{
#' \item Uses optimized parameters from training
#' \item Predicts theta_p by minimizing the prediction error
#' }
#' \item \strong{Result Aggregation}: Combines results from parallel chunks and saves as RDS files
#' \item \strong{Output Organization}: Returns structured list of predictions and losses
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Processing time scales with the number of TF-TG pairs and trained models
#' \item The \code{each} parameter controls memory usage; reduce if encountering memory issues
#' \item Results are automatically saved as RDS files for persistence
#' \item Use \code{select_cell_type} and \code{select_branch} for targeted analysis
#' \item Parallel processing significantly accelerates computation for large datasets
#' }
#'
#' @seealso \code{\link[scDEDS]{R_cal}}, \code{\link[scDEDS]{Hill_cal}}, \code{\link[scDEDS]{S_cal}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_data_for_GRN_infer = base::readRDS("./5 Infer GRN/interest_cell_type_branch_data_for_GRN_infer.rds")
#' interest_cell_type_tao = base::readRDS("./3 get iGRN/interest_cell_type_tao.rds")
#' interest_cell_type_branch_model_train = base::readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
#' interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
#' each = 1000
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#'
#' # way 1: one step (If the running memory is sufficient)
#' interest_cell_type_pGRN = infer_GRN(
#'   interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'   interest_cell_type_tao = interest_cell_type_tao,
#'   interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'   interest_cell_type_group = interest_cell_type_group,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   select_cell_type = NULL,
#'   select_branch = NULL,
#'   each = each,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#'
#' # way 2: more steps but stable
#' # If the process is terminated due to memory issues, you can manually continue running it branch by branch (saving the results to interest_cell_type_pGRN).
#' interest_cell_type_pGRN = base::list()
#' for (cell_type in base::names(interest_cell_type_branch_data_for_GRN_infer)) {
#'   message("Infer cell type ", cell_type, ".")
#'   interest_cell_type_pGRN[[cell_type]] = base::list()
#'   for (n in 1:(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)) {
#'     message("Infer cell type ", cell_type, " branch ", n, ".")
#'     interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = infer_GRN(
#'       interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
#'       interest_cell_type_tao = interest_cell_type_tao,
#'       interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'       interest_cell_type_group = interest_cell_type_group,
#'       interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'       select_cell_type = cell_type,
#'       select_branch = n,
#'       each = each,
#'       ncores = ncores
#'     )[[cell_type]][[base::paste0("branch", n)]]
#'   }
#' }
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
infer_GRN = function(
    interest_cell_type_branch_data_for_GRN_infer = interest_cell_type_branch_data_for_GRN_infer,
    interest_cell_type_tao = interest_cell_type_tao,
    interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
    interest_cell_type_group = interest_cell_type_group,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    select_cell_type = NULL,
    select_branch = NULL,
    each = 1000000,
    ncores = 1
)
{
  t_start = base::Sys.time()
  message("Run: Inferring GRN ", t_start, ".")
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

  interest_cell_type_pGRN = base::list()
  if (base::is.null(select_cell_type)) {
    range_cell_type = base::names(interest_cell_type_branch_data_for_GRN_infer)
  } else {
    range_cell_type = select_cell_type
  }
  for (cell_type in range_cell_type) {
    interest_cell_type_pGRN[[cell_type]] = base::list()
    tao = as.numeric(interest_cell_type_tao[cell_type])

    if (base::is.null(select_branch)) {
      range_branch = 1:(base::length(interest_cell_type_branch_data_for_GRN_infer[[cell_type]]) - 1)
    } else {
      range_branch = select_branch
    }
    for (n in range_branch) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      re_pred = base::rownames(interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]])
      N_pred_all = base::length(re_pred)
      message(N_pred_all, " potential regulatory relationships will be predicted.")

      re_train = base::sort(base::names(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]))
      message(base::length(re_train), " trained regulatory relationships will be predicted.")

      message("Partitioning data.")
      n_part = base::ceiling(N_pred_all / each)
      each_last = N_pred_all - each * (n_part - 1)

      re_pred_list = base::list()
      for (s in 1:n_part) {
        if (s != n_part) {
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):(each * s)]
        } else (
          re_pred_list[[s]] = re_pred[(each * (s - 1) + 1):N_pred_all]
        )
      }

      pGRN_list = base::list()
      for (s in 1:n_part) {
        pGRN_list[[s]] = base::list()
        pGRN_list[[s]][["theta_p"]] = base::data.frame(base::matrix(
          0, nrow = base::length(re_pred_list[[s]]), ncol = base::length(re_train),
          dimnames = base::list(re_pred_list[[s]], re_train))
        )
        pGRN_list[[s]][["loss"]] = base::data.frame(base::matrix(
          Inf, nrow = base::length(re_pred_list[[s]]), ncol = base::length(re_train),
          dimnames = base::list(re_pred_list[[s]], re_train))
        )
      }

      message("Inferring GRN of branch ", n, " of cell type ", cell_type, ".")
      state_data = interest_cell_type_branch_data_for_GRN_infer[[cell_type]][[base::paste0("branch", n)]]
      TFE_col = base::grep("TFE_", base::colnames(state_data), value = TRUE)
      TGA_col = base::grep("TGA_", base::colnames(state_data), value = TRUE)
      TGE_col = base::grep("TGE_", base::colnames(state_data), value = TRUE)
      Length = base::length(TFE_col)
      seqLength = base::seq_len(base::length(TFE_col) - 1)
      train_params = lapply(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]], \(sublist) sublist[["best_params"]])
      Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]]

      process_data_slice = function(pGRN = base::list(theta_p = NA, loss = NA),
                                    RE_PRED = re_pred_list[[s]],
                                    ncols = base::length(re_train),
                                    STATE_DATA = state_data[re_pred_list[[s]], -(1:3)],
                                    TFE_col = TFE_col,
                                    TGA_col = TGA_col,
                                    TGE_col = TGE_col,
                                    Length = Length,
                                    seqLength = seqLength,
                                    re_train = re_train,
                                    train_params = train_params,
                                    tao = tao,
                                    Tslot_K = Tslot_K) {
        for (re in RE_PRED) {
          theta_p = base::data.frame(base::matrix(0, nrow = 1, ncol = ncols, dimnames = base::list(re, re_train)))
          loss = base::data.frame(base::matrix(Inf, nrow = 1, ncol = ncols, dimnames = base::list(re, re_train)))

          TFE = base::as.numeric(STATE_DATA[re, TFE_col])
          TGA = base::as.numeric(STATE_DATA[re, TGA_col])
          TGE = base::as.numeric(STATE_DATA[re, TGE_col])

          for (re_tr in re_train) {
            pa_best = train_params[[re_tr]]
            predict_result = stats::optim(par = tao, fn = function(theta) {
              R = scDEDS::R_cal(
                theta_TF_TG = theta, tao = tao, r1 = 200,
                r2 = (pa_best["r.r2_TG"] + pa_best["r.r2_TF"])/2,
                r3 = 1, r4 = 1000, r5 = 1
              )
              alpha1 = (pa_best[base::paste0("alpha.alpha1_TG_K-1.", seqLength)] + pa_best[base::paste0("alpha.alpha1_TF_K-1.", seqLength)])/2
              alpha2 = (pa_best[base::paste0("alpha.alpha2_TG_K-1.", seqLength)] + pa_best[base::paste0("alpha.alpha2_TF_K-1.", seqLength)])/2
              beta1 = (pa_best[base::paste0("beta.beta1_TG_K-1.", seqLength)] + pa_best[base::paste0("beta.beta1_TF_K-1.", seqLength)])/2
              beta2 = (pa_best[base::paste0("beta.beta2_TG_K-1.", seqLength)] + pa_best[base::paste0("beta.beta2_TF_K-1.", seqLength)])/2
              beta3 = pa_best[base::paste0("beta.beta3_TG_K-1.", seqLength)]
              s1 = (pa_best["s.s1_TG"] + pa_best["s.s1_TF"])/2
              s2 = (pa_best["s.s2_TG"] + pa_best["s.s2_TF"])/2
              u11 = pa_best["u.u11_TG"]
              u21 = pa_best["u.u21_TG"]
              u311 = pa_best["u.u311_TG"]
              u321 = pa_best["u.u321_TG"]
              u331 = pa_best["u.u331_TG"]
              u12 = pa_best["u.u12_TG"]
              u22 = pa_best["u.u22_TG"]
              u312 = pa_best["u.u312_TG"]
              u322 = pa_best["u.u322_TG"]
              u332 = pa_best["u.u332_TG"]
              TGE_tilde = pa_best[base::paste0("E~.E~_TG_K-1.", seqLength)]
              TFE_tilde = pa_best[base::paste0("E~.E~_TF_K-1.", seqLength)]
              U1_tilde = scDEDS::Hill_cal(x = TGE_tilde, Dissociation_Constant = u11, Hill_Coefficient = u12)
              U2_tilde = scDEDS::Hill_cal(x = TFE_tilde, Dissociation_Constant = u21, Hill_Coefficient = u22)
              U1 = scDEDS::Hill_cal(x = TGE[-Length], Dissociation_Constant = u11, Hill_Coefficient = u12)
              U2 = scDEDS::Hill_cal(x = TFE[-Length], Dissociation_Constant = u21, Hill_Coefficient = u22)
              U31 = scDEDS::Hill_cal(x = TGA[-Length], Dissociation_Constant = u311, Hill_Coefficient = u312)
              U32 = scDEDS::Hill_cal(x = TGE[-Length], Dissociation_Constant = u321, Hill_Coefficient = u322)
              U33 = scDEDS::Hill_cal(x = TGE_tilde[-Length], Dissociation_Constant = u331, Hill_Coefficient = u332)
              v1 = pa_best["v.v1_TG"]
              v2 = pa_best["v.v2_TG"]
              v3 = pa_best["v.v3_TG"]
              base::sum(base::sapply(2:Length, function(K) {
                prev_idx = K - 1
                2 * (TFE[prev_idx] + alpha1[prev_idx] * R * scDEDS::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde[prev_idx]) + beta1[prev_idx] - TFE[K])^2 +
                  2 * (TGA[prev_idx] + alpha2[prev_idx] * R * scDEDS::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde[prev_idx]) + beta2[prev_idx] - TGA[K])^2 +
                  (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33[prev_idx]) * Tslot_K[K] + beta3[prev_idx] - TGE[K])^2
              }))
            }, lower = 0.01, upper = 0.99, method = "Brent")
            theta_p[re, re_tr] = predict_result$par
            loss[re, re_tr] = predict_result$value
          }
          loss = loss[!base::grepl(".", base::names(loss), fixed = TRUE)]
          theta_p = theta_p[!base::grepl(".", base::names(theta_p), fixed = TRUE)]

          pGRN[[re]] = base::list(theta_p = theta_p, loss = loss)
        }

        pGRN[["theta_p"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pGRN[4:base::length(pGRN)], function(x) x[["theta_p"]]
        )))
        base::rownames(pGRN[["theta_p"]]) = base::names(pGRN[4:base::length(pGRN)])

        pGRN[["loss"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pGRN[4:base::length(pGRN)], function(x) x[["loss"]]
        )))
        base::rownames(pGRN[["loss"]]) = base::names(pGRN[4:base::length(pGRN)])

        return(pGRN)
      }

      prediction_list = parallel::mclapply(
        X = 1:n_part,
        FUN = function(s) {process_data_slice(
          pGRN = base::list(theta_p = NA, loss = NA),
          RE_PRED = re_pred_list[[s]],
          ncols = base::length(re_train),
          STATE_DATA = state_data[re_pred_list[[s]], -(1:3)],
          TFE_col = TFE_col,
          TGA_col = TGA_col,
          TGE_col = TGE_col,
          Length = Length,
          seqLength = seqLength,
          re_train = re_train,
          train_params = train_params,
          tao = tao,
          Tslot_K = Tslot_K
        )},
        mc.cores = ncores
      )

      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]] = base::list()

      message("Summarizing theta_p of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]] =
        base::do.call("rbind", base::lapply(prediction_list, function(x) x[["theta_p"]]))

      message("Summarizing loss of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["loss"]] =
        base::do.call("rbind", base::lapply(prediction_list, function(x) x[["loss"]]))

      message("Saving the predtive results of this branch.")
      base::saveRDS(
        interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]],
        file = base::paste0("./interest_cell_type_pGRN.", base::gsub(" ", "_", cell_type), ".branch", n, ".rds")
      )
    }
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring GRN ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_pGRN)
}
