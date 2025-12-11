#' @title Select Best Predictions for Gene Regulatory Networks
#'
#' @description
#' Filters and selects the best predictions from inferred gene regulatory networks (GRN) by applying multiple quality thresholds.
#' The function processes each cell type and branch independently, selects the optimal regulatory strength (theta_p) and fitting loss for each TF-TG pair, and integrates results across branches using weighted averaging based on cell counts.
#'
#' @param interest_cell_type_pGRN The output of function infer_GRN.
#' @param interest_cell_type_branch_model_train The output of model_train.
#' @param interest_cell_type_tao The output of function get_cell_type_tao.
#' @param interest_cell_type_genes_pseudotime_info The output of function get_genes_pseudotime_info.
#' @param fit_loss_quantile_threshold Numeric between 0 and 1. Quantile threshold for filtering predictions based on loss. Default is 0.2 (20th percentile).
#' @param theta_p_quantile_threshold Numeric between 0 and 1. Quantile threshold for filtering predictions based on regulatory strength. Default is 0.8 (80th percentile).
#' @param fit_loss_upper Positive numeric. Upper bound for acceptable fitting loss. Default is 0.2.
#' @param theta_p_lower Numeric between 0 and 1. Lower bound for acceptable regulatory strength. Default is 0.2.
#' @param each See in ?infer_GRN.
#' @param ncores See in ?get_interest_cell_type_data.
#'
#' @returns
#' The input list \code{interest_cell_type_pGRN} with two additional elements for each cell type:
#' \itemize{
#' \item For each branch: A data frame \code{best_result} containing the best theta_p, loss, and the model used for each TF-TG pair.
#' \item A \code{pred_summary} data frame that integrates results across branches, with weighted averages for theta_p and fit_loss.
#' }
#'
#' @details
#' The selection process follows a multi-step strategy:
#' \enumerate{
#' \item \strong{Threshold Calculation}: Determines loss and theta_p thresholds using quantiles and fixed bounds.
#' \item \strong{Parallel Processing}: Partitions data for efficient parallel computation.
#' \item \strong{Best Prediction Selection}: For each TF-TG pair, selects the best prediction based on a hierarchy of criteria:
#' \itemize{
#' \item Uses predictions with loss below the training loss maximum to find the best prediction.
#' \item If fails, uses predictions with loss below the calculated loss threshold to find the best prediction.
#' \item If fails, uses predictions with theta_p above the calculated theta_p threshold to find the best prediction.
#' \item If fails, selects the prediction with the minimum loss to find the best prediction.
#' }
#' \item \strong{Result Integration}: Combines branch-specific results using weighted averaging based on cell counts.
#' \item \strong{Quality Control}: Removes predictions with theta_p outside the range from 0.011 to 0.989.
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Parallel processing speeds up the selection, especially for large datasets.
#' \item The \code{each} parameter controls memory usage; reduce if memory issues occur.
#' \item The function modifies the input list by adding new elements; original data is preserved.
#' \item Adjust quantile thresholds to balance sensitivity and specificity.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' interest_cell_type_branch_model_train = base::readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
#' interest_cell_type_tao = base::readRDS("./3 get iGRN/interest_cell_type_tao.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
#' fit_loss_quantile_threshold1 = 0.2
#' theta_p_quantile_threshold1 = 0.8
#' fit_loss_upper1 = 0.2
#' theta_p_lower1 = 0.2
#' each = 1000
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_pGRN = select_best_prediction(
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
#'   interest_cell_type_tao = interest_cell_type_tao,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   fit_loss_quantile_threshold = fit_loss_quantile_threshold1,
#'   theta_p_quantile_threshold = theta_p_quantile_threshold,
#'   fit_loss_upper = fit_loss_upper1,
#'   theta_p_lower = theta_p_lower1,
#'   each = each,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_pGRN, file = "./5 Infer GRN/interest_cell_type_pGRN.rds")
#' }
select_best_prediction = function(
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    interest_cell_type_branch_model_train = interest_cell_type_branch_model_train,
    interest_cell_type_tao = interest_cell_type_tao,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    fit_loss_quantile_threshold = 0.2,
    theta_p_quantile_threshold = 0.8,
    fit_loss_upper = 0.2,
    theta_p_lower = 0.2,
    each = 1000000,
    ncores = 1
)
{
  t_start = base::Sys.time()
  message("Run: Filtering GRN ", t_start, ".")
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

  for (cell_type in base::names(interest_cell_type_pGRN)) {
    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = list()
    train_loss_max = base::max(interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]]$fit_loss)

    for (n in 1:(base::length(interest_cell_type_branch_model_train[[cell_type]]) - 2)) {
      message("Process branch ", n, " in cell type ", cell_type, ".")
      re_pred = base::rownames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]])
      N_pred_all = base::length(re_pred)
      message("There are ", N_pred_all, " potential regulatory relationships in branch ", n, " of cell type ", cell_type, ".")

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
        pGRN_list[[s]][["theta_p"]] = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["theta_p"]][re_pred_list[[s]], ]
        pGRN_list[[s]][["loss"]] = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["loss"]][re_pred_list[[s]], ]
        pGRN_list[[s]][["best_result"]] = base::data.frame(
          theta_p = base::rep(0, base::length(re_pred_list[[s]])),
          loss = base::rep(Inf, base::length(re_pred_list[[s]])),
          fit = base::rep(NA, base::length(re_pred_list[[s]])),
          row.names = re_pred_list[[s]]
        )
      }

      message("Selecting best theta for each TF-TG of branch ", n, " of cell type ", cell_type, ".")
      loss_max = base::min(
        stats::quantile(
          base::unlist(interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][["loss"]]),
          probs = fit_loss_quantile_threshold,
          na.rm = TRUE
        ),
        fit_loss_upper
      )
      theta_p_min = base::max(
        stats::quantile(
          base::unlist(interest_cell_type_pGRN[[cell_type]][[paste0("branch", n)]][["theta_p"]]),
          probs = theta_p_quantile_threshold,
          na.rm = TRUE
        ),
        interest_cell_type_tao[cell_type],
        theta_p_lower
      )

      process_data_slice = function(pGRN = pGRN_list[[s]],
                                    train_loss_max = train_loss_max,
                                    loss_max = loss_max,
                                    theta_p_min = theta_p_min) {
        for (re in base::rownames(pGRN[[1]])) {
          theta_p = pGRN[["theta_p"]][re, , drop = FALSE]
          loss = pGRN[["loss"]][re, , drop = FALSE]
          best_result = base::data.frame(base::matrix(NA, nrow = 1, ncol = 3, dimnames = base::list(re, c("theta_p", "loss", "fit"))))
          best_result[re, 1:2] = c(0, Inf)

          theta_p_selected1 = base::max(theta_p[, loss < train_loss_max, drop = FALSE])
          theta_p_selected2 = base::max(theta_p[, loss < loss_max, drop = FALSE])
          theta_p_selected3 = theta_p[, theta_p > theta_p_min, drop = FALSE]
          if (theta_p_selected1 != -Inf) {
            best_result[re, "theta_p"] = theta_p_selected1
            best_result[re, "fit"] = base::names(loss)[which(theta_p == theta_p_selected1)][1]
            best_result[re, "loss"] = loss[re, best_result[re, "fit"]]
          } else if (theta_p_selected2 != -Inf) {
            best_result[re, "theta_p"] = theta_p_selected2
            best_result[re, "fit"] = base::names(loss)[which(theta_p == theta_p_selected2)][1]
            best_result[re, "loss"] = loss[re, best_result[re, "fit"]]
          } else if (base::ncol(theta_p_selected3) > 0) {
            best_result[re, "loss"] = base::min(loss[base::names(theta_p_selected3)])
            best_result[re, "fit"] = base::names(loss)[base::which(loss == best_result[re, "loss"])][1]
            best_result[re, "theta_p"] = theta_p[re, best_result[re, "fit"]]
          } else {
            best_result[re, "loss"] = base::min(loss)
            best_result[re, "fit"] = base::names(loss)[base::which(loss == best_result[re, "loss"])][1]
            best_result[re, "theta_p"] = theta_p[re, best_result[re, "fit"]]
          }

          pGRN[[re]] = base::list(best_result = best_result)
        }

        pGRN[["best_result"]] = base::as.data.frame(data.table::rbindlist(base::lapply(
          pGRN[4:base::length(pGRN)], function(x) x[["best_result"]]
        )))
        base::rownames(pGRN[["best_result"]]) = base::names(pGRN[4:base::length(pGRN)])
        pGRN[["best_result"]] = pGRN[["best_result"]][
          pGRN[["best_result"]]$theta_p >= 0.011 &
            pGRN[["best_result"]]$theta_p <= 0.989 &
            !base::is.infinite(pGRN[["best_result"]]$loss), , drop = FALSE]

        return(pGRN[["best_result"]])
      }

      prediction_list = parallel::mclapply(
        X = 1:n_part,
        FUN = function(s) {process_data_slice(
          pGRN = pGRN_list[[s]],
          train_loss_max = train_loss_max,
          loss_max = loss_max,
          theta_p_min = theta_p_min
        )},
        mc.cores = ncores
      )

      message("Summarizing best_result of this branch.")
      interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]] =
        base::do.call("rbind", prediction_list)
      base::colnames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]) = base::paste0(
        base::colnames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]), "_", "branch", n
      )
    }

    message("Summarizing the predtive results of each branch.")
    rows = base::sort(base::unique(base::Reduce(base::c, base::lapply(
      interest_cell_type_pGRN[[cell_type]][base::paste0("branch", 1:(base::length(interest_cell_type_branch_model_train[[cell_type]]) - 2))],
      function(x) base::rownames(x[["best_result"]])
    ))))
    cols = base::Reduce(base::c, base::lapply(
      interest_cell_type_pGRN[[cell_type]][base::paste0("branch", 1:(base::length(interest_cell_type_branch_model_train[[cell_type]]) - 2))],
      function(x) base::colnames(x[["best_result"]])
    ))
    pred_summary = base::data.frame(base::matrix(NA, nrow = base::length(rows), ncol = base::length(cols), dimnames = base::list(rows, cols)))
    for (n in 1:(base::length(interest_cell_type_branch_model_train[[cell_type]]) - 2)) {
      pred_summary[
        base::rownames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]),
        base::colnames(interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]])
      ] = interest_cell_type_pGRN[[cell_type]][[base::paste0("branch", n)]][["best_result"]]
    }

    weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
    pred_summary$theta_p = base::apply(pred_summary[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(pred_summary), value = TRUE)], 1, function(x) {
      non_zero_idx = which(!base::is.na(x))
      if (base::length(non_zero_idx) == 0) 0 else {
        base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
      }
    })
    pred_summary$fit_loss = base::apply(pred_summary[, base::grep("^loss_branch([0-9]{1,2})$", base::names(pred_summary), value = TRUE)], 1, function(x) {
      non_zero_idx = which(!base::is.na(x))
      if (base::length(non_zero_idx) == 0) Inf else {
        base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
      }
    })
    interest_cell_type_pGRN[[cell_type]][["pred_summary"]] = pred_summary
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Filtering GRN ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_pGRN)
}
