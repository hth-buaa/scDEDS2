#' @title Model Training for Gene Regulatory Networks Prediction
#'
#' @description
#' This function trains prediction models for gene regulatory networks using a hybrid optimization approach combining genetic algorithms and gradient ascent.
#' It processes each cell type and branch independently, training individual models for each regulator-target pair (regulon).
#' The training incorporates multiple early stopping criteria and parallel processing capabilities for efficient computation.
#'
#' @param interest_cell_type_branch_training_set The output of function select_training_set.
#' @param interest_cell_type_branch_init_params The output of function set_init_params.
#' @param interest_cell_type_group The output of function cell_grouping.
#' @param interest_cell_type_tao The output of function get_cell_type_tao.
#' @param interest_cell_type_genes_pseudotime_info The output of function get_genes_pseudotime_info.
#' @param max_epochs Positive integer. Maximum number of training epochs (iterations). Default is 1000.
#' @param early_stop Positive integer. Number of consecutive epochs without significant improvement to trigger early stopping. Default is 10.
#' @param eps_theta Tiny positive numeric. Threshold for early stopping when predicted regulatory strength (theta_p) is close to observed strength (theta_i). Default is 1e-3.
#' @param eps_loss Tiny positive numeric. Threshold for early stopping when loss change becomes negligible. Default is 1e-5.
#' @param popSize_ga Positive integer. Population size for the genetic algorithm, see popSize in ?GA::ga. Default is 50.
#' @param maxiter_ga Positive integer. Maximum number of generations for the genetic algorithm per iteration, see maxiter in ?GA::ga. Default is 150.
#' @param pcrossover_ga Numeric between 0 and 1. Probability of crossover between chromosome pairs in genetic algorithm, see pcrossover in ?GA::ga. Default is 0.9.
#' @param pmutation_ga Numeric between 0 and 1. Probability of mutation in parent chromosomes in genetic algorithm, see more in ?GA::ga. Default is 0.5.
#' @param parallel_ga Logical. Whether to enable parallel computation for the genetic algorithm, see more in ?GA::ga. Default is FALSE.
#' @param seed_ga Integer. Random seed for reproducibility of genetic algorithm results, see more in ?GA::ga. Default is 123.
#' @param iterations_grad Positive integer. Number of gradient ascent steps per iteration. Default is 30.
#' @param ncores_grad Positive integer. Number of cores for parallel gradient computation. Default is 1.
#' @param theta_difference_threshold Positive numeric. Threshold for filtering results based on difference between predicted and observed regulatory strengths. Default is 0.1.
#' @param fit_loss_threshold Positive numeric. Threshold for filtering results based on fitting loss. Default is 0.1.
#' @param ncores See in ?get_interest_cell_type_data.
#'
#' @returns
#' A nested list organized by cell type and branch, containing:
#' \itemize{
#' \item For each branch: Detailed training results for each regulon, including:
#' \itemize{
#' \item \code{re}: Regulon name
#' \item \code{train_data}: Training data used
#' \item \code{Tslot_K}: Psudotime intervals
#' \item \code{tao}: Regulatory threshold τ. If regulatory strength θ > τ, a regulatory relationship is confirmed (positive sample); o therwise, it is considered a negative sample.
#' \item \code{f_value_train}: Objective function values
#' \item \code{params_data_frame}: Parameters of each training set
#' \item \code{loss_train}: Training loss history
#' \item \code{best_epoch}: Epoch with lowest loss
#' \item \code{best_params}: Optimized parameters
#' \item \code{best_pred_result_re}: Optimized fit train regution data.
#' \item \code{stop_cause}: Reason for training termination
#' \item \code{train_time}: Duration of training
#' \item \code{train_epoch}: Train epoch
#' }
#' \item \code{best_pred_result}: Comprehensive dataframe with predictions for all regulons across branches
#' \item \code{best_pred_result_filtered}: Filtered version of predictions based on quality thresholds
#' }
#'
#' @details
#' The training process follows a hybrid optimization strategy:
#' \enumerate{
#' \item \strong{Initialization}: Parameters are initialized based on \code{interest_cell_type_branch_init_params}
#' \item \strong{Genetic Algorithm Phase}: Explores parameter space using GA with specified crossover/mutation probabilities
#' \item \strong{Gradient Ascent Phase}: Refines GA solutions using gradient-based optimization
#' \item \strong{Early Stopping}: Monitors multiple criteria:
#' \itemize{
#' \item Maximum epochs reached
#' \item Overfitting detection (loss increases)
#' \item Predicted vs. observed regulatory strength convergence
#' \item Insignificant loss improvement
#' }
#' \item \strong{Result Integration}: Combines branch-specific predictions using weighted averaging based on cell counts
#' \item \strong{Quality Filtering}: Filters predictions based on accuracy and loss thresholds
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Using larger pcrossover_ga, pmutation_ga, and maxiter_ga to jump out of the local optimum
#' \item Training each regulon independently enables parallel processing across regulons
#' \item Early stopping parameters should be adjusted based on dataset size and complexity
#' \item Higher crossover/mutation probabilities increase exploration but may slow convergence
#' \item Memory requirements scale with the number of regulons and parameter dimensions
#' }
#'
#' @seealso \code{\link[GA]{ga}}, \code{\link[scDEDS]{f_train_cal}}, \code{\link[scDEDS]{loss_cal}}, \code{\link[scDEDS]{theta_p_cal}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_training_set = base::readRDS("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")
#' interest_cell_type_branch_init_params = base::readRDS("./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds")
#' interest_cell_type_group = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_group.rds")
#' interest_cell_type_tao = base::readRDS("./3 get iGRN/interest_cell_type_tao.rds")
#' interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
#'
#' max_epochs = 50
#' early_stop = 5
#'
#' popSize_ga = 50
#' maxiter_ga = 150
#' pcrossover_ga = 0.9
#' pmutation_ga = 0.5
#' parallel_ga = FALSE
#' seed_ga = 123 # Reproducible
#' # seed_ga = stats::runif(1) # Not Reproducible
#'
#' iterations_grad = 30
#' ncores_grad = 1
#'
#' theta_difference_threshold = 0.1
#' fit_loss_threshold = 0.1
#'
#' ncores = parallel::detectCores() - 1 # in Linux
#' # ncores = 1 # in Windows
#'
#' interest_cell_type_branch_model_train = model_train(
#'   interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
#'   interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
#'   interest_cell_type_group = interest_cell_type_group,
#'   interest_cell_type_tao = interest_cell_type_tao,
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   max_epochs = max_epochs, early_stop = early_stop,
#'   popSize_ga = popSize_ga, maxiter_ga = maxiter_ga, pcrossover_ga = pcrossover_ga, pmutation_ga = pmutation_ga, parallel_ga = parallel_ga, seed_ga = 123,
#'   iterations_grad = iterations_grad, ncores_grad = ncores_grad,
#'   theta_difference_threshold = theta_difference_threshold, fit_loss_threshold = fit_loss_threshold, ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_branch_model_train, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_model_train.rds")
#' }
model_train = function (
    interest_cell_type_branch_training_set = interest_cell_type_branch_training_set,
    interest_cell_type_branch_init_params = interest_cell_type_branch_init_params,
    interest_cell_type_group = interest_cell_type_group,
    interest_cell_type_tao = interest_cell_type_tao,
    interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
    max_epochs = 1000, early_stop = 10,
    eps_theta = 1e-3, eps_loss = 1e-5,
    popSize_ga = 50, maxiter_ga = 150, pcrossover_ga = 0.9, pmutation_ga = 0.5, parallel_ga = FALSE, seed_ga = 123,
    iterations_grad = 30, ncores_grad = 1,
    theta_difference_threshold = 0.1, fit_loss_threshold = 0.1, ncores = 1
)
{
  t_start = base::Sys.time()
  message("Run: Training models for each branch ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "4.2 BUild Prediction Model"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  interest_cell_type_branch_model_train = base::list()
  for (cell_type in base::names(interest_cell_type_branch_training_set)) {
    message("Starting training the prediction models for each branch in cell type ", cell_type, ".")

    message("Initializing model training iteration result storage.")
    interest_cell_type_branch_model_train[[cell_type]] = base::list()

    message("Initializing optimal prediction result storage.")
    best_pred_result = interest_cell_type_branch_training_set[[cell_type]][["all"]]
    best_pred_result["theta_p"] = 0
    best_pred_result[base::paste0("theta_p_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = 0
    best_pred_result$theta_i_bin = base::as.integer(best_pred_result$theta_i > 0)
    best_pred_result["theta_p_bin"] = 0
    best_pred_result[base::paste0("theta_p_bin_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = 0
    best_pred_result["fit_loss"] = Inf
    best_pred_result[base::paste0("fit_loss_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = Inf
    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = best_pred_result

    message("Iterating through each branch and training each branch independently.")
    for (n in 1:(base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1)) {
      message("Starting training the prediction model for branch ", n, " in cell type ", cell_type, ".")
      message("There are ", nrow(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]), " regulons for training.")
      n_re = 0
      interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] = base::list()

      # Train each regulon.
      message("Training each regulon.")
      model_train_process = function(re,
                                     n = n,
                                     n_re = base::match(re, base::names(interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]])),
                                     N_re = base::nrow(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]),
                                     N_branch = base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1,
                                     branch_params_list = interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]][[re]],
                                     train_data = base::as.matrix(interest_cell_type_branch_training_set[[cell_type]][[n]][re, , drop = FALSE]),
                                     Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]],
                                     tao = interest_cell_type_tao[cell_type],
                                     popSize_ga = popSize_ga,
                                     max_epochs = max_epochs,
                                     pcrossover_ga = pcrossover_ga,
                                     pmutation_ga = pmutation_ga,
                                     maxiter_ga = maxiter_ga,
                                     parallel_ga = parallel_ga,
                                     seed_ga = seed_ga,
                                     best_pred_result_re = best_pred_result[re, , drop = FALSE]) {
        result = base::list()
        result[["re"]] = re

        signal = base::paste0(n_re, "/", N_re)
        message("Progress ", signal, "...")
        Signal = base::paste0("[Branch", n, "/", N_branch, "-TrainRegulon", signal, "] ")
        message(Signal, "Extracting initial parameter values, upper bounds, lower bounds, and names.")
        params = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$initialization)
        }))
        params_upper = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$upper)
        }))
        params_lower = base::unlist(base::lapply(branch_params_list, function(category) {
          base::lapply(category, function(param) param$lower)
        }))
        params_names = base::names(params)

        # message(Signal, "Extracting training set, training set, and test set data.")
        message(Signal, "Extracting training set.")
        result[["train_data"]] = train_data

        message(Signal, "Extracting the pseudotime interval length for all cell groups.")
        result[["Tslot_K"]] = Tslot_K

        message(Signal, "Extracting the binary classification threshold for regulatory strength.")
        result[["tao"]] = tao

        message(Signal, "Initializing the objective function value in model training.")
        f_value_train = 1e+5 + 1

        message(Signal, "Initializing the genetic algorithm population.")
        prev_population = base::matrix(
          data = base::rep(base::as.numeric(params), popSize_ga),
          nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
        )

        message(Signal, "Initializing the model parameters storage for each training round.")
        para_data_frame = base::as.data.frame(base::matrix(
          data = 0, nrow = max_epochs, ncol = base::length(params),
          dimnames = base::list(1:max_epochs, params_names)
        ))

        message(Signal, "Initializing the loss storage for each training round.")
        loss_train = base::rep(0, max_epochs)
        names(loss_train) = 1:max_epochs
        f_value_train_vec = base::rep(0, max_epochs)
        names(f_value_train_vec) = 1:max_epochs

        t1_re = base::Sys.time()
        for(epoch in 1:max_epochs) {
          SIGNAL = base::paste0(base::gsub("] ", "", base::paste0(Signal, "-Epoch", epoch), fixed = TRUE), "] ")
          if (epoch %% 10 == 1) {message(SIGNAL, "Training Model.")}

          # Genetic Algorithm Training
          result_ga = GA::ga(
            type = "real-valued", fitness = scDEDS::f_train_cal, lower = params_lower, upper = params_upper,
            popSize = popSize_ga, pcrossover = pcrossover_ga, pmutation = pmutation_ga,
            elitism = base::ceiling(base::min(0.2 * maxiter_ga, popSize_ga / 5)),
            maxiter = maxiter_ga, suggestions = prev_population, parallel = parallel_ga, seed = seed_ga, monitor = FALSE,
            params_names = params_names, train_data = train_data, tao = tao, Tslot_K = Tslot_K
          )
          if (parallel_ga) {
            prev_population = base::matrix(
              data = base::rep(base::as.numeric(result_ga@solution[1, ]), popSize_ga),
              nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
            )
            result_ga = GA::ga(
              type = "real-valued", fitness = scDEDS::f_train_cal, lower = params_lower, upper = params_upper,
              popSize = popSize_ga, elitism = base::ceiling(base::min(0.2 * maxiter_ga, popSize_ga / 5)),
              maxiter = 1, suggestions = prev_population, parallel = FALSE, seed = seed_ga, monitor = FALSE,
              params_names = params_names, train_data = train_data, tao = tao, Tslot_K = Tslot_K
            ); base::gc()
          }
          f_value_train = -result_ga@fitnessValue

          # Gradient Descent Method Training
          result_gradient = gradient_ascent_train(
            iterations = iterations_grad, alpha_lower = 1e-5, alpha_upper = 1e-3, alpha_guess = 1e-4,
            init_theta = base::as.numeric(result_ga@solution[1, ]),
            h_max = 1e-9, h_min_percent = 0.01, ncores = ncores_grad,
            params_lower = params_lower, params_upper = params_upper, monitor = FALSE,
            params_names = params_names, train_data = train_data, tao = tao, Tslot_K = Tslot_K
          )
          f_value_train = -result_gradient$objective_value
          f_value_train_vec[epoch] = f_value_train
          result[["f_value_train"]] = f_value_train_vec
          prev_population = base::matrix(
            data = base::rep(result_gradient$parameters, popSize_ga),
            nrow = popSize_ga, ncol = base::length(params), byrow = TRUE
          )
          para_data_frame[epoch, ] = result_gradient$parameters

          result[["para_data_frame"]] = para_data_frame

          # Calcalating loss.
          loss_train[epoch] = scDEDS::loss_cal(
            data = train_data,
            params = result_gradient$parameters, params_names = params_names,
            tao = tao, Tslot_K = Tslot_K
          )

          # Early Stopping Strategy
          if (epoch >= max_epochs) {
            message(SIGNAL, "Reaching the maximum number of training epochs.")
            loss_vec = loss_train[1:base::ifelse(
              base::length(base::which(loss_train != 0)) > 0, base::max(base::which(loss_train != 0)), 0
            )]
            result[["loss_train"]] = loss_vec

            best_epoch = base::as.numeric(base::names(base::which.min(loss_vec)))
            result[["best_epoch"]] = best_epoch

            best_params = base::as.numeric(para_data_frame[best_epoch, , drop = FALSE])
            base::names(best_params) = params_names
            result[["best_params"]] = best_params

            best_pred_result_re[re, base::paste0("fit_loss_branch", n)] = loss_vec[best_epoch]
            best_pred_result_re[re, base::paste0("theta_p_branch", n)] = scDEDS::theta_p_cal(
              data = train_data,
              best_params = best_params, params_names = params_names,
              tao = tao, Tslot_K = Tslot_K
            )
            best_pred_result_re[re, base::paste0("theta_p_bin_branch", n)] = base::ifelse(best_pred_result_re[re, base::paste0("theta_p_branch", n)] >= tao, 1, 0)
            result[["best_pred_result_re"]] = best_pred_result_re

            stop_cause = "maximum training epochs reached."
            message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
            result[["stop_cause"]] = stop_cause
            break
          } else if (epoch >= early_stop) {
            loss_vec = loss_train[1:base::ifelse(
              base::length(base::which(loss_train != 0)) > 0, base::max(base::which(loss_train != 0)), 0
            )]
            result[["loss_train"]] = loss_vec

            best_epoch = base::as.numeric(base::names(base::which.min(loss_vec)))
            result[["best_epoch"]] = best_epoch

            best_params = base::as.numeric(para_data_frame[best_epoch, , drop = FALSE])
            base::names(best_params) = params_names
            result[["best_params"]] = best_params

            best_pred_result_re[re, base::paste0("fit_loss_branch", n)] = loss_vec[best_epoch]
            best_pred_result_re[re, base::paste0("theta_p_branch", n)] = scDEDS::theta_p_cal(
              data = train_data,
              best_params = best_params, params_names = params_names,
              tao = tao, Tslot_K = Tslot_K
            )
            best_pred_result_re[re, base::paste0("theta_p_bin_branch", n)] = base::ifelse(best_pred_result_re[re, base::paste0("theta_p_branch", n)] >= tao, 1, 0)
            result[["best_pred_result_re"]] = best_pred_result_re

            if (best_epoch <= base::length(loss_vec) - early_stop + 1) {
              stop_cause = "overfit"
              message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            } else if (base::abs(best_pred_result_re[re, base::paste0("theta_p_branch", n)] - best_pred_result_re[re, "theta_i"]) < eps_theta) {
              stop_cause = "theta_p_branch is close to theta_i"
              message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            } else if (base::max(loss_vec[(base::length(loss_vec) - early_stop + 1):base::length(loss_vec)]) - loss_vec[best_epoch] < eps_loss) {
              stop_cause = "loss change insignificant"
              message(SIGNAL, "Meeting the early stop condition: ", stop_cause, ".")
              result[["stop_cause"]] = stop_cause
              break
            }
          }
        }
        t2_re = base::Sys.time(); base::print(t2_re - t1_re)

        result[["train_time"]] = base::sprintf("%.2f %s", base::as.numeric(t2_re - t1_re), base::units(t2_re - t1_re))
        result[["train_epoch"]] = epoch

        return(result)
      }

      interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]] = parallel::mclapply(
        X = base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]),
        FUN = function(re) {
          model_train_process(
            re,
            n = n,
            n_re = base::match(re, base::names(interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]])),
            N_re = base::nrow(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]]),
            N_branch = base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1,
            branch_params_list = interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]][[re]],
            train_data = base::as.matrix(interest_cell_type_branch_training_set[[cell_type]][[n]][re, , drop = FALSE]),
            Tslot_K = interest_cell_type_group[[cell_type]][["Branches_Tslot_k"]][[n]],
            tao = interest_cell_type_tao[cell_type],
            popSize_ga = popSize_ga,
            max_epochs = max_epochs,
            pcrossover_ga = pcrossover_ga,
            pmutation_ga = pmutation_ga,
            maxiter_ga = maxiter_ga,
            parallel_ga = parallel_ga,
            seed_ga = seed_ga,
            best_pred_result_re = best_pred_result[re, , drop = FALSE]
          )
        },
        mc.cores = ncores
      )
      base::names(interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]]) =
        base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])

      for (re in base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])) {
        best_pred_result[re, base::paste0("theta_p_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]][re, base::paste0("theta_p_branch", n)]
        best_pred_result[re, base::paste0("theta_p_bin_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]][re, base::paste0("theta_p_bin_branch", n)]
        best_pred_result[re, base::paste0("fit_loss_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_pred_result_re"]][re, base::paste0("fit_loss_branch", n)]
      }
    }
    message("All branch trainings have been completed.")

    message("Integrating the results from all branches.")
    best_pred_result$theta_p_bin = base::as.integer(base::rowSums(
      best_pred_result[, base::grep("^theta_p_bin_branch\\d+$", base::names(best_pred_result))]) > 0)

    weights = base::unlist(interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_n_cell"]])
    best_pred_result$theta_p = base::apply(best_pred_result[, base::grep("^theta_p_branch([0-9]{1,2})$", base::names(best_pred_result), value = TRUE)], 1, function(x) {
      non_zero_idx = which(x != 0)
      if (base::length(non_zero_idx) == 0) 0 else {
        base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
      }
    })
    best_pred_result$fit_loss = base::apply(best_pred_result[, base::grep("^fit_loss_branch([0-9]{1,2})$", base::names(best_pred_result), value = TRUE)], 1, function(x) {
      non_zero_idx = which(x != Inf)
      if (base::length(non_zero_idx) == 0) Inf else {
        base::sum(x[non_zero_idx] * weights[non_zero_idx]) / base::sum(weights[non_zero_idx])
      }
    })

    best_pred_result = base::cbind(
      best_pred_result[, 1:4],
      base::abs(best_pred_result$theta_i - best_pred_result$theta_p),
      best_pred_result[, 5:base::ncol(best_pred_result)]
    )
    base::colnames(best_pred_result)[5] = "theta_difference"

    best_pred_result[base::paste0("theta_difference_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = Inf
    best_pred_result[base::paste0("train_time_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = 0
    best_pred_result[base::paste0("train_epoch_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = 0
    best_pred_result[base::paste0("best_epoch_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = 0
    best_pred_result[base::paste0("stop_cause_branch", 1:(base::length(interest_cell_type_branch_training_set[[cell_type]])-1))] = NA
    for (n in 1:(base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1)) {
      for (re in base::rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])) {
        best_pred_result[re, base::paste0("theta_difference_branch", n)] =
          base::abs(best_pred_result[re, "theta_i"] - best_pred_result[re, base::paste0("theta_p_branch", n)])
        best_pred_result[re, base::paste0("train_time_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["train_time"]]
        best_pred_result[re, base::paste0("train_epoch_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["train_epoch"]]
        best_pred_result[re, base::paste0("best_epoch_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["best_epoch"]]
        best_pred_result[re, base::paste0("stop_cause_branch", n)] =
          interest_cell_type_branch_model_train[[cell_type]][[base::paste0("branch", n)]][[re]][["stop_cause"]]
      }
    }

    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result"]] = best_pred_result

    message("Filtering training results.")
    interest_cell_type_branch_model_train[[cell_type]][["best_pred_result_filtered"]] = best_pred_result[
      best_pred_result$theta_difference < theta_difference_threshold & best_pred_result$fit_loss < fit_loss_threshold,
    ]
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Training models for each branch ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_branch_model_train)
}
