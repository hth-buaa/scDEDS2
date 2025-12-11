#' @title Numerical Gradient Calculation for Model Training Objective Function
#'
#' @description
#' This function computes the numerical gradient of the training objective function using finite difference approximation.
#' It employs parallel processing to efficiently calculate partial derivatives with respect to each model parameter,
#' supporting gradient-based optimization algorithms for model training.
#'
#' @param params See in ?f_train_cal.
#' @param h_max Numeric. Maximum step size for finite difference approximation. Default is 0.01.
#' @param h_min_percent Numeric. Minimum step size as a percentage of parameter range (between upper and lower bounds). Default is 0.01 (1%).
#' @param ncores
#' See in ?get_interest_cell_type_data.
#' While this code can be executed on Windows systems, it is still strongly recommended to run it on Linux systems, as the execution time without parallel computing would be astronomically long.
#' The actual runtime depends on the total available RAM, the number of CPU cores, the previously determined number of cell groups, and the number of transcription factors (TFs) and target genes (TGs) obtained.
#' Larger core counts, fewer cell groups, and fewer genes correspond to shorter runtime, while the total available RAM dictates the maximum number of cores that can participate in parallel computing.
#' @param params_names See in ?f_train_cal.
#' @param train_data See in ?f_train_cal.
#' @param tao See in ?f_train_cal.
#' @param Tslot_K See in ?f_train_cal.
#' @param params_upper Numeric vector. Upper bounds for each parameter, used to determine the maximum step size in finite difference approximation.
#' @param params_lower Lower bounds for each parameter, used to determine the minimum step size in finite difference approximation.
#'
#' @returns
#' A numeric vector representing the gradient of the objective function with respect to each parameter.
#' Using the finite difference method to approximate the gradient.
#'
#' @details
#' The function implements the following computational strategy:
#' \enumerate{
#' \item Step Size Determination: For each parameter, computes an adaptive step size that is the minimum of \code{h_max} and \code{h_min_percent} of the parameter's range (difference between upper and lower bounds)
#' \item Baseline Evaluation: Computes the objective function value at the current parameters
#' \item Parallel Gradient Calculation: Uses parallel processing to compute partial derivatives for each parameter by perturbing one parameter at a time
#' \item Robustness Handling: Replaces NA or infinite gradient values with 0 for numerical stability
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item Depends on \code{scDEDS::f_train_cal} for objective function calculation
#' \item Uses forward difference approximation
#' \item Parallel processing significantly speeds up computation for high-dimensional parameter spaces
#' \item Adaptive step sizing helps balance numerical precision and computational stability
#' \item NA and infinite values are replaced with 0, which may mask numerical issues but can ensure smooth program operation
#' }
#'
#' @seealso \code{\link{f_train_cal}}, \code{\link[parallel]{mclapply}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' h_max = 0.01
#' h_min_percent = 0.01
#' ncores = ceiling(detectCores() / 2) # In Linux.
#' # The settings for the other parameters are inside the function scDEDS::gradient_ascent_train.
#' temp_grad = f_train_cal_gr(
#'   params,
#'   h_max = h_max,
#'   h_min_percent = h_min_percent,
#'   ncores = ncores,
#'   params_names = params_names,
#'   train_data = train_data,
#'   tao = tao,
#'   Tslot_K = Tslot_K,
#'   params_upper = params_upper,
#'   params_lower = params_lower
#' )
#' }
f_train_cal_gr = function (params,
                           h_max = 0.01,
                           h_min_percent = 0.01,
                           ncores = 1,
                           params_names = params_names,
                           train_data = train_data,
                           tao = tao,
                           Tslot_K = Tslot_K,
                           params_upper = params_upper,
                           params_lower = params_lower)
{
  vec = base::as.numeric(params)
  n = base::length(vec)
  f0 = scDEDS::f_train_cal(
    params = vec,
    params_names = params_names,
    train_data = train_data,
    tao = tao,
    Tslot_K = Tslot_K
  )
  h_vec = base::pmin(base::rep(h_max, n), (params_upper - params_lower) * h_min_percent)
  grad = base::unlist(parallel::mclapply(1:n, function(i) {
    par_perturb = vec
    par_perturb[i] = vec[i] + h_vec[i]
    gr = (scDEDS::f_train_cal(
      params = par_perturb,
      params_names = params_names,
      train_data = train_data,
      tao = tao,
      Tslot_K = Tslot_K
    ) - f0)/h_vec[i]
    rm(list = base::setdiff(base::ls(), "gr"))
    return(gr)
  }, mc.cores = ncores, mc.silent = TRUE, mc.allow.recursive = FALSE))
  grad[base::is.na(grad) | base::is.infinite(grad)] = 0
  return(grad)
}
