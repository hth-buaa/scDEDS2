#' @title This function is the minimized objective function in the model trained using the training set data.
#'
#' @description
#' Computes the objective function value to be minimized during model training.
#' The function evaluates the model fit by calculating weighted prediction errors for gene expression dynamics (TFE, TGA, TGE) and includes regularization terms to prevent overfitting.
#' The output is scaled and negated for compatibility with standard minimization algorithms.
#'
#' @param params Vector. A flattened vector containing all model parameters for a specific branch, with parameter names defined before.
#' @param params_names Character vector. Names corresponding to each element in the params vector, used to assign parameters within the function.
#' @param train_data Data frame. Training dataset containing only one regulatory pair (TF-TG) with columns for theta_i(regulatory strength) and time-series values for TFE, TGA, and TGE features.
#' @param tao See in ?R_cal.
#' @param Tslot_K Numeric vector. Psudotime intervals between consecutive measurements, used in the TGE dynamics calculation.
#'
#' @returns
#' A numeric value representing the negative weighted objective function (multiplied by -10000 for optimization purposes).
#' Lower values indicate better model fit. The objective function combines:
#' \itemize{
#' \item Mean squared prediction errors for TFE, TGA, and TGE values across time points
#' \item Variance regularization terms for beta parameters and expression states
#' \item Weighting factors that emphasize different components of the model
#' }
#'
#' @details
#' The function performs the following complex operations for each TG-TF pair in the training data:
#' \enumerate{
#' \item Parameter Organization: Assigns names to the parameter vector
#' \item TG-TF Pair Processing: Extracts and processes each regulator-target pair from training data
#' \item Objective Computation: Calculates weighted squared errors between predicted and observed values The mathematical formulation incorporates, where the error terms include squared differences for TFE, TGA, and TGE predictions, and regularization terms include variances of beta parameters and expression states.
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item The function returns negative values for compatibility with minimization algorithms
#' \item The function is computationally intensive due to multiple nested operations
#' \item Regularization terms help prevent overfitting by penalizing parameter variations
#' \item The weighting scheme (-10000 times) is designed for optimization stability
#' }
#'
#' @seealso \code{\link{R_cal}}, \code{\link{Hill_cal}}, \code{\link{S_cal}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # The settings for parameters are inside the function scDEDS::f_train_cal_gr, scDEDS::gradient_ascent_train, or scDEDS::gradient_ascent_train model_train.
#' f0 = f_train_cal(
#'   params,
#'   params_names = params_names,
#'   train_data = train_data,
#'   tao = tao,
#'   Tslot_K = Tslot_K
#' )
#' }
f_train_cal = function (params,
                        params_names = params_names,
                        train_data = train_data,
                        tao = tao,
                        Tslot_K = Tslot_K)
{
  base::names(params) = params_names
  TFE = base::as.numeric(train_data[, base::grep("TFE_", base::colnames(train_data), value = TRUE)])
  TGA = base::as.numeric(train_data[, base::grep("TGA_", base::colnames(train_data), value = TRUE)])
  TGE = base::as.numeric(train_data[, base::grep("TGE_", base::colnames(train_data), value = TRUE)])
  Length = base::length(TFE)
  seqLength =base::seq_len(Length - 1)
  R = base::as.numeric(scDEDS::R_cal(
    theta_TF_TG = base::as.numeric(train_data[, "theta_i"]),
    tao = tao,
    r1 = 200,
    r2 = (params["r.r2_TG"] + params["r.r2_TF"])/2,
    r3 = 1,
    r4 = 1000,
    r5 = 1
  ))
  alpha1 = (params[base::paste0("alpha.alpha1_TG_K-1.", seqLength)] + params[base::paste0("alpha.alpha1_TF_K-1.", seqLength)])/2
  alpha2 = (params[base::paste0("alpha.alpha2_TG_K-1.", seqLength)] + params[base::paste0("alpha.alpha2_TF_K-1.", seqLength)])/2
  beta1 = (params[base::paste0("beta.beta1_TG_K-1.", seqLength)] + params[base::paste0("beta.beta1_TF_K-1.", seqLength)])/2
  beta2 = (params[base::paste0("beta.beta2_TG_K-1.", seqLength)] + params[base::paste0("beta.beta2_TF_K-1.", seqLength)])/2
  beta3 = params[base::paste0("beta.beta3_TG_K-1.", seqLength)]
  s1 = (params["s.s1_TG"] + params["s.s1_TF"])/2
  s2 = (params["s.s2_TG"] + params["s.s2_TF"])/2
  u11 = params["u.u11_TG"]
  u21 = params["u.u21_TG"]
  u311 = params["u.u311_TG"]
  u321 = params["u.u321_TG"]
  u331 = params["u.u331_TG"]
  u12 = params["u.u12_TG"]
  u22 = params["u.u22_TG"]
  u312 = params["u.u312_TG"]
  u322 = params["u.u322_TG"]
  u332 = params["u.u332_TG"]
  TGE_tilde = params[base::paste0("E~.E~_TG_K-1.", seqLength)]
  TFE_tilde = params[base::paste0("E~.E~_TF_K-1.", seqLength)]
  U1_tilde = scDEDS::Hill_cal(x = TGE_tilde, Dissociation_Constant = u11, Hill_Coefficient = u12)
  U2_tilde = scDEDS::Hill_cal(x = TFE_tilde, Dissociation_Constant = u21, Hill_Coefficient = u22)
  U1 = scDEDS::Hill_cal(x = TGE[-Length], Dissociation_Constant = u11, Hill_Coefficient = u12)
  U2 = scDEDS::Hill_cal(x = TFE[-Length], Dissociation_Constant = u21, Hill_Coefficient = u22)
  U31 = scDEDS::Hill_cal(x = TGA[-Length], Dissociation_Constant = u311, Hill_Coefficient = u312)
  U32 = scDEDS::Hill_cal(x = TGE[-Length], Dissociation_Constant = u321, Hill_Coefficient = u322)
  U33 = scDEDS::Hill_cal(x = TGE_tilde[-Length], Dissociation_Constant = u331, Hill_Coefficient = u332)
  v1 = params["v.v1_TG"]
  v2 = params["v.v2_TG"]
  v3 = params["v.v3_TG"]
  -10000 * (base::mean(base::sapply(2:Length, function(K) {
    prev_idx = K - 1
    2 * (TFE[prev_idx] + alpha1[prev_idx] * R * scDEDS::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde[prev_idx]) + beta1[prev_idx] - TFE[K])^2 +
      2 * (TGA[prev_idx] + alpha2[prev_idx] * R * scDEDS::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde[prev_idx]) + beta2[prev_idx] - TGA[K])^2 +
      (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33[prev_idx]) * Tslot_K[K] + beta3[prev_idx] - TGE[K])^2
  })) + 1 * (base::mean(beta1^2) + base::mean(beta2^2) + base::mean(beta3^2) + 1 * stats::var(TGE_tilde) + 1 * stats::var(TFE_tilde)))
}
