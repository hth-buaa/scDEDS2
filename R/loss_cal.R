#' @title Calculate Prediction Loss for a Specific TG-TF Pair
#'
#' @description
#' This function computes the prediction loss (mean squared error) for a specific target gene-transcription factor (TG-TF) pair using optimized model parameters.
#' It evaluates how well the trained model predicts the expression dynamics of TFE, TGA, and TGE across multiple time points.
#'
#' @param data
#' Data frame. The dataset containing expression values for one TG-TF pair,
#' typically including columns for TFE, TGA, TGE across multiple time points, and the regulatory strength theta_s.
#' Typically choose one from three options: training set data, validation set data, or test set data.
#' @param params See in ?f_train_cal.
#' @param params_names See in ?f_train_cal.
#' @param tao See in ?f_train_cal.
#' @param Tslot_K See in ?f_train_cal.
#'
#' @returns
#' A numeric value representing the mean squared prediction error across all time points for the specified TG-TF pair.
#' The loss combines prediction errors for the evolutino process of TFE (Transcription Factor Expression), TGA (Target Gene Accessibility), TGE (Target Gene Expression).
#' Lower values indicate better prediction performance.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#' \item Data Extraction: Extracts TFE, TGA, and TGE expression values across time points
#' \item Parameter Retrieval: Retrieves all relevant model parameters for the specific TG-TF pair
#' \item R Function Applications
#' \item Hill Function Applications: Computes multiple Hill functions for different regulatory components
#' \item Loss Computation: Calculates mean squared prediction errors across all time points
#' }
#' The loss function incorporates weighted squared errors for TFE, TGA, and TGE predictions, with TFE and TGA errors weighted twice as heavily as TGE errors.
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item The function depends on optimized parameters from a previous training process
#' \item The TG-TF identifier must follow the specific format for proper parsing
#' \item The data parameter must contain the required columns with proper naming conventions
#' \item The function is typically used for model evaluation on validation or test data
#' }
#'
#' @seealso \code{\link{R_cal}}, \code{\link{Hill_cal}}, \code{\link{S_cal}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # The settings for the other parameters are inside the function scDEDS::model_train.
#' print(loss_cal(
#'   data = train_data,
#'   params = result_gradient$parameters,
#'   params_names = params_names,
#'   tao = tao,
#'   Tslot_K = Tslot_K
#' ))
#' }
loss_cal = function (
    data = train_data,
    params = result_gradient$parameters,
    params_names = params_names,
    tao = tao,
    Tslot_K = Tslot_K
)
{
  base::names(params) = params_names
  TFE = base::as.numeric(data[, base::grep("TFE_", base::colnames(data), value = TRUE)])
  TGA = base::as.numeric(data[, base::grep("TGA_", base::colnames(data), value = TRUE)])
  TGE = base::as.numeric(data[, base::grep("TGE_", base::colnames(data), value = TRUE)])
  Length = base::length(TFE)
  seqLength = base::seq_len(Length - 1)
  R = scDEDS::R_cal(
    theta_TF_TG = base::as.numeric(data[, "theta_i"]),
    tao = tao,
    r1 = 200,
    r2 = (params["r.r2_TG"] + params["r.r2_TF"])/2,
    r3 = 1,
    r4 = 1000,
    r5 = 1
  )
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
  base::sum(base::sapply(2:Length, function(K) {
    prev_idx = K - 1
    2 * (TFE[prev_idx] + alpha1[prev_idx] * R * scDEDS::S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde[prev_idx]) + beta1[prev_idx] - TFE[K])^2 +
      2 * (TGA[prev_idx] + alpha2[prev_idx] * R * scDEDS::S_cal(s = s2, U = U2[prev_idx], U_tilde = U2_tilde[prev_idx]) + beta2[prev_idx] - TGA[K])^2 +
      (TGE[prev_idx] + R * (v1 * U31[prev_idx] - v2 * U32[prev_idx] - v3 * U33[prev_idx]) * Tslot_K[K] + beta3[prev_idx] - TGE[K])^2
  }))
}
