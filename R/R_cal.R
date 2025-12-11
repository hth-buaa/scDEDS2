#' @title The function R in the evolutionary equations of the paper
#'
#' @param theta_TF_TG Numeric (vector). Regulatory strength of transcription factor TF on target gene TG
#' @param tao
#' Numeric from 0 to 1. Regulatory threshold τ.
#' If regulatory strength θ > τ, a regulatory relationship is confirmed (positive sample); o
#' therwise, it is considered a negative sample.
#' @param r1 Numeric greater than 0. A hyperparameter. Generally greater than 100. Default is 200.
#' @param r2 Numeric greater than 0. A hyperparameter.
#' @param r3 Numeric greater than 0. A hyperparameter. Default is 1.
#' @param r4 Numeric greater than 0. A hyperparameter. Generally greater than 100. Default is 1000.
#' @param r5 Numeric greater than 0. A hyperparameter. Default is 1.
#'
#' @returns Numeric (vector) from 0 to 1. The computation results of function R.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' R = R_cal(
#'   theta_TF_TG = c(0.18, 0.23),
#'   tao = 0.323,
#'   r1 = 200,
#'   r2 = 1,
#'   r3 = 1,
#'   r4 = 1000,
#'   r5 = 1
#' )
#' }
R_cal = function(theta_TF_TG, tao, r1 = 200, r2, r3 = 1, r4 = 1000, r5 = 1)
{
  R_unnormalize = theta_TF_TG / (r1 * tao) + r5 / (
    (1 + 10^(r2 * (r3 - theta_TF_TG))) * (1 + 10^(r4 * r2 * (r3 - theta_TF_TG)))
  )
  r6 = 1 / (r1 * tao) + r5 / (
    (1 + 10^(r2 * (r3 - 1))) * (1 + 10^(r4 * r2 * (r3 - 1)))
  )
  R = R_unnormalize / r6
  return(R)
}
