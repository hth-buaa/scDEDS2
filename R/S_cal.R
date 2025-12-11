#' @title The function S1 or S2 in the evolutionary equations of the paper
#'
#' @param s Numeric greater than 0. cGenerally range from 0.2 to 5.
#' @param U Numeric ranging from 0 to 1. An intermediate variable.
#' @param U_tilde Numeric ranging from 0 to 1. Numeric ranging from 0 to 1.
#'
#' @returns Numeric. The computation results of function  S1 or S2.
#' @export
#'
#' @examples
#' \dontrun{
#' print(S_cal(s = s1, U = U1[prev_idx], U_tilde = U1_tilde[prev_idx]))
#' }
S_cal = function(s, U, U_tilde)
{
  Um = 1 - U
  Um_tilde = 1 - U_tilde
  Numerator = s * (U * Um)^(s-1) * (U_tilde - U) * ((1 - U) * (1 - U_tilde) + U_tilde * U)
  Denominator = (U_tilde * Um_tilde)^(s+1) * exp(s * (2 * U_tilde - 1) * (U_tilde - U) / (U_tilde * Um_tilde))
  S = Numerator / Denominator
  return(S)
}
