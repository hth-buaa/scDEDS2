#' @title The Hill-type equation/function in the evolutionary equations of the paper
#'
#' @param x Numeric (vector) greater than 0. Concentration or expression level of the regulatory molecule.
#' @param Dissociation_Constant
#' Numeric greater than 0. The dissociation constant (K_d) representingthe concentration at which half-maximal binding occurs.
#' Lower values indicate higher binding affinity.
#' @param Hill_Coefficient
#' Numeric greater than 0. The Hill coefficient (n) representing the degree of cooperativity in binding.
#' Values greater than 1 indicate positive cooperativity, values equal to 1 indicate non-cooperative binding, and values less than 1 indicate negative cooperativity.
#'
#' @returns
#' A numeric value(s) between 0 and 1 representing the fractional occupancy or binding probability of the regulatory molecule.
#' A value of 0 indicates no binding, 1 indicates complete saturation, and intermediate values represent partial binding.
#' @export
#'
#' @examples
#' \dontrun{
#' U1_tilde = Hill_cal(
#'   x = TGE_tilde,
#'   Dissociation_Constant = u11,
#'   Hill_Coefficient = u12
#' )
#' }
Hill_cal = function(x, Dissociation_Constant, Hill_Coefficient)
{
  x = as.numeric(x)
  Dissociation_Constant = as.numeric(Dissociation_Constant)
  Hill_Coefficient = as.numeric(Hill_Coefficient)
  U_scale = x ^ Hill_Coefficient
  U = U_scale / (Dissociation_Constant + U_scale)
  return(U)
}
