#' @title Initialize Model Parameters for Gene Regulatory Network Prediction
#'
#' @description
#' This function initializes the parameter structures for the gene regulatory network prediction model scDEDS.
#'
#' @param interest_cell_type_branch_training_set The output of function select_training_set.
#'
#' @returns A nested list structure containing initialized parameters for each cell type and branch.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#' \item Creates a hierarchical parameter structure for each cell type, branch, and regulatory interaction (TF-TG pair).
#' \item Initializes parameters for the scDEDS model, including:
#' \itemize{
#' \item \code{r}: Regulatory interaction parameters (r2_TG, r2_TF) with initial values, upper, and lower bounds.
#' \item \code{alpha}: Cell group state association parameters (alpha1, alpha2 for TG and TF) for each time point K.
#' \item \code{beta}: Linear coefficients (beta1, beta2, beta3 for TG and TF) for each time point K.
#' \item \code{s}: Scaling factors (s1_TG, s1_TF, s2_TG, s2_TF).
#' \item \code{u}: Non-linear transformation parameters (u11_TG, u21_TG, u311_TG, u321_TG, u331_TG, u12_TG, u22_TG, u312_TG, u322_TG, u332_TG).
#' \item \code{E~}: Baseline expression levels (E~TG_K-1, E~TF_K-1) for each time point K, initialized as the mean of TGE and TFE from the training set, respectively.
#' \item \code{v}: Additional scaling parameters (v1_TG, v2_TG, v3_TG).
#' }
#' \item For each TF-TG pair, the function uses the mean TFE and TGE values from the training set to initialize the E~ parameters.
#' \item All parameters are assigned initial values, upper bounds, and lower bounds to constrain the optimization process.
#' \item The function organizes the parameters in a nested list structure, first by cell type, then by branch, and finally by regulatory interaction (TF-TG pair).
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_branch_training_set = base::readRDS("./4.1 Build Prediction Model - Select Training Set/interest_cell_type_branch_training_set.rds")
#' interest_cell_type_branch_init_params = set_init_params(
#'   interest_cell_type_branch_training_set = interest_cell_type_branch_training_set
#' )
#' base::saveRDS(interest_cell_type_branch_init_params, file = "./4.2 BUild Prediction Model/interest_cell_type_branch_init_params.rds")
#' }
set_init_params = function (interest_cell_type_branch_training_set = interest_cell_type_branch_training_set) {
  t_start = base::Sys.time()
  message("Run: Setting initial parameter values for all branches ", t_start, ".")
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

  set_init_params_process = function(Ks, TFE_mean, TGE_mean) {
    message("Initializing with a list.")
    init_params = base::list()

    message("Setting r.")
    init_params[["r"]] = base::list()

    message("Setting r2 in r.")
    init_params[["r"]][["r2_TG"]] = base::list()
    init_params[["r"]][["r2_TG"]][["initialization"]] = 1
    init_params[["r"]][["r2_TG"]][["upper"]] = 1.113383
    init_params[["r"]][["r2_TG"]][["lower"]] = 0.370506

    init_params[["r"]][["r2_TF"]] = base::list()
    init_params[["r"]][["r2_TF"]][["initialization"]] = 1
    init_params[["r"]][["r2_TF"]][["upper"]] = 1.113383
    init_params[["r"]][["r2_TF"]][["lower"]] = 0.370506

    message("Setting alpha.")
    init_params[["alpha"]] = base::list()

    message("Setting alpha1 in alpha.")
    init_params[["alpha"]][["alpha1_TG_K-1"]] = base::list()
    init_params[["alpha"]][["alpha1_TG_K-1"]][["initialization"]] = base::rep(1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha1_TG_K-1"]][["initialization"]]) = Ks - 1
    init_params[["alpha"]][["alpha1_TG_K-1"]][["upper"]] = base::rep(50, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha1_TG_K-1"]][["upper"]]) = Ks - 1
    init_params[["alpha"]][["alpha1_TG_K-1"]][["lower"]] = base::rep(0.1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha1_TG_K-1"]][["lower"]]) = Ks - 1

    init_params[["alpha"]][["alpha1_TF_K-1"]] = base::list()
    init_params[["alpha"]][["alpha1_TF_K-1"]][["initialization"]] = base::rep(1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha1_TF_K-1"]][["initialization"]]) = Ks - 1
    init_params[["alpha"]][["alpha1_TF_K-1"]][["upper"]] = base::rep(50, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha1_TF_K-1"]][["upper"]]) = Ks - 1
    init_params[["alpha"]][["alpha1_TF_K-1"]][["lower"]] = base::rep(0.1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha1_TF_K-1"]][["lower"]]) = Ks - 1

    message("Setting alpha2 in alpha.")
    init_params[["alpha"]][["alpha2_TG_K-1"]] = base::list()
    init_params[["alpha"]][["alpha2_TG_K-1"]][["initialization"]] = base::rep(1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha2_TG_K-1"]][["initialization"]]) = Ks - 1
    init_params[["alpha"]][["alpha2_TG_K-1"]][["upper"]] = base::rep(50, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha2_TG_K-1"]][["upper"]]) = Ks - 1
    init_params[["alpha"]][["alpha2_TG_K-1"]][["lower"]] = base::rep(0.1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha2_TG_K-1"]][["lower"]]) = Ks - 1

    init_params[["alpha"]][["alpha2_TF_K-1"]] = base::list()
    init_params[["alpha"]][["alpha2_TF_K-1"]][["initialization"]] = base::rep(1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha2_TF_K-1"]][["initialization"]]) = Ks - 1
    init_params[["alpha"]][["alpha2_TF_K-1"]][["upper"]] = base::rep(50, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha2_TF_K-1"]][["upper"]]) = Ks - 1
    init_params[["alpha"]][["alpha2_TF_K-1"]][["lower"]] = base::rep(0.1, base::length(Ks))
    base::names(init_params[["alpha"]][["alpha2_TF_K-1"]][["lower"]]) = Ks - 1

    message("Setting beta.")
    init_params[["beta"]] = base::list()

    message("Setting beta1 in beta.")
    init_params[["beta"]][["beta1_TG_K-1"]] = base::list()
    init_params[["beta"]][["beta1_TG_K-1"]][["initialization"]] = base::rep(0, base::length(Ks))
    base::names(init_params[["beta"]][["beta1_TG_K-1"]][["initialization"]]) = Ks - 1
    init_params[["beta"]][["beta1_TG_K-1"]][["upper"]] = base::rep(10, base::length(Ks))
    base::names(init_params[["beta"]][["beta1_TG_K-1"]][["upper"]]) = Ks - 1
    init_params[["beta"]][["beta1_TG_K-1"]][["lower"]] = base::rep(-10, base::length(Ks))
    base::names(init_params[["beta"]][["beta1_TG_K-1"]][["lower"]]) = Ks - 1

    init_params[["beta"]][["beta1_TF_K-1"]] = base::list()
    init_params[["beta"]][["beta1_TF_K-1"]][["initialization"]] = base::rep(0, base::length(Ks))
    base::names(init_params[["beta"]][["beta1_TF_K-1"]][["initialization"]]) = Ks - 1
    init_params[["beta"]][["beta1_TF_K-1"]][["upper"]] = base::rep(10, base::length(Ks))
    base::names(init_params[["beta"]][["beta1_TF_K-1"]][["upper"]]) = Ks - 1
    init_params[["beta"]][["beta1_TF_K-1"]][["lower"]] = base::rep(-10, base::length(Ks))
    base::names(init_params[["beta"]][["beta1_TF_K-1"]][["lower"]]) = Ks - 1

    message("Setting beta2 in beta.")
    init_params[["beta"]][["beta2_TG_K-1"]] = base::list()
    init_params[["beta"]][["beta2_TG_K-1"]][["initialization"]] = base::rep(0, base::length(Ks))
    base::names(init_params[["beta"]][["beta2_TG_K-1"]][["initialization"]]) = Ks - 1
    init_params[["beta"]][["beta2_TG_K-1"]][["upper"]] = base::rep(10, base::length(Ks))
    base::names(init_params[["beta"]][["beta2_TG_K-1"]][["upper"]]) = Ks - 1
    init_params[["beta"]][["beta2_TG_K-1"]][["lower"]] = base::rep(-10, base::length(Ks))
    base::names(init_params[["beta"]][["beta2_TG_K-1"]][["lower"]]) = Ks - 1

    init_params[["beta"]][["beta2_TF_K-1"]] = base::list()
    init_params[["beta"]][["beta2_TF_K-1"]][["initialization"]] = base::rep(0, base::length(Ks))
    base::names(init_params[["beta"]][["beta2_TF_K-1"]][["initialization"]]) = Ks - 1
    init_params[["beta"]][["beta2_TF_K-1"]][["upper"]] = base::rep(10, base::length(Ks))
    base::names(init_params[["beta"]][["beta2_TF_K-1"]][["upper"]]) = Ks - 1
    init_params[["beta"]][["beta2_TF_K-1"]][["lower"]] = base::rep(-10, base::length(Ks))
    base::names(init_params[["beta"]][["beta2_TF_K-1"]][["lower"]]) = Ks - 1

    message("Setting beta3 in beta.")
    init_params[["beta"]][["beta3_TG_K-1"]] = base::list()
    init_params[["beta"]][["beta3_TG_K-1"]][["initialization"]] = base::rep(0, base::length(Ks))
    base::names(init_params[["beta"]][["beta3_TG_K-1"]][["initialization"]]) = Ks - 1
    init_params[["beta"]][["beta3_TG_K-1"]][["upper"]] = base::rep(10, base::length(Ks))
    base::names(init_params[["beta"]][["beta3_TG_K-1"]][["upper"]]) = Ks - 1
    init_params[["beta"]][["beta3_TG_K-1"]][["lower"]] = base::rep(-10, base::length(Ks))
    base::names(init_params[["beta"]][["beta3_TG_K-1"]][["lower"]]) = Ks - 1

    message("Setting s.")
    init_params[["s"]] = base::list()

    message("Setting s1 in s.")
    init_params[["s"]][["s1_TG"]] = base::list()
    init_params[["s"]][["s1_TG"]][["initialization"]] = 1
    init_params[["s"]][["s1_TG"]][["upper"]] = 4.9
    init_params[["s"]][["s1_TG"]][["lower"]] = 0.1

    init_params[["s"]][["s1_TF"]] = base::list()
    init_params[["s"]][["s1_TF"]][["initialization"]] = 1
    init_params[["s"]][["s1_TF"]][["upper"]] = 4.9
    init_params[["s"]][["s1_TF"]][["lower"]] = 0.1

    message("Setting s2 in s.")
    init_params[["s"]][["s2_TG"]] = base::list()
    init_params[["s"]][["s2_TG"]][["initialization"]] = 1
    init_params[["s"]][["s2_TG"]][["upper"]] = 4.9
    init_params[["s"]][["s2_TG"]][["lower"]] = 0.1

    init_params[["s"]][["s2_TF"]] = base::list()
    init_params[["s"]][["s2_TF"]][["initialization"]] = 1
    init_params[["s"]][["s2_TF"]][["upper"]] = 4.9
    init_params[["s"]][["s2_TF"]][["lower"]] = 0.1

    message("Setting u.")
    init_params[["u"]] = base::list()

    message("Setting u1 in u.")
    init_params[["u"]][["u11_TG"]] = base::list()
    init_params[["u"]][["u11_TG"]][["initialization"]] = 1
    init_params[["u"]][["u11_TG"]][["upper"]] = 777600
    init_params[["u"]][["u11_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u21_TG"]] = base::list()
    init_params[["u"]][["u21_TG"]][["initialization"]] = 1
    init_params[["u"]][["u21_TG"]][["upper"]] = 777600
    init_params[["u"]][["u21_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u311_TG"]] = base::list()
    init_params[["u"]][["u311_TG"]][["initialization"]] = 1
    init_params[["u"]][["u311_TG"]][["upper"]] = 777600
    init_params[["u"]][["u311_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u321_TG"]] = base::list()
    init_params[["u"]][["u321_TG"]][["initialization"]] = 1
    init_params[["u"]][["u321_TG"]][["upper"]] = 777600
    init_params[["u"]][["u321_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u331_TG"]] = base::list()
    init_params[["u"]][["u331_TG"]][["initialization"]] = 1
    init_params[["u"]][["u331_TG"]][["upper"]] = 777600
    init_params[["u"]][["u331_TG"]][["lower"]] = 0.1

    message("Setting u2 in u.")
    init_params[["u"]][["u12_TG"]] = base::list()
    init_params[["u"]][["u12_TG"]][["initialization"]] = 1
    init_params[["u"]][["u12_TG"]][["upper"]] = 4.9
    init_params[["u"]][["u12_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u22_TG"]] = base::list()
    init_params[["u"]][["u22_TG"]][["initialization"]] = 1
    init_params[["u"]][["u22_TG"]][["upper"]] = 4.9
    init_params[["u"]][["u22_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u312_TG"]] = base::list()
    init_params[["u"]][["u312_TG"]][["initialization"]] = 1
    init_params[["u"]][["u312_TG"]][["upper"]] = 4.9
    init_params[["u"]][["u312_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u322_TG"]] = base::list()
    init_params[["u"]][["u322_TG"]][["initialization"]] = 1
    init_params[["u"]][["u322_TG"]][["upper"]] = 4.9
    init_params[["u"]][["u322_TG"]][["lower"]] = 0.1

    init_params[["u"]][["u332_TG"]] = base::list()
    init_params[["u"]][["u332_TG"]][["initialization"]] = 1
    init_params[["u"]][["u332_TG"]][["upper"]] = 4.9
    init_params[["u"]][["u332_TG"]][["lower"]] = 0.1

    message("Setting E~.")
    init_params[["E~"]] = base::list()

    message("Setting E~_TG in E~.")
    init_params[["E~"]][["E~_TG_K-1"]] = base::list()
    init_params[["E~"]][["E~_TG_K-1"]][["initialization"]] = base::rep(TGE_mean, base::length(Ks))
    base::names(init_params[["E~"]][["E~_TG_K-1"]][["initialization"]]) = Ks - 1
    init_params[["E~"]][["E~_TG_K-1"]][["upper"]] = base::rep(9, base::length(Ks))
    base::names(init_params[["E~"]][["E~_TG_K-1"]][["upper"]]) = Ks - 1
    init_params[["E~"]][["E~_TG_K-1"]][["lower"]] = base::rep(1, base::length(Ks))
    base::names(init_params[["E~"]][["E~_TG_K-1"]][["lower"]]) = Ks - 1

    message("Setting E~_TF in E~.")
    init_params[["E~"]][["E~_TF_K-1"]] = base::list()
    init_params[["E~"]][["E~_TF_K-1"]][["initialization"]] = base::rep(TFE_mean, base::length(Ks))
    base::names(init_params[["E~"]][["E~_TF_K-1"]][["initialization"]]) = Ks - 1
    init_params[["E~"]][["E~_TF_K-1"]][["upper"]] = base::rep(9, base::length(Ks))
    base::names(init_params[["E~"]][["E~_TF_K-1"]][["upper"]]) = Ks - 1
    init_params[["E~"]][["E~_TF_K-1"]][["lower"]] = base::rep(1, base::length(Ks))
    base::names(init_params[["E~"]][["E~_TF_K-1"]][["lower"]]) = Ks - 1

    message("Setting v.")
    init_params[["v"]] = base::list()

    message("Setting v1 in v.")
    init_params[["v"]][["v1_TG"]] = base::list()
    init_params[["v"]][["v1_TG"]][["initialization"]] = 1
    init_params[["v"]][["v1_TG"]][["upper"]] = 200
    init_params[["v"]][["v1_TG"]][["lower"]] = 0.1

    message("Setting v2 in v.")
    init_params[["v"]][["v2_TG"]] = base::list()
    init_params[["v"]][["v2_TG"]][["initialization"]] = 1
    init_params[["v"]][["v2_TG"]][["upper"]] = 200
    init_params[["v"]][["v2_TG"]][["lower"]] = 0.1

    message("Setting v3 in v.")
    init_params[["v"]][["v3_TG"]] = base::list()
    init_params[["v"]][["v3_TG"]][["initialization"]] = 1
    init_params[["v"]][["v3_TG"]][["upper"]] = 200
    init_params[["v"]][["v3_TG"]][["lower"]] = 0.1

    return(init_params)
  }
  interest_cell_type_branch_init_params = base::list()
  for (cell_type in base::names(interest_cell_type_branch_training_set)) {
    message("Setting initial parameter values for the data of cell type ", cell_type, ".")
    interest_cell_type_branch_init_params[[cell_type]] = base::list()
    for (n in 1:(base::length(interest_cell_type_branch_training_set[[cell_type]]) - 1)) {
      interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]] = base::list()
      for (re in rownames(interest_cell_type_branch_training_set[[cell_type]][[base::paste0("branch", n)]])) {
        interest_cell_type_branch_init_params[[cell_type]][[base::paste0("branch", n)]][[re]] = set_init_params_process(
          Ks = 2:base::sum(base::grepl("TFE_", base::colnames(interest_cell_type_branch_training_set[[cell_type]][[n]]))),
          TFE_mean = base::mean(base::as.numeric(dplyr::select(interest_cell_type_branch_training_set[[cell_type]][[n]], dplyr::matches("TFE"))[re, ])),
          TGE_mean = base::mean(base::as.numeric(dplyr::select(interest_cell_type_branch_training_set[[cell_type]][[n]], dplyr::matches("TGE"))[re, ]))
        )
      }
    }
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Setting initial parameter values for all branches ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_branch_init_params)
}
