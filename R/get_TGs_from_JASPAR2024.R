#' @title Retrieve Transcription Factors Genes (TFs, TFGs) from JASPAR2024 Database
#'
#' @description
#' This function extracts transcription factors Genes (TFs, TFGs) from the JASPAR2024 database that are present in the provided single-cell RNA sequencing data.
#' It filters TFs based on species and collection parameters, and returns a vector of TF names available in both JASPAR2024 and the input data.
#'
#' @param scRNAseq See in ?identify_TGs.
#' @param species
#' Character. Species to query in JASPAR2024 database. Default is "Homo sapiens".
#' Parameter species and genome must be compatible.
#' Execute this command to verify: \code{rbioapi::rba_jaspar_species(2024)$species}.
#' @param collection
#' Collection Character vector. Specifies which JASPAR sub-collections to search through.
#' Available options include: "CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED".
#' Default includes all.
#'
#' @returns A character vector of transcription factor names that exist in both JASPAR2024 database and the provided scRNAseq data.
#'
#' #' @details
#' The function performs the following operations:
#' \enumerate{
#' \item Connects to JASPAR2024 database and retrieves position weight matrices (PWMs)
#' \item Filters TFs to retain only those present in the scRNA-seq data
#' \item Creates a dedicated directory for output organization
#' \item Returns a list of validated TF names
#' }
#'
#' @note
#' This function temporarily changes the working directory to save results but reverts to the original directory upon completion.
#' The function assumes that all TFs have potential regulatory relationships with each target gene.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Retrieve human TFs from core JASPAR collections
#' scRNAseq = base::readRDS("./scRNAseq.rds")
#' TFs_in_JASPAR2024 = get_TGs_from_JASPAR2024(
#'   scRNAseq = scRNAseq,
#'   species = "Homo sapiens",
#'   collection = "CORE"
#' )
#' base::saveRDS(TFs_in_JASPAR2024, file = "./1.2 Data Preprocessing - Get TGs From JASPAR2024/TFs_in_JASPAR2024.rds")
#' }
get_TGs_from_JASPAR2024 = function(
    scRNAseq = scRNAseq,
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                   "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED")
  )
{
  ### Start.
  t_start = base::Sys.time()
  original_dir = base::getwd()
  new_folder = "1.2 Data Preprocessing - Get TGs From JASPAR2024"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  message("Run: Retrieve TFGs from JASPAR2024 Database ", t_start, ".")

  ### Get TFDB PWMs from JASPAR.
  message("Retrieve known TFs from JASPAR2024, assuming all TFs have potential regulatory relationships with each TG.")
  pwm = TFBSTools::getMatrixSet(
    x = JASPAR2024::db(JASPAR2024::JASPAR2024()),
    opts = base::list(all = FALSE, species = species, collection = collection, all_versions = FALSE, matrixtype = "PWM")
  )

  ### Filter TFs in scRNA-seq data.
  message("Filter TFs in scRNA-seq data.")
  pwm@listData = base::Filter(function(x) x@name %in% base::rownames(scRNAseq), pwm@listData)
  TFs = base::sapply(pwm@listData, function(x) x@name)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Retrieve TFGs from JASPAR2024 Database ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(TFs)
}
