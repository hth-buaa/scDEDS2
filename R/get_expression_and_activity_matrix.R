#' @title Obtain Expression and Gene Activity Matrices for GRN Construction
#'
#' @description
#' This function processes single-cell multi-omics data to extract expression matrices for target genes (TGs) and transcription factors (TFs),
#' and computes gene activity matrices from scATAC-seq data for GRN (Gene Regulatory Network) construction.
#'
#' @param TGs
#' Character vector of target gene symbols to include in the analysis.
#' Default is results_identify_TGs$TGs. results_identify_TGs is the output of function identify_TGs.
#' @param TFs
#' Character vector of transcription factor symbols to include in the analysis.
#' Default is TFs_in_JASPAR2024. TFs_in_JASPAR2024 is the output of function get_TGs_from_JASPAR2024.
#' @param scRNAseq See in ?identify_TGs.
#' @param scATACseq See in ?identify_TGs.
#' @param promoter_range
#' Numeric. The genomic range (in base pairs) to extend around transcription start sites for calculating gene activity scores.
#' Must be set prior to code execution and remain unchanged during runtime.
#' @param path
#' Absolute path to the fragment tbi. file (both .tsv.gz and .tbi index must exist in the same directory).
#' For example, "/mnt/d/pbmc-equation/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz".
#'
#' @returns
#' A list containing five elements:
#' \itemize{
#' \item \code{scRNAseq_for_GRN}: Subsetted scRNA-seq data containing only TGs and TFs
#' \item \code{scATACseq_for_GRN}: Subsetted scATAC-seq data with gene activity matrix (only TGs and TFs)
#' \item \code{gene_for_GRN}: charactor vector of TGs and TFs
#' \item \code{TGs}: Filtered target genes present in both datasets
#' \item \code{TFs}: Filtered transcription factors present in both datasets
#' }
#'
#' @details
#' The function performs the following key operations:
#' \enumerate{
#' \item Creates a dedicated output directory for organization
#' \item Extracts expression matrices for TGs and TFs from scRNA-seq data
#' \item Annotates chromatin fragments and computes gene activity scores from scATAC-seq data
#' \item Integrates both datasets and filters genes to ensure consistency
#' \item Returns processed data ready for GRN inference
#' }
#'
#' @note
#' Important requirements:
#' \itemize{
#' \item Requires both fragment file and its index (.tbi) in the specified path
#' \item Uses EnsDb.Hsapiens.v86 for genomic annotations (modify for other species)
#' \item Temporarily changes working directory but restores upon completion
#' \item activity is calculated based on ATAC-seq fragment counts in gene bodies and promoter regions
#' \item This function must run in Linux, as \code{GenomeInfoDb::seqlevelsStyle(annotations) = 'UCSC'} is guaranteed to work in Linux but has very low success rate on Windows.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run in linux.
#' library(EnsDb.Hsapiens.v86)
#' results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
#' TFs_in_JASPAR2024 = base::readRDS("./1.2 Data Preprocessing - Get TGs From JASPAR2024/TFs_in_JASPAR2024.rds")
#' scRNAseq = base::readRDS("./scRNAseq.rds")
#' scATACseq = base::readRDS("./scATACseq.rds")
#' promoter_range = 50
#' Basic_Info = get_expression_and_activity_matrix(
#'   TGs = results_identify_TGs$TGs,
#'   TFs = TFs_in_JASPAR2024,
#'   scRNAseq = scRNAseq,
#'   scATACseq = scATACseq,
#'   promoter_range = promoter_range,
#'   path = "/mnt/d/pbmc-equation/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
#' )
#' base::saveRDS(Basic_Info, file = "./1.3 Data Preprocessing - Get Expression And Activity Matrix/Basic_Info.rds")
#' }
get_expression_and_activity_matrix = function(
    TGs = results_identify_TGs$TGs,
    TFs = TFs_in_JASPAR2024,
    scRNAseq = scRNAseq,
    scATACseq = scATACseq,
    promoter_range = 50,
    path
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Obtaining Counts Matrices for TGs/TFs Expression and TGs Activity ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "1.3 Data Preprocessing - Get Expression And Activity Matrix"
  if (!base::dir.exists(new_folder)) {base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Extracting gene expression matrix.
  gene_for_GRN = base::union(TGs, TFs)
  message("Extracting gene expression matrix.")
  scRNAseq_for_GRN = base::subset(scRNAseq, features = gene_for_GRN)

  ### Generating gene activity matrix (based on ATAC-seq fragment counts in gene bodies and promoter regions per cell)
  message("Annotating chromatin fragments in scATAC-seq data with genomic annotation.")
  annotations = Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
  GenomeInfoDb::seqlevelsStyle(annotations) = 'UCSC'
  Signac::Annotation(scATACseq) = annotations

  message("Adding fragments. Must have both _fragments.tsv.gz and _fragments.tsv.gz.tbi in the current directory.")
  fragment = Signac::CreateFragmentObject(path = path)
  fragment@cells = SeuratObject::Cells(scATACseq)
  Signac::Fragments(scATACseq) = NULL
  Signac::Fragments(scATACseq) = fragment

  message("Adding an 'ACTIVITY' matrix to scATAC-seq data in Seurat object.")
  gene.activities = Signac::GeneActivity(scATACseq,
                                         features = base::rownames(scRNAseq),
                                         extend.upstream = promoter_range,
                                         extend.downstream = promoter_range)
  scATACseq[["ACTIVITY"]] = SeuratObject::CreateAssayObject(counts = gene.activities)
  SeuratObject::DefaultAssay(scATACseq) = "ACTIVITY"
  scATACseq_for_GRN = base::subset(scATACseq, features = gene_for_GRN)

  ### Update data.
  message("Update data.")
  gene_for_GRN = base::intersect(base::rownames(scRNAseq_for_GRN), base::rownames(scATACseq_for_GRN))
  TGs = base::intersect(TGs, gene_for_GRN)
  TFs = base::intersect(TFs, gene_for_GRN)
  scRNAseq_for_GRN = base::subset(scRNAseq, features = gene_for_GRN)
  scATACseq_for_GRN = base::subset(scATACseq, features = gene_for_GRN)
  SeuratObject::DefaultAssay(scATACseq_for_GRN) = "ACTIVITY"

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Obtaining Counts Matrices for TGs/TFs Expression and TGs Activity ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(base::list("scRNAseq_for_GRN" = scRNAseq_for_GRN,
                    "scATACseq_for_GRN" = scATACseq_for_GRN,
                    "gene_for_GRN" = gene_for_GRN,
                    "TGs" = TGs,
                    "TFs" = TFs))
}
