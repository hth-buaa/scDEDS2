#' @title Infer Peak-TG Relationships
#'
#' @description
#' Infers regulatory relationships between chromatin accessibility peaks and target genes (TGs) by integrating scRNA-seq and scATAC-seq data with pre-computed gene regulatory network (GRN) predictions.
#' The function processes each cell type independently at two stringency levels (loose and strict) to identify target genes present in the filtered GRN, then retrieves promoter-annotated peaks associated with these genes and calculates multiple importance scores for each peak-TG pair.
#'
#' @param scRNAseq See in ?identify_TGs.
#' @param scATACseq See in ?identify_TGs.
#' @param interest_cell_type_pGRN The output of function filter_GRN.
#' @param results_identify_TGs The output of function identify_TGs.
#' @param weight_Accessibility_percentile A positive numeric value specifying the weight for the accessibility percentile in importance score calculations. Default is 1.
#' @param weight_Expression_percentile A positive numeric value specifying the weight for the expression percentile in importance score calculations. Default is 1.
#' @param weight_Accessibility_non_zero_ratio_percentile A positivenumeric value specifying the weight for the non-zero ratio percentile of accessibility in importance score calculations. Default is 5.
#' @param weight_Expression_non_zero_ratio_percentile A positivenumeric value specifying the weight for the non-zero ratio percentile of expression in importance score calculations. Default is 5.
#'
#' @returns
#' A nested list named \code{interest_cell_type_peak_TG_pred}, organized by cell type and filtering level (loose/strict). For each cell type and level, the list contains:
#' \itemize{
#' \item \code{TGs}: A character vector of unique target gene symbols that are present in both the filtered GRN predictions and the promoter peak annotations.
#' \itemize{
#' \item \code{peak}: Mean accessibility of the peak in non-zero cells.
#' \item \code{SYMBOL}: Mean accessibility of the peak in non-zero cells.
#' \item \code{Accessibility}: Mean accessibility of the peak in non-zero cells.
#' \item \code{Expression}: Mean expression of the TG in non-zero cells.
#' \item \code{Accessibility_percentile}, \code{Expression_percentile}: Percentile ranks of accessibility and expression.
#' \item \code{Accessibility_non_zero_ratio}, \code{Expression_non_zero_ratio}: Proportions of non-zero cells for the peak and TG.
#' \item \code{Accessibility_non_zero_ratio_percentile}, \code{Expression_non_zero_ratio_percentile}: Percentile ranks of Accessibility_non_zero_ratio and Expression_non_zero_ratio.
#' \item \code{Accessibility_non_zero_ratio}, \code{Expression_non_zero_ratio}: Proportions of non-zero cells for the peak and TG.
#' \item \code{peak_importance_score}, \code{TG_importance_score}, \code{peak_TG_importance_score}: Composite scores combining weighted percentiles.
#' \item \code{peak_importance_score_percentile}, \code{TG_importance_score_percentile}, \code{peak_TG_importance_score_percentile}: Percentile ranks of peak_importance_score, TG_importance_score, and peak_TG_importance_score.
#' \item \code{TG_count}, \code{peak_count}: The number of TG corresponding to each peak and the number of peak corresponding to each TG.
#' \item \code{TG_count_percentile}, \code{peak_count_percentile}: Percentile ranks of TG_count and peak_count
#' }
#' }
#' The structure is: \code{interest_cell_type_peak_TG_pred[[cell_type]][[level]][["TGs"]]} and \code{interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]]}.
#'
#' @details
#' This function performs the following key steps for each cell type and stringency level:
#' \enumerate{
#' \item \strong{Extract Target Genes}: Retrieves unique TG symbols from the GRN predictions (rows are named as "TF_to_TG").
#' \item \strong{Match Promoter Peaks}: Subsets the promoter-annotated peaks from \code{results_identify_TGs} to include only peaks associated with the extracted TGs. It merges gene symbols from the \code{SYMBOL} and \code{flank_geneSymbols} columns, splits entries by semicolon, and expands the data frame so each peak-TG pair has a single row.
#' \item \strong{Calculate Accessibility and Expression Metrics}: Computes mean accessibility and expression for non-zero cells, along with non-zero ratios (proportion of cells with non-zero values). Percentile ranks are calculated for these metrics.
#' \item \strong{Compute Importance Scores}: Generates composite scores using user-defined weights to prioritize peak-TG pairs based on accessibility and expression characteristics.
#' \item \strong{Filter and Structure Results}: Removes pairs where accessibility or expression is zero, and ensures TGs are consistent with the GRN predictions.
#' }
#' The "loose" level typically includes a broader set of TGs, while "strict" provides high-confidence targets. This enables integrative analysis of chromatin accessibility and GRN predictions to identify potential regulatory elements.
#'
#' @note
#' Important considerations when using this function:
#' \itemize{
#' \item The input \code{results_identify_TGs} must contain a properly formatted \code{peak_anno_promoter} component. Ensure promoter annotations are complete and use consistent gene symbols.
#' \item The function requires both scRNA-seq and scATAC-seq data to be normalized and processed (e.g., using \code{Seurat::NormalizeData} for RNA and \code{Signac::RunTFIDF} for ATAC). The data should cover the same set of cells or be integrated appropriately.
#' \item The weights for percentiles (e.g., \code{weight_Accessibility_percentile}) allow customization of importance scores. Adjust these based on the relative importance of accessibility versus expression in your analysis.
#' }
#' Afterwards, users can filter predicted TGs based on TG_importance_store or TG_importance_score_percentile, filter predicted peaks based on peak_importance_store or peak_importance_score_percentile, and filter predicted peak TG associations based on peak_TG_importance_store or peak_TG_importance_score_percentile.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' scRNAseq = readRDS("./scRNAseq.rds")
#' scATACseq = readRDS("./scATACseq.rds")
#' interest_cell_type_pGRN = base::readRDS("./5 Infer GRN/interest_cell_type_pGRN.rds")
#' results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
#' weight_Accessibility_percentile = 1
#' weight_Expression_percentile = 1
#' weight_Accessibility_non_zero_ratio_percentile = 5
#' weight_Expression_non_zero_ratio_percentile = 5
#' interest_cell_type_peak_TG_pred = infer_peak_to_TG (
#'   scRNAseq = scRNAseq,
#'   scATACseq = scATACseq,
#'   interest_cell_type_pGRN = interest_cell_type_pGRN,
#'   results_identify_TGs = results_identify_TGs,
#'   weight_Accessibility_percentile = weight_Accessibility_percentile,
#'   weight_Expression_percentile = weight_Expression_percentile,
#'   weight_Accessibility_non_zero_ratio_percentile = weight_Accessibility_non_zero_ratio_percentile,
#'   weight_Expression_non_zero_ratio_percentile = weight_Expression_non_zero_ratio_percentile
#' )
#' base::saveRDS(interest_cell_type_peak_TG_pred, file = "./5 Infer GRN/interest_cell_type_peak_TG_pred.rds")
#' }
infer_peak_to_TG = function(
    scRNAseq,
    scATACseq,
    interest_cell_type_pGRN = interest_cell_type_pGRN,
    results_identify_TGs = results_identify_TGs,
    weight_Accessibility_percentile = 1,
    weight_Expression_percentile = 1,
    weight_Accessibility_non_zero_ratio_percentile = 5,
    weight_Expression_non_zero_ratio_percentile = 5
)
{
  t_start = base::Sys.time()
  message("Run: Inferring Peak-TG Relationships ", t_start, ".")
  original_dir = base::getwd()
  new_folder = "5 Infer GRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {
    message("Folder already exists: ", new_folder, ".")
  }
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  interest_cell_type_peak_TG_pred = base::list()
  for (cell_type in base::names(interest_cell_type_pGRN)) {
    message("Inferring peak-TG relationships of cell type ", cell_type, ".")
    interest_cell_type_peak_TG_pred[[cell_type]] = base::list()

    for (level in c("loose", "strict")) {
      message("Inferring peak-TG relationships of cell type ", cell_type, " in ", level, " level.")
      interest_cell_type_peak_TG_pred[[cell_type]][[level]] = base::list()

      message("Getting TGs.")
      TGs = base::unique(base::sapply(base::strsplit(base::rownames(
        interest_cell_type_pGRN[[cell_type]][[base::paste0("pred_summary_", level)]]
      ), "_to_"), function(x) x[2]))

      message("Getting peak-TG relationships.")
      peak_TG = base::subset(
        results_identify_TGs[["peak_anno_promoter"]][, c("SYMBOL", "flank_geneSymbols"), drop = FALSE],
        SYMBOL %in% TGs
      ) %>%
        tibble::rownames_to_column("rowname") %>%
        dplyr::mutate(
          flank_geneSymbols_clean = stringr::str_replace_all(flank_geneSymbols, "\\s+", ""),
          flank_geneSymbols_clean = base::ifelse(flank_geneSymbols_clean == "", NA, flank_geneSymbols_clean)
        ) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(merged_genes = {
          symbol_gene = SYMBOL
          if (base::is.na(flank_geneSymbols_clean)) {gene_set = symbol_gene} else {
            flank_genes = stringr::str_split(flank_geneSymbols_clean, ";")[[1]]
            gene_set = base::unique(c(symbol_gene, flank_genes))
          }
          gene_set = gene_set[!base::is.na(gene_set) & gene_set != ""]
          if (base::length(gene_set) == 0) {return(NA_character_)} else {base::paste(gene_set, collapse = ";")}
        }) %>%
        dplyr::select(-flank_geneSymbols_clean) %>%
        dplyr::ungroup()  %>%
        dplyr::select(-SYMBOL, -flank_geneSymbols) %>%
        dplyr::rename(peak = rowname, SYMBOL = merged_genes) %>%
        tidyr::separate_rows(., SYMBOL, sep = ";") %>%
        dplyr::arrange(SYMBOL, peak) %>%
        base::as.data.frame()

      message("Getting the accessibility of each peaks.")
      scATACseq = Signac::RunTFIDF(scATACseq)
      atac_data = SeuratObject::GetAssayData(scATACseq, layer = "data")
      matched_peaks = peak_TG$peak[peak_TG$peak %in% base::rownames(atac_data)]
      peak_subset = atac_data[matched_peaks, ]
      non_zero_means_peaks = base::apply(peak_subset, 1, function(x) {
        non_zero_vals = x[x > 0]
        if(base::length(non_zero_vals) > 0) {base::mean(non_zero_vals)} else {0}
      })
      peak_TG$Accessibility = non_zero_means_peaks[peak_TG$peak]
      peak_TG = stats::na.omit(peak_TG)

      message("Getting the expression of each genes.")
      scRNAseq = Seurat::NormalizeData(scRNAseq)
      rna_data = SeuratObject::GetAssayData(scRNAseq, layer = "data")
      matched_genes = peak_TG$SYMBOL[peak_TG$SYMBOL %in% base::rownames(rna_data)]
      gene_subset = rna_data[matched_genes, ]
      non_zero_means_genes = base::apply(gene_subset, 1, function(x) {
        non_zero_vals = x[x > 0]
        if(base::length(non_zero_vals) > 0) {base::mean(non_zero_vals)} else {0}
      })
      peak_TG$Expression = non_zero_means_genes[peak_TG$SYMBOL]
      peak_TG = stats::na.omit(peak_TG)

      message("Filtering results.")
      peak_TG = peak_TG[!(peak_TG$Accessibility == 0 | peak_TG$Expression == 0), ]
      peak_TG = peak_TG %>% dplyr::filter(SYMBOL %in% base::as.character(TGs))

      message("Calculating peak accessibility percentile and TG expression percentile.")
      peak_TG$Accessibility_percentile = base::rank(peak_TG$Accessibility) / base::max(base::rank(peak_TG$Accessibility))
      peak_TG$Expression_percentile = base::rank(peak_TG$Expression) / base::max(base::rank(peak_TG$Expression))

      message("Calculating non zero ratio of peak accessibility.")
      non_zero_ratios_peaks = base::apply(peak_subset, 1, function(x) {
        non_zero_vals = x[x > 0]
        sum(non_zero_vals) / base::length(x)
      })
      peak_TG$Accessibility_non_zero_ratio = non_zero_ratios_peaks[peak_TG$peak]

      message("Calculating non zero ratio of TG expression.")
      non_zero_retios_genes = base::apply(gene_subset, 1, function(x) {
        non_zero_vals = x[x > 0]
        sum(non_zero_vals) / base::length(x)
      })
      peak_TG$Expression_non_zero_ratio = non_zero_retios_genes[peak_TG$SYMBOL]

      message("Calculating non zero ratio percentile of peak accessibility and non zero ratio percentile of TG expression.")
      peak_TG$Accessibility_non_zero_ratio_percentile = base::rank(peak_TG$Accessibility_non_zero_ratio) / base::max(base::rank(peak_TG$Accessibility_non_zero_ratio))
      peak_TG$Expression_non_zero_ratio_percentile = base::rank(peak_TG$Expression_non_zero_ratio) / base::max(base::rank(peak_TG$Expression_non_zero_ratio))

      message("Calculating peak importance score and peak importance score percentile.")
      peak_TG$peak_importance_score =
        weight_Accessibility_percentile * peak_TG$Accessibility_percentile +
        weight_Accessibility_non_zero_ratio_percentile * peak_TG$Accessibility_non_zero_ratio_percentile
      peak_TG$peak_importance_score_percentile = base::rank(peak_TG$peak_importance_score) / base::max(base::rank(peak_TG$peak_importance_score))

      message("Calculating TG importance score and TG importance score percentile.")
      peak_TG$TG_importance_score =
        weight_Expression_percentile * peak_TG$Expression_percentile +
        weight_Expression_non_zero_ratio_percentile * peak_TG$Expression_non_zero_ratio_percentile
      peak_TG$TG_importance_score_percentile = base::rank(peak_TG$TG_importance_score) / base::max(base::rank(peak_TG$TG_importance_score))

      message("Calculating peak-TG importance score and peak-TG importance score percentile.")
      peak_TG$peak_TG_importance_score =
        weight_Accessibility_percentile * peak_TG$Accessibility_percentile +
        weight_Expression_percentile * peak_TG$Expression_percentile +
        weight_Accessibility_non_zero_ratio_percentile * peak_TG$Accessibility_non_zero_ratio_percentile +
        weight_Expression_non_zero_ratio_percentile * peak_TG$Expression_non_zero_ratio_percentile
      peak_TG$peak_TG_importance_score_percentile = base::rank(peak_TG$peak_TG_importance_score) / base::max(base::rank(peak_TG$peak_TG_importance_score))

      message("Counting the number of TG corresponding to each peak and the number of peak corresponding to each TG, and calculating their percentile.")
      peak_TG = peak_TG %>%
        dplyr::group_by(peak) %>%
        dplyr::mutate(TG_count = dplyr::n_distinct(SYMBOL)) %>%
        dplyr::group_by(SYMBOL) %>%
        dplyr::mutate(peak_count = dplyr::n_distinct(peak)) %>%
        dplyr::ungroup()
      peak_TG$TG_count_percentile = base::rank(peak_TG$TG_count) / base::max(base::rank(peak_TG$TG_count))
      peak_TG$peak_count_percentile = base::rank(peak_TG$peak_count) / base::max(base::rank(peak_TG$peak_count))

      message("Saving results.")
      interest_cell_type_peak_TG_pred[[cell_type]][[level]][["TGs"]] = base::sort(base::unique(base::intersect(TGs, peak_TG$SYMBOL)))
      interest_cell_type_peak_TG_pred[[cell_type]][[level]][["peak_TG"]] = peak_TG
    }
  }

  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Inferring Peak-TG Relationships ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_peak_TG_pred)
}

