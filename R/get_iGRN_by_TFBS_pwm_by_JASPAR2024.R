#' @title Construct Standard Gene Regulatory Networks (iGRNs) from Transcription Factor Binding Sites (TFBS) Position Weight Matrices PWMs
#'
#' @description
#' This function constructs standard gene regulatory networks (iGRNs) by predicting transcription factor binding sites (TFBS) using position weight matrices (PWMs) from the JASPAR2024 database.
#' It maps transcription factors (TFs) to target genes (TGs) based on PWM matching scores in promoter regions,
#' generating TF-TG association matrices for each cell type.
#'
#' @param interest_cell_type_genes_pseudotime_info The output of function get_genes_pseudotime_info
#' @param promoter_range See in ?identify_TGs.
#' @param results_identify_TGs The output of function identify_TGs.
#' @param genome
#' A BSgenome object representing the reference genome to use for sequence extraction.
#' Default is \code{BSgenome.Hsapiens.UCSC.hg38} (R package BSgenome.Hsapiens.UCSC.hg38 is needed).
#' Also \code{BSgenome.Hsapiens.UCSC.hg19} is valid  (R package BSgenome.Hsapiens.UCSC.hg19 is needed).
#' @param min_score_for_matchPWM
#' Character. Minimum score threshold for PWM matching.
#' Can be a percentage string between 0 and 1 (e.g., "80%"). Default is "80%".
#' @param species See in ?get_TGs_from_JASPAR2024.
#' @param collection See in ?get_TGs_from_JASPAR2024.
#' @param output_predicted_TFBS Logical. Whether to output detailed TFBS prediction information. Default is FALSE.
#' @param ncores See in ?get_interest_cell_type_data.
#'
#' @returns
#' A named list (by cell type) where each element is a matrix (its class is data.frame) representing the TF-TG association network:
#' \itemize{
#' \item Rows: Target genes (TGs)
#' \item Columns: Transcription factors (TFs)
#' \item Values: Normalized PWM match scores (ranging from 0 to approximately 0.95)
#' }
#' The matrices are sparse, with most values being 0 indicating no predicted regulatory relationship.
#'
#' @details
#' The function performs the following key operations:
#' \enumerate{
#' \item Maps TGs to their promoter regions using chromatin accessibility data
#' \item Retrieves PWMs for known TFs from JASPAR2024 database
#' \item Normalizes PWM matrices to represent base probabilities
#' \item Uses parallel processing to match PWMs to promoter sequences in Linux
#' \item Constructs TF-TG association matrices based on match scores
#' \item Normalized results are divided by 1.05 times the maximum match score to prevent match scores that are too close to 1 from causing instability in subsequent training
#' }
#'
#' @section Algorithm Details:
#' For each cell type, the function:
#' \itemize{
#' \item Extracts promoter sequences for all target genes
#' \item Matches each TF's PWM against all promoter sequences
#' \item Records the maximum match score for each TF-TG pair
#' \item Filters TFs to those present in the expression data
#' \item Constructs sparse association matrices
#' }
#'
#' @note
#' Important considerations:
#' \itemize{
#' \item \code{min_score_for_matchPWM} is a critical parameter that affects network sparsity
#' \item Only works on Linux systems for parallel processing due to \code{mclapply} limitations
#' \item PWM matching is computationally intensive - use parallel processing when possible
#' \item The function generates a CSV file with detailed TFBS predictions in the output directory
#' \item Normalization (division by 1.05 times max score) ensures values range from 0 to 0.95
#' \item Requires BSgenome and JASPAR2024 packages for sequence data and PWMs
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' interest_cell_type_genes_pseudotime_info = base::readRDS("./2.2 Data Processing - Cell Grouping/interest_cell_type_genes_pseudotime_info.rds")
#' promoter_range = 50
#' results_identify_TGs = base::readRDS("./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' min_score_for_matchPWM = "80%"
#' ncores = ceiling(parallel::detectCores() / 2) # in Linux
#' # ncores = 1 # in Windows
#' interest_cell_type_iGRN_all_TGTF_pairs = get_iGRN_by_TFBS_pwm_by_JASPAR2024(
#'   interest_cell_type_genes_pseudotime_info = interest_cell_type_genes_pseudotime_info,
#'   promoter_range = promoter_range,
#'   results_identify_TGs = results_identify_TGs,
#'   genome = BSgenome.Hsapiens.UCSC.hg38,
#'   min_score_for_matchPWM = min_score_for_matchPWM,
#'   species = "Homo sapiens",
#'   collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
#'   output_predicted_TFBS = FALSE,
#'   ncores = ncores
#' )
#' base::saveRDS(interest_cell_type_iGRN_all_TGTF_pairs, file = "./3 get iGRN/interest_cell_type_iGRN_all_TGTF_pairs.rds")
#' }
get_iGRN_by_TFBS_pwm_by_JASPAR2024 = function(
    interest_cell_type_genes_pseudotime_info  = interest_cell_type_genes_pseudotime_info,
    promoter_range = 50,
    results_identify_TGs = results_identify_TGs,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    min_score_for_matchPWM = "80%",
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                   "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    output_predicted_TFBS = FALSE,
    ncores = 1
)
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Constructing Standard Gene Regulatory Networks (GRNs) from Transcription Factor Binding Site (TFBS) Position Weight Matrices (PWMs) ", t_start, ".")
  # message("Retrieve known TFs from JASPAR2024 and assume all TFs have potential regulatory relationships with each TG.")
  message("Only considering TFs and TGs that appear in pseudotime analysis.")
  message("Generating a TG-TF association matrix with known TFs filtered by TF-TFBS match strength (rows: TGs, columns: TFs, regulatory strength based on TF-TFBS match scores).")
  original_dir = base::getwd()
  new_folder = "3 get iGRN"
  if (!base::dir.exists(new_folder)) {
    base::dir.create(new_folder, recursive = TRUE)
    message("Folder already creates: ", new_folder, ".")
  } else {message("Folder already exists: ", new_folder, ".")}
  base::setwd(new_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")

  ### Defineing inner function.
  get_iGRN_by_TFBS_pwm_by_JASPAR2024_process = function(
    TGs,
    TFs,
    promoter_range,
    peak_anno_promoter,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    min_score_for_matchPWM = "80%",
    ncores = 1,
    species = "Homo sapiens",
    collection = c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII",
                   "FAM", "PBM", "PBM_HOMEO", "PBM_HLH", "UNVALIDATED"),
    output_predicted_TFBS = FALSE)
  {
    # Constructing TG-Promoter mapping (splitting rows without strand information into two rows: positive strand and negative strand).
    map_from_TGs_to_promoter = function(peak_anno_promoter)
    {
      message("Constructing TG-Promoter mapping (splitting rows without strand information into two rows: positive strand and negative strand).")
      TGs_and_promoter = base::data.frame(
        ensembl_gene_id = peak_anno_promoter$ENSEMBL,
        chromosome_name = peak_anno_promoter$seqnames,
        start_position = peak_anno_promoter$start,
        end_position = peak_anno_promoter$end,
        strand = peak_anno_promoter$strand,
        TGs = peak_anno_promoter$SYMBOL,
        promoter_peak = base::rownames(peak_anno_promoter),
        width = peak_anno_promoter$width,
        distanceToTSS = peak_anno_promoter$distanceToTSS,
        stringsAsFactors = FALSE
      )
      TGs_and_promoter$strand = base::as.character(TGs_and_promoter$strand)
      TGs_and_promoter_withnout_strand = TGs_and_promoter[TGs_and_promoter$strand == "*", ]
      TGs_and_promoter$strand[TGs_and_promoter$strand == "*"] = "1"
      TGs_and_promoter_withnout_strand$strand = base::as.character(TGs_and_promoter_withnout_strand$strand)
      TGs_and_promoter_withnout_strand$strand[TGs_and_promoter_withnout_strand$strand == "*"] = "-1"
      if (base::nrow(TGs_and_promoter_withnout_strand) > 0) {
        TGs_and_promoter = base::rbind(TGs_and_promoter, TGs_and_promoter_withnout_strand)
      }
      TGs_and_promoter$chromosome_name = base::gsub("^chr", "", TGs_and_promoter$chromosome_name)
      TGs_and_promoter$chromosome_name = base::as.character(TGs_and_promoter$chromosome_name)
      return(TGs_and_promoter)
    }
    TGs_and_promoter = map_from_TGs_to_promoter(peak_anno_promoter)
    TGs_and_promoter = TGs_and_promoter[TGs_and_promoter$TGs %in% TGs, , drop = FALSE]

    # Defining a function to extract promoter sequences.
    get_promoter_sequence = function(chromosome, start, end, strand,
                                     upstream_distance = promoter_range,
                                     genome = genome)
    {
      if (strand == 1) {
        promoter_start = base::max(start - upstream_distance, 1)
        promoter_end = start - 1
      } else {
        promoter_start = end + 1
        promoter_end = base::min(
          end + upstream_distance,
          GenomeInfoDb::seqlengths(genome)[[base::paste0("chr", chromosome)]]
        )
      }
      promoter_seq =  Biostrings::getSeq(genome, base::paste0("chr", chromosome), promoter_start, promoter_end)
      return(promoter_seq)
    }

    # Retrieving PWMs of TFBS from JASPAR database (typical TF binding length: 5-20 bp)
    message("Retrieving PWMs of TFBS from JASPAR database (typical TF binding length: 5-20 bp).")
    pwm = TFBSTools::getMatrixSet(
      x = JASPAR2024::db(JASPAR2024::JASPAR2024()),
      opts = base::list(all = F, species = species, collection = collection,
                        all_versions = FALSE, matrixtype = "PWM")
    )

    # Retaining only TFs present in pseudotime analysis data.
    message("Retaining only TFs present in  pseudotime analysis data.")
    pwm@listData = base::Filter(function(x) x@name %in% TFs, pwm@listData)

    # Normalize each column of the position weight matrix to sum to 1 (rows: nucleotide bases, columns: binding site positions, values: base occurrence probabilities, column sum constraint: 1).
    message("Normalize each column of the position weight matrix to sum to 1 (rows: nucleotide bases, columns: binding site positions, values: base occurrence probabilities, column sum constraint: 1).")
    for (i in base::seq_along(pwm@listData)) {
      pwm@listData[[i]]@profileMatrix = 2 ^ (pwm@listData[[i]]@profileMatrix) / 4
    }

    # Calculating the sequence of each promoter.
    # Using parallel::mclapply for parallel computing (Only effective on Linux systems).
    message("Calculating the sequence of each promoter.")
    message("The number of CPU cores for parallel::mclapply parallel computing is: ", ncores)
    promoter_seq = base::list()
    promoter_seq <- parallel::mclapply(1:base::dim(TGs_and_promoter)[1], function(i) {
      get_promoter_sequence(TGs_and_promoter[i, "chromosome_name"],
                            TGs_and_promoter[i, "start_position"],
                            TGs_and_promoter[i, "end_position"],
                            TGs_and_promoter[i, "strand"],
                            genome = genome)
    }, mc.cores = ncores)

    # Processing each PWM and matching to TFGs based on PWM patterns.
    # Using parallel::mclapply for parallel computing (Only effective on Linux systems).
    message("Processing each PWM and matching to TFGs based on PWM patterns.")
    message("The minimum score threshold for PWM matching to TFGs is: ", min_score_for_matchPWM)
    message("The number of CPU cores for parallel::mclapply parallel computing is: ", ncores)
    results_PWM_match_TF_by_parallel = base::list()
    length_pwm = base::length(base::names(pwm))
    for (name in base::names(pwm)) {
      pos = base::match(name, base::names(pwm))
      if (pos %% 10 == 1 | pos == length_pwm) {
        message(pos, "th in ", length_pwm)
      }
      results = parallel::mclapply(1:base::dim(TGs_and_promoter)[1], function(i) {
        predicted_TFBS = Biostrings::matchPWM(pwm = pwm[[name]]@profileMatrix,
                                              subject = promoter_seq[[i]],
                                              min.score = min_score_for_matchPWM,
                                              with.score = TRUE)
        gene_name = TGs_and_promoter$TGs[i]
        if (!base::is.null(predicted_TFBS) &&
            base::length(gene_name) != 0 &&
            !base::is.na(mean(predicted_TFBS@elementMetadata@listData[["score"]]))) {
          if (output_predicted_TFBS == T) {
            result_line = base::paste(gene_name,
                                      pwm[[name]]@name,
                                      predicted_TFBS,
                                      base::max(predicted_TFBS@elementMetadata@listData[["score"]]),
                                      sep = ",")
          } else {
            result_line = base::paste(gene_name,
                                      pwm[[name]]@name,
                                      "-",
                                      base::max(predicted_TFBS@elementMetadata@listData[["score"]]),
                                      sep = ",")
          }
          return(result_line)
        } else {
          return(NULL)
        }
      }, mc.cores = ncores)
      results = base::unlist(base::Filter(base::Negate(is.null), results))

      results_PWM_match_TF_by_parallel[[name]] = results
    }

    message("Saving parallel computation results to 3 get iGRN/results_PWM_match_TF_by_parallel.rds")
    base::saveRDS(results_PWM_match_TF_by_parallel, file = "results_PWM_match_TF_by_parallel.rds")

    # Writing results to CSV file.
    message("Saving parallel computation results to 3 get iGRN/TFBS.csv")
    file_conn = base::file("TFBS.csv", "w")
    base::writeLines("TGs, Potential_TFs, Predicted_TFBS, match_score", file_conn)
    for (result_list in results_PWM_match_TF_by_parallel) {
      base::writeLines(as.character(result_list), file_conn)
    }
    base::close(file_conn)

    # Generating TF-TG association matrix (rows: TGs, columns: TFs, values: match scores normalized by 1.05 times max match score).
    message("Generating TF-TG association matrix (rows: TGs, columns: TFs, values: match scores normalized by 1.05 times max match score).")
    potential_TFs = data.table::fread("TFBS.csv")
    potential_TFs = potential_TFs[potential_TFs$Predicted_TFBS != "" & !base::is.na(potential_TFs$Predicted_TFBS), ]
    # potential_TFs = potential_TFs[potential_TFs$Potential_TFs %in% base::rownames(scRNAseq), ]
    potential_TFs = potential_TFs[potential_TFs$Potential_TFs %in% TFs, ]
    potential_TFs = potential_TFs[, -c(3)]
    # library(dplyr)
    potential_TFs = potential_TFs %>%
      dplyr::group_by(TGs, Potential_TFs) %>%
      dplyr::summarise(match_score = base::mean(match_score)) %>%
      dplyr::ungroup()
    M_TGs_TFs = tidyr::pivot_wider(potential_TFs, names_from = Potential_TFs, values_from = match_score, values_fill = 0)
    M_TGs_TFs = tibble::column_to_rownames(M_TGs_TFs, "TGs")
    M_TGs_TFs = M_TGs_TFs / (base::max(M_TGs_TFs) * 1.05)

    message("Generating an association matrix containing ", nrow(M_TGs_TFs), " TGs and", ncol(M_TGs_TFs), " TFs with a sparsity of ", base::round(base::sum(M_TGs_TFs == 0) / (base::nrow(M_TGs_TFs) * base::ncol(M_TGs_TFs)), 3))
    return(M_TGs_TFs)
  }

  ### Running inner function.
  interest_cell_type_iGRN_all_TGTF_pairs = parallel::mclapply(
    X = base::names(interest_cell_type_genes_pseudotime_info),
    FUN = function(cell_type) {
      message("Preparing to construct iGRN for cell type ", cell_type, ".")
      get_iGRN_by_TFBS_pwm_by_JASPAR2024_process(
        TGs = base::sort(base::unique(base::Reduce(
          base::union,
          interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TGs"]]
        ))),
        TFs = base::sort(base::unique(base::Reduce(
          base::union,
          interest_cell_type_genes_pseudotime_info[[cell_type]][["Branches_TFs"]]
        ))),
        promoter_range = promoter_range,
        peak_anno_promoter = results_identify_TGs[["peak_anno_promoter"]],
        genome = genome,
        min_score_for_matchPWM = min_score_for_matchPWM,
        ncores = ncores,
        species = species,
        collection = collection,
        output_predicted_TFBS = output_predicted_TFBS
      )},
    mc.cores = 1)
  base::names(interest_cell_type_iGRN_all_TGTF_pairs) = base::names(interest_cell_type_genes_pseudotime_info)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Constructing Standard Gene Regulatory Networks (GRNs) from Transcription Factor Binding Site (TFBS) Position Weight Matrices (PWMs) ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(interest_cell_type_iGRN_all_TGTF_pairs)
}
