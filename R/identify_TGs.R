#' @title Identify Target Genes via Chromatin Accessibility and Annotation
#' @description
#' This function performs genomic annotation on chromatin accessibility peaks (from scATAC-seq) to identify potential target genes (TGs) with open promoter.
#' It uses the ChIPseeker package for peak annotation
#' and generates output files including annotation plots (peak_anno.png), annotation results (peak_anno.csv), and identified TGs (TGs.csv).
#'
#' @note
#' The scRNA-seq data and scATAC-seq data need to be loaded as objects named scRNAseq and scATACseq,
#' respectively, and they must come from the same cells with identical cell labels.
#'
#' @param genome_for_anno
#' A TxDb object for genome annotation (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene, TxDb.Hsapiens.UCSC.hg19.knownGene).
#' Default is TxDb.Hsapiens.UCSC.hg38.knownGene.
#' @param promoter_range A positive numeric object. Defines the promoter region around TSS (in base pairs). Default is 50 bp.
#' @param scRNAseq A Seurat object containing scRNA-seq data (the initially loaded data) and gene names as rownames.
#' @param scATACseq A Seurat object containing scATAC-seq data (the initially loaded data) with files _fragments.tsv.gz.tbi and _fragments.tsv.gz (in Signac framework).
#' @param annoDb "org.Hs.eg.db". Specifies the organism-level annotation package.
#'
#' @returns
#' A list containing two elements:
#'      peak_anno_promoter: A data frame of annotated peaks.
#'     TGs: A character vector of unique target gene symbols identified.
#' @export
#'
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' # Load scRNA-seq data and scATAC-seq data as object scRNAseq and scATACseq with the same cells and cell labels.
#' # scRNAseq = SeuratData::LoadData("pbmcMultiome", "pbmc.rna")
#' # scATACseq = SeuratData::LoadData("pbmcMultiome", "pbmc.atac")
#' scRNAseq = base::readRDS("./scRNAseq.rds")
#  scATACseq = base::readRDS("./scATACseq.rds")
#' promoter_range = 50
#' results_identify_TGs = identify_TGs(
#'   genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   promoter_range = promoter_range,
#'   scRNAseq = scRNAseq,
#'   scATACseq = scATACseq,
#'   annoDb = "org.Hs.eg.db"
#' )
#' base::saveRDS(results_identify_TGs, file = "./1.1 Data Preprocessing - Identify TGs By Annotation/results_identify_TGs.rds")
#' }
identify_TGs = function(genome_for_anno = TxDb.Hsapiens.UCSC.hg38.knownGene,
                        promoter_range = 50,
                        scRNAseq = scRNAseq,
                        scATACseq = scATACseq,
                        annoDb = "org.Hs.eg.db")
{
  ### Start.
  t_start = base::Sys.time()
  message("Run: Annotating Chromatin Fragments and Identifying TGs ", t_start, ".")
  original_dir = base::getwd()
  target_folder = "1.1 Data Preprocessing - Identify TGs By Annotation"
  if (base::dir.exists(target_folder)) {
    base::unlink(target_folder, recursive = TRUE, force = TRUE)
    message("The original folder has been deleted: ", target_folder, ".")
  }
  base::dir.create(target_folder, showWarnings = FALSE, recursive = TRUE)
  message("The new folder has been created: ", target_folder, ".")
  base::setwd(target_folder)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  message("Annotate chromatin fragments using the annotatePeakfunction from the ChIPseeker package.")
  message("The genome annotation database used is ", base::unique(genome_for_anno$user_genome), ".")

  ### Peak annotation.
  peak_anno = ChIPseeker::annotatePeak(
    peak = scATACseq@assays[["ATAC"]]@ranges,
    tssRegion = c(-promoter_range, promoter_range),
    TxDb = genome_for_anno,
    annoDb = annoDb,
    assignGenomicAnnotation = TRUE,
    flankDistance = promoter_range,
    sameStrand = FALSE,
    ignoreOverlap = TRUE
  )

  ### Savethe annotation results.
  message("Save the annotation results.")
  ggplot2::ggsave("peak_anno.png", plot = ChIPseeker::plotAnnoBar(peak_anno), width = 10, height = 6)
  peak_anno = base::as.data.frame(peak_anno) %>% dplyr::filter(!is.na(ENSEMBL)) %>% dplyr::filter(!is.na(SYMBOL)) %>% dplyr::filter(!is.na(annotation))
  utils::write.csv(peak_anno, file = "peak_anno.csv", row.names = FALSE)
  base::rownames(peak_anno) = IRanges::paste(peak_anno$seqnames, peak_anno$start, peak_anno$end, sep = "-")

  ### Obtain TGs.
  message("Save TGs.")
  peak_anno_promoter = peak_anno %>% dplyr::filter(annotation == "Promoter")
  TGs = base::unique(base::intersect(peak_anno_promoter$SYMBOL, base::rownames(scRNAseq)))
  message(length(TGs), " target genes have been detected, ", base::sum(TGs %in% base::rownames(scRNAseq)), " of which are present in the RNA data.")
  tryCatch(
    utils::write.csv(TGs, file = "TGs.csv", row.names = FALSE), error = function(e) {
      stop("The file could not be saved: ", e$message)
    }
  )

  ### Update the annotation file.
  peak_anno_promoter = peak_anno_promoter %>% dplyr::filter(SYMBOL %in% TGs)

  ### End.
  base::setwd(original_dir)
  message("The current working directory has been switched to: ", base::getwd(), ".")
  t_end = base::Sys.time()
  message("Finish: Annotating Chromatin Fragments and Identifying TGs ", t_end, ".")
  message("Running time: ")
  base::print(t_end - t_start)
  return(base::list("peak_anno_promoter" = peak_anno_promoter, "TGs" = TGs))
}
