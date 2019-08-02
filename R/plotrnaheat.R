#' plotrnaheat function
#'
#' Uses ComplexHeatmap to draw heatmap
#'
#' @param mae MultiAssayExperiment object containing RNAseq assay and Cohort assay from rnasubset
#' @param deg data frame of differential expressed genes from rnadeg function
#' @param gene character(1) Gene symbol of gene of interest
#' @param n numeric(1) number of genes to plot for the heatmap, default = 100
#'
#' @import ComplexHeatmap
#' @import circlize
#' @importFrom edgeR cpm
#'
#' @return ComplexHeatmap of top 100 significant differentially expressed genes
#'
#' @examples
#' #using data from the curatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_ThresholdedByGene", "RPPAArray"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' egfr_deg <- diffrna(lusc_t.egfr)
#' plotrnaheat(lusc_t.egfr, egfr_deg, "EGFR", 100)
#'
#' @export
#'
plotrnaheat <- function(mae, deg, gene, n = 100) {
  stopifnot(any(grepl("Cohort", names(mae))))

  rna_assay_num <- grep("RNASeq", names(mae))
  cohort_assay_num <- grep("Cohort", names(mae))
  lvl <- mae[[cohort_assay_num]]
  good_levels <- names(lvl[1, lvl[1,] != "medium"])

  rna <- assay(mae[[rna_assay_num]])
  rna2 <- rna[deg$genes[1:n], good_levels]

  h1 <- t(apply(log(rna2+1, 2), 1, scale))
  rownames(h1) <- deg$genes[1:n]

  df <- data.frame(as.character(lvl[1, lvl[1,] != "medium"]),
                   as.numeric(lvl[2, lvl[1,] != "medium"]), stringsAsFactors = FALSE)
  colnames(df) <- c(paste0(gene, "_group"), paste0(gene, "_expression_rank"))
  cellcol <- c("#ca0020", "#0571b0")
  names(cellcol) <- c("high", "low")
  col1 <- list(Cell = cellcol,
               expression = colorRamp2(c(1, ncol(lvl)), c("white", "purple")))
  names(col1) <- c(paste0(gene, "_group"), paste0(gene, "_expression_rank"))
  top_ha <- HeatmapAnnotation(df = df,
                              col = col1,
                              show_annotation_name = TRUE)
  Heatmap(h1, top_annotation = top_ha, name = "color scale",
          #col = colorRamp2(c(min(cd.t), 0, max(cd.t)), c("blue", "white", "red")),
          show_column_names = T,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 6),
          column_dend_reorder = as.numeric(df[, 1]))


    }
