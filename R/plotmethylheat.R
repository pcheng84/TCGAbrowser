#' plotmethylheat function
#'
#' Uses ComplexHeatmap to draw heatmap
#'
#' @param mae MultiAssayExperiment object containing SimpleMethyl assay and Cohort assay from rnasubset
#' @param dmi data frame of differential methylated islands from diffmethyl function
#' @param gene character(1) Gene symbol of gene of interest
#' @param n numeric(1) number of genes to plot for the heatmap, default = 100
#'
#' @import ComplexHeatmap
#' @import circlize
#'
#' @return ComplexHeatmap of top 100 significant differentially methylated islands
#'
#' @examples
#' #using data from the curatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_ThresholdedByGene", "RPPAArray", "Methylation_methyl450"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' #prepare methylation data for differential methylation
#' lusc_t.egfr2 <- methylprep(lusc_t.egfr)
#' #find differential methlyated islands between the high and low groups
#' egfr_diffmeth <- diffmethyl(lusc_t.egfr2)
#' plotmethylheat(lusc_t.egfr2, egfr_diffmeth, "EGFR", 100)
#'
#' @export
#'
plotmethylheat <- function(mae, dmi, gene, n = 100) {
  stopifnot(any(grepl("Cohort", names(mae))))
  stopifnot(any(grepl("SimpleMethyl", names(mae))))

  meth_assay <- grep("SimpleMethyl", names(mae))
  exp_assay <- grep("Cohort", names(mae))

  mae2 <- intersectColumns(mae[, , c(meth_assay, exp_assay)])

  annot <- dcast(as.data.frame(sampleMap(mae2)), primary ~ assay, value.var = "colname")
  lvl2 <- data.frame(Cohort = colnames(mae2[[2]]), Level = mae2[[2]][1,])
  lvl3 <- merge(annot, lvl2, by = "Cohort")

  meth <- assay(mae2[rownames(dmi)[1:n], lvl3[lvl3$Level != "medium", "primary"], 1])
  grps <- factor(lvl3[lvl3$Level != "medium", "Level"], levels = c("low", "high"))

  df <- data.frame(Cell = lvl3[lvl3$Level != "medium", "Level"],
                   stringsAsFactors = FALSE)
  colnames(df) <- paste0(gene, "_group")
  cellcol <- c("#ca0020", "#0571b0")
  names(cellcol) <- c("high", "low")
  col1 <- list(Cell = cellcol)
  names(col1) <- paste0(gene, "_group")
  top_ha <- HeatmapAnnotation(df = df,
                              col = col1,
                              show_annotation_name = TRUE)
  Heatmap(meth, top_annotation = top_ha, name = "color scale",
          #col = colorRamp2(c(min(cd.t), 0, max(cd.t)), c("blue", "white", "red")),
          show_column_names = T,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 6),
          column_dend_reorder = as.numeric(df[, 1]))

}
