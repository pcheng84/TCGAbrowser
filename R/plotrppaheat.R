#' plotrppaheat function
#'
#' Uses ComplexHeatmap to draw heatmap for RPPA data
#'
#' @param mae MultiAssayExperiment object. Must contain assays Mutation (from curatedTCGA) and Cohort (from function rnasubset())
#' @param gene character(1) Gene symbol of gene of interest
#'
#' @import data.table
#' @import ComplexHeatmap
#' @import circlize
#'
#' @return ComplexHeatmap of top 100 significant differentially expressed genes
#'
#' @examples
#' #using data from the curatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "RPPAArray"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' egfr_diffrppa <- diffrppa(lusc_t.egfr, 0.1, 0.4)
#' plotrppaheat(lusc_t.egfr, egfr_diffrppa, "EGFR", 50)
#'
#' @export
#'
plotrppaheat <- function(mae, dpe, gene, n) {
  #make sure MultiAssayExperiment object contains Cohort assay and RPPA assay
  stopifnot(any(grepl("Cohort", names(mae))))
  stopifnot(any(grepl("RPPA", names(mae))))

  #Get expression levels (which were appended by function rnasubset()), ID those with high / low expression
  rppa_assay <- grep("RPPA", names(mae))
  cohort_assay <- grep("Cohort", names(mae))

  mae2 <- intersectColumns(mae[, , c(rppa_assay, cohort_assay)])

  #create annotation table for sample matching between RPPA and Cohort samples
  annot <- dcast(as.data.frame(sampleMap(mae2)), primary ~ assay, value.var = "colname")
  lvl2 <- data.frame(Cohort = colnames(mae2[[2]]), Level = mae2[[2]][1,])
  lvl3 <- merge(annot, lvl2, by = "Cohort")

  #Extract RPPA data for all genes, subset out the "medium" expression group
  rppa <- assay(mae2[rownames(dpe)[1:n], lvl3[lvl3$Level != "medium", "primary"], 1])
  grps <- factor(lvl3[lvl3$Level != "medium", "Level"], levels = c("low", "high"))

  #Make heatmap annotation data frame
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
  Heatmap(as.matrix(rppa), top_annotation = top_ha, name = "color scale",
          #col = colorRamp2(c(min(cd.t), 0, max(cd.t)), c("blue", "white", "red")),
          show_column_names = T,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 6),
          column_dend_reorder = as.numeric(df[, 1]))

}
