#' rppaheat function
#'
#' Uses ComplexHeatmap to draw heatmap for RPPA data
#'
#' @param mae MultiAssayExperiment object. Must contain assays Mutation (from curatedTCGA) and ExpressionLevel (from function rnasubset())
#' @param gene character(1) Gene symbol of gene of interest
#'
#' @import data.table
#' @import ComplexHeatmap
#' @import circlize
#' @importFrom edgeR cpm
#' @importFrom TCGAutils TCGAbarcode
#'
#' @return ComplexHeatmap of top 100 significant differentially expressed genes
#'
#' @examples
#' mae <- curatedTCGAData("SKCM", c("RNASeq2GeneNorm", "Mutation", "GISTIC", "RPPAArray") , FALSE)
#' gene <- "SOX10"
#' sox10 <- rnasubset(mae, gene, 10)
#' rppaheat(sox10, gene)
#' @export
#'
rppaheat <- function(mae, gene) {
  #Find which assays contain expression levels, RPPA data, and RNASeq data
  exp_assay <- grep("Expression", names(mae))
  rppa_assay <- grep("RPPA", names(mae))
  rna_assay <- grep("RNASeq2GeneNorm", names(mae))

  #Get names of samples with high and low expression, convert them to common format
  high <- TCGAbarcode(colnames(mae[[exp_assay]])[mae[[exp_assay]] == "high"])
  low  <- TCGAbarcode(colnames(mae[[exp_assay]])[mae[[exp_assay]] == "low"])

  #Find the patients in RPPA assay corresponding to above
  overlap_high <- TCGAbarcode(colnames(mae[[rppa_assay]])) %in% high
  overlap_high <- colnames(mae[[rppa_assay]])[overlap_high]

  overlap_low <- TCGAbarcode(colnames(mae[[rppa_assay]])) %in% low
  overlap_low <- colnames(mae[[rppa_assay]])[overlap_low]
  all_overlaps <- c(overlap_high, overlap_low)

  #Extract RPPA data, limit to samples ID'd above, convert to matrix, get complete cases, convert names to common format
  counts <- mae[[rppa_assay]]
  counts <- counts[, all_overlaps]
  counts <- as.matrix(assay(counts))
  counts <- counts[complete.cases(counts), ]
  colnames(counts) <- TCGAbarcode(colnames(counts))

  #Make heatmap annotation data frame, 2 columns, first is high/low, second is gene expression
  df <- data.frame(c(rep("high", times = length(overlap_high)),
                     rep("low" , times = length(overlap_low))),
                   c(assay(mae[gene, TCGAbarcode(overlap_high), rna_assay]),
                     assay(mae[gene, TCGAbarcode(overlap_low ), rna_assay])))
  colnames(df) <- c(paste0(gene, "_group"), paste0(gene, "_expression"))

  #make colors, name them high/low
  cellcol <- c("#ca0020", "#0571b0")
  names(cellcol) <- c("high", "low")
  col1 <- list(Cell = cellcol,
               expression = colorRamp2(range(df[, 2]), c("white", "purple")))
  names(col1) <- c(paste0(gene, "_group"), paste0(gene, "_expression"))

  #Complete heatmap annotation object
  top_ha <- HeatmapAnnotation(df = df, col = col1)

  #Draw and decorate heatmap
  x <- Heatmap(counts, top_annotation = top_ha, name = "color scale",
               show_column_names = T,
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 6),
               column_dend_reorder = as.numeric(df[, 1]))

  draw(x)
  for(an in colnames(df)) {
    decorate_annotation(an, {
      # annotation names on the right
      grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
      # annotation names on the left
      #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
    })
  }

}
