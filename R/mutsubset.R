#' mutsubset function
#'
#' Adds mutation column to patient table
#'
#' @param pat data frame with patient data, each patient is one row
#' @param rna data frame with RNAseq count values, genes in rows, patients in columns
#' @param mut data frame with binary mutation, one gene per row, one patient per column
#' @param gene character(1) Gene symbol
#'
#' @import data.table
#'
#' @return Returns a data frame with mutation of gene as WT or mut and gene expression
#'
#' @examples
#' data(skcm)
#' mutsubset(pat, rna, mut, "SOX10")
#'
#' @export
mutsubset <- function(pat, rna, mut, gene) {
  #retrieves patients with mutation
  mutpat <- intersect(pat$bcr_patient_barcode, colnames(mut))
  setkey(pat, bcr_patient_barcode)
  pat2 <- pat[mutpat, ,mult = "first"]
  setkey(rna, Gene)
  pat2[, gene := as.numeric(mut[gene, mutpat, with = FALSE])]
  pat2[, gene2 := factor(gene, levels = 1:0, labels =c("Mutated", "WT"))]
  pat2$level <- as.numeric(rna[gene, pat2$name, with = FALSE])
  setkey(pat2, level)
  pat2$exprs_rank <- 1:nrow(pat2)
  pat2
}
