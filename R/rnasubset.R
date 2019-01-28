#' rnasubset function
#'
#' Adds gene expression column and gene level columns to patient table
#'
#' @param pat data frame with patient data, each patient is one row
#' @param rna data frame with RNAseq count values, genes in rows, patients in columns
#' @param gene character(1) Gene symbol
#' @param percent numeric(1) percentile of patients to compare 1-50
#'
#' @import data.table
#'
#' @return Returns a data frame with additional columns for gene of interest expression levels and identifier for high, middle and low groups
#'
#' @examples
#' data(skcm)
#' rnasubset(pat, rna, "SOX10", 10)
#'
#' @export
rnasubset <- function(pat, rna, gene, percent) {
  #retrieves patients with RNAseq data
  # replace pat, rna with MAE object, return MAE with new "high / low" column in coldata
  pat.rna <- rna[, c("Gene", colnames(rna)[colnames(rna) %in% pat$name]), with=F]
  setkey(pat.rna, Gene)
  high <- round((ncol(pat.rna) - 1) - percent/100 * (ncol(pat.rna)-1), digits=0 )
  low <- round((ncol(pat.rna) - 1) * percent/100, digits=0)

  #selecting row for gene and only retrieving values
  pat.rna.gene <- melt.data.table(pat.rna[gene, setdiff(colnames(pat.rna), "Gene"), with=F], id.vars = NULL, measure.vars = colnames(pat.rna)[-1], variable.name = "name", value.name="level")
  setkey(pat.rna.gene, level)
  pat.rna.gene$exprs_rank <- 1:nrow(pat.rna.gene)
  #pat.rna.gene[, name := factor(name, levels=unique(name))]
  #pat.rna.gene[, ':=' (high = level > level[eval(high)], low = level <= level[eval(low)])]
  #pat.rna.gene[, high := level > level[eval(high)]]
  #pat.rna.gene[, low := level <= level[eval(low)]]
  pat.rna.gene[, `:=` (gene2 = factor(c(rep("low", times = eval(low)), rep("middle", times = eval(high - low)), rep("high", times = eval(ncol(pat.rna) - 1 - high)))))]
  phenosgene <- merge(pat, pat.rna.gene, by= "name")
  phenosgene
}
