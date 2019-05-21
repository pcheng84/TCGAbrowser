#' mutsubset function
#'
#' Adds an assay called Cohort to a MultiAssayExperiment object to indicate which samples are mutated or wild-type for a gene of interest.
#'
#' @param mae MultiAssayExperiment object containing a simplifed Mutation assay generated from `qreduceTCGA`, must have assay with "Mutation" in the name
#' @param gene character(1) A gene to subset the mutation data
#'
#' @import MultiAssayExperiment
#'
#' @return Returns a MultiAssayExperiment object with an assay called Cohort that labels the samples according to "mutated" and "wild-type" for the gene of interest
#'
#' @examples
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTICT", "RPPAArray"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_egfr <- mutsubset(lusc_t, "EGFR")
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
