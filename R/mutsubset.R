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
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_ThresholdedByGene", "RPPAArray"), FALSE)
#' genome(lusc[[2]]) <- vapply(genome(lusc[[2]]), translateBuild, character(1L))
#' seqlevelsStyle(lusc[[2]]) <- "NCBI"
#' lusc2 <- qreduceTCGA(lusc)

#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc2, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_egfr <- mutsubset(lusc_t, "EGFR")
#'
#' @export
mutsubset <- function(mae, gene) {

  #check if mae is MultiAssayExperiment
  stopifnot(class(mae) == "MultiAssayExperiment")

  #Find which assay contains Mutation data (assay names have suffixes and prefixes for each cancer)
  assay_num <- grep("Mutation", names(mae))

  #check if gene exists in mutation dataset
  stopifnot(toupper(gene) %in% rownames(mae[[assay_num]]))

  #Get mutation data as a data.frame for calculations
  gene1 <- longFormat(mae[gene, , assay_num])

  #Mark which samples have the mutated gene and which are wild-type, save as matrix
  level <- matrix(nrow = 1, ncol = nrow(gene1), dimnames = list("level", gene1$colname))
  level["level", factor(gene1$value, levels = c(1,0), labels =c("Mutated", "WT"))]


  #Append expression level matrix to original MAE object
  mae2 <- c(mae, Cohort = level, mapFrom = assay_num)

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
