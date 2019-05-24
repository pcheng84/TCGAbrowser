#' rppadeg function
#'
#' Uses limma to find differentially expressed proteins between the high and low groups
#'
#' @param mae MultiAssayExperiment object containing RNAseq assay and Cohort assay from rnasubset
#'
#'
#' @import data.table
#' @import edgeR
#' @import limma
#'
#' @return data frame of significant differentially expressed proteins between the two groups defined in rnasubset or mutsubset
#'
#' @examples
#' #using data from the cureatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTICT", "RPPAArray"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' egfr_deg <- rppadeg(lusc_t.egfr)
#'
#' @export
#'
rppadeg <- function(mae) {
  stopifnot(any(grepl("Cohort", names(mae))))

  #Get expression levels (which were appended by function rnasubset()), ID those with high / low expression
  level_assay_num <- grep("Cohort", names(mae))
  lvl <- mae[[level_assay_num]]
  good_levels <- names(lvl[1, lvl[1,] != "medium"])

  #Extract RPPA data for all genes, subset out the "medium" expression group
  rppa_assay_num <- grep("RPPA", names(mae))
  rppa <- assay(mae[[rppa_assay_num]])
  rppa <- rppa[, TCGAbarcode(names(rppa)) %in% TCGAbarcode(good_levels)]

  #Associate expression levels with RPPA names
  levels <- lvl
  colnames(levels) <- TCGAbarcode(colnames(levels))
  levels <- levels[, TCGAbarcode(names(rppa))]

  #Convert RPPA data to matrix and create model matrix
  rppa.mat <- as.matrix(rppa)
  design <- model.matrix(~ factor(levels["level", ]))

  # Run limma and retrun results
  fit <- lmFit(rppa.mat, design)
  fit2 <- eBayes(fit)
  fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.05, lfc=1,sort="p")
  fit3
}
