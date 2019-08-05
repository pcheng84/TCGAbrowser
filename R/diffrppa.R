#' diffrppa function
#'
#' Uses limma to find differentially expressed proteins between the high and low groups
#'
#' @param mae MultiAssayExperiment object containing RNAseq assay and Cohort assay from rnasubset
#'
#'
#' @import data.table
#' @import edgeR
#' @import limma
#' @importFrom TCGAutils TCGAbarcode
#'
#' @return data frame of significant differentially expressed proteins between the two groups defined in rnasubset or mutsubset
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
#' egfr_diffrppa <- diffrppa(lusc_t.egfr)
#'
#' @export
#'
diffrppa <- function(mae) {
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
  rppa <- assay(mae2[, lvl3[lvl3$Level != "medium", "primary"], 1])
  grps <- factor(lvl3[lvl3$Level != "medium", "Level"], levels = c("low", "high"))

  #create model matrix
  design <- model.matrix(~grps)
  colnames(design) <- c("low", "high")

  #Convert RPPA data to matrix and create model matrix
  rppa.mat <- as.matrix(rppa)
  design <- model.matrix(~grps)

  # Run limma and retrun results
  fit <- lmFit(rppa.mat, design)
  fit2 <- eBayes(fit)
  fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.05, lfc=1,sort="p")
  fit3
}
