#' diffmethyl function
#'
#' Uses limma to find differentially methylated islands between the high and low groups
#'
#' @param mae MultiAssayExperiment object containing Methylation assay and Cohort assay from rnasubset
#'
#' @import data.table
#' @import edgeR
#' @import limma
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#'
#' @return data frame of significant differentially methylated islands between the two groups defined in rnasubset or mutsubset
#'
#' @examples
#' #using data from the cureatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_AllbyGene", "RPPAArray", "Methylation_methyl450"), FALSE)
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' #subset LUSC data by expression of EGFR
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' #prepare methylation data for differential methylation
#' methylprep(lusc_t.egfr)
#' #find differential methlyated islands between the high and low groups
#' egfr_diffmeth <- diffmethyl(lusc_t.egfr)
#'
#' @export
#'
diffmethyl <- function(mae) {
  stopifnot(any(grepl("Cohort", names(mae))))
  stopifnot(any(grepl("SimpleMethyl", names(mae))))

  meth_assay <- grep("SimpleMethyl", names(mae))
  exp_assay <- grep("Cohort", names(mae))

  mae2 <- intersectColumns(mae[, , c(meth_assay, exp_assay)])

  annot <- dcast(as.data.frame(sampleMap(mae2)), primary ~ assay, value.var = "colname")
  lvl2 <- data.frame(Cohort = colnames(mae2[[2]]), Level = mae2[[2]][1,])
  lvl3 <- merge(annot, lvl2, by = "Cohort")

  lvl.high <- lvl3[lvl3$Level == "high", grep("Mutation", colnames(lvl3))]
  lvl.low <- lvl3[lvl3$Level == "low", grep("Mutation", colnames(lvl3))]

  meth <- assay(mae2[[1]])

  #Get expression levels (which were appended by function rnasubset()), ID those with high / low expression
  level_assay_num <- grep("Cohort", names(mae))
  lvl <- mae[[level_assay_num]]
  good_levels <- names(lvl[1, lvl[1,] != "medium"])

  #Associate expression levels with methylation names
  levels <- lvl
  colnames(levels) <- TCGAbarcode(colnames(levels))
  levels <- levels[, TCGAbarcode(colnames(meth))]

  #convert median methylation data to matrix and create model matrix
  meth.mat <- as.matrix(methyltest7[, 2:ncol(methyltest7)])
  row.names(meth.mat) <- methyltest7$Island_Name
  design <- model.matrix(~ factor(levels["level", ]))

  #Perform DGE analysis on asin transformed data and return result
  fit <- lmFit(asin(meth.mat), design)
  fit2 <- eBayes(fit, trend = TRUE)
  fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.05, lfc = 0.1, sort="p")
  fit3

}
