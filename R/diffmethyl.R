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
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_ThresholdedByGene", "RPPAArray", "Methylation_methyl450"), FALSE)
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' #subset LUSC data by expression of EGFR
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' #prepare methylation data for differential methylation
#' lusc_t.egfr2 <- methylprep(lusc_t.egfr)
#' #find differential methlyated islands between the high and low groups
#' egfr_diffmeth <- diffmethyl(lusc_t.egfr2)
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

  meth <- assay(mae2[, lvl3[lvl3$Level != "medium", "primary"], 1])
  grps <- factor(lvl3[lvl3$Level != "medium", "Level"], levels = c("low", "high"))
  #create model matrix
  design <- model.matrix(~grps)
  colnames(design) <- c("low", "high")
  #Perform DGE analysis using lmFit and return table with intercept and pvalue
  fit <- lmFit(meth, design)

  fit2 <- eBayes(fit)
  t1 <- data.frame(intercept = fit2$coefficients[,1],
                    f = fit2[["F"]],
                    pval = fit2[["F.p.value"]],
                   adj.pval = p.adjust(fit2[["F.p.value"]], method = "BH"))
  t1

}
