#' diffcp function
#'
#' Computes gene with differential copy number between two subgroups selected by `rnasubset` or `mutsubset` by Chi-squared test
#'
#' @param mae MultiAssayExperiment object. Must contain assays GISTIC or a matrix of copy number values, assay should contain "GISTIC" in the name,
#' and Cohort (from `rnasubset` or `mutsubset`)
#'
#' @importFrom MultiAssayExperiment intersectColumns
#' @import data.table
#' @importFrom reshape2 dcast
#' @return data frame with frequency of mutation in high and low group with p-value from Chi-squared test
#'
#' @examples
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_AllbyGene", "RPPAArray"), FALSE)
#' genome(lusc[[2]]) <- vapply(genome(lusc[[2]]), translateBuild, character(1L))
#' seqlevelsStyle(lusc[[2]]) <- "NCBI"
#' lusc2 <- qreduceTCGA(lusc)
#'
#' #remove normal samples and just focus on primaries
#' tnmae <- splitAssays(lusc2, c("01", "11"))
#'
#' lusc.t <- tnmae[, , grep("^01", names(tnmae))]
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc2, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_egfr <- rnasubset(lusc_t, "EGFR", 10)
#' egfr_dcp <- diffcp(lusc_egfr)
#'
#' @export
#'

diffcp <- function(mae) {
  #check only one mutation and one cohort assay in the multiassayexperiment object
  stopifnot(class(mae) == "MultiAssayExperiment", any(grepl("GISTIC", names(mae))), any(grepl("Cohort", names(mae))))
  if(length(grep("GISTIC", names(mae))) > 1)
    stop("Only one copy number assay is allowed")
  if(length(grep("Cohort", names(mae))) > 1)
    stop("Only one cohort comparison is allowed")

  # Find which assays contain mutation and expression level
  cp_assay <- grep("GISTIC", names(mae))
  exp_assay <- grep("Cohort", names(mae))

  mae2 <- intersectColumns(mae[, , c(cp_assay, exp_assay)])

  annot <- dcast(as.data.frame(sampleMap(mae2)), primary ~ assay, value.var = "colname")
  lvl2 <- data.frame(Cohort = colnames(mae2[[2]]), Level = mae2[[2]][1,])
  lvl3 <- merge(annot, lvl2, by = "Cohort")

  lvl.high <- lvl3[lvl3$Level == "high", grep("GISTIC", colnames(lvl3))]
  lvl.low <- lvl3[lvl3$Level == "low", grep("GISTIC", colnames(lvl3))]

  cp <- data.table(Gene = rownames(assay(mae2[[1]])),
                   assay(mae2[[1]]))

  cpmelt <- melt(cp, id.vars = "Gene", variable.name = "name")
  setkey(cpmelt, name)
  # simplify copy numer to 1, 0 and -1
  cpmelt[value == 2, value := 1]
  cpmelt[value == -2, value := -1]

  sum.high <- xtabs( ~ value + Gene,data = cpmelt[lvl.high])
  sum.low <- xtabs( ~ value + Gene, data = cpmelt[lvl.low])

  cpsum <- rbind(sum.high, sum.low)
  #remove genes with no losses or gains
  cpsum2 <- cpsum[,colSums(cpsum[c(1, 4), ]) != 0 | colSums(cpsum[c(3, 6),]) != 0]

  #creates data.table to be filled with p values
  cp2 <- data.table(Gene = colnames(cpsum2))

  #removes warning message from chi-squared test
  cp2[, p.value.high := apply(cpsum2, 2, function(x) {
    oopts <- options(warn = -1)
    on.exit(oopts)
    (chisq.test(matrix(x[c(2,3,5,6)], 2, byrow = TRUE))$p.value)
  })]

  cp2[, p.value.low := apply(cpsum2, 2, function(x) {
    oopts <- options(warn = -1)
    on.exit(oopts)
    (chisq.test(matrix(x[c(1,2,4,5)], 2, byrow = TRUE))$p.value)
  })]

  cp2[, FDR.high := p.adjust(p.value.high, method="BH")]
  cp2[, FDR.low := p.adjust(p.value.low, method="BH")]
  setkey(cp2, p.value.high)

  #will add locus id later
  #cp2[, Locus := cp.name$'Locus ID'[match(cp2$Gene, cp.name$'Gene Symbol')]]
}
