#' rnasubset function
#'
#' Annotates MultiAssayExperiment object to indicate which samples have high, medium, or low expression of a gene of interest.
#'
#' @param mae MultiAssayExperiment object containing RNAseq data, must have assay with "RNAseq" in the name
#' @param gene character(1) A gene to subset the RNAseq data
#' @param percent numeric(1) percentile of patients to compare 1-50
#'
#' @import MultiAssayExperiment
#'
#' @return Returns a MultiAssayExperiment object with an assay called Cohort that labels the samples according to "high", "middle" and "low" expression for the gene of interest
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
#' lusc_egfr <- rnasubset(lusc_t, "EGFR", 10)
#'
#' @export
rnasubset <- function(mae, gene, percent) {

  #check if mae is MultiAssayExperiment
  stopifnot(class(mae) == "MultiAssayExperiment")

  #Find which assay contains RNA data (assay names have suffixes and prefixes for each cancer)
  assay_num <- grep("RNASeq2GeneNorm", names(mae))

  #check if gene exists in RNAseq dataset
  stopifnot(toupper(gene) %in% rownames(mae[[assay_num]]))

  #Get RNA data as a data.frame for calculations
  gene1 <- longFormat(mae[gene, , assay_num])

  percent <- 10
  #Find RNA expression values corresponding to top and bottom percentile
  percents <- quantile(gene1$value, c(percent / 100, (1 - (percent / 100))))

  #Mark which samples have high / med / low expression, save as matrix
  level <- matrix(nrow = 3, ncol = nrow(gene1), dimnames = list(c("level", "rank", "counts"), gene1$colname))
  level["level", gene1$value <= percents[1]] <- "low"
  level["level", gene1$value >= percents[2]] <- "high"
  level["level", gene1$value > percents[1] & gene1$value < percents[2]] <- "medium"
  level["rank",] <- as.numeric(rank(gene1$value))
  level["counts", ] <- gene1$value

  #Append expression level matrix to original MAE object
  mae2 <- c(mae, Cohort = level, mapFrom = assay_num)

  #return MAE object with new assay indicating expression level
  mae2
}
