#' rnasubset function
#'
#' Annotates MultiAssayExperiment object to indicate which samples have high, medium, or low expression of a gene of interest.
#'
#' @param mae MultiAssayExperiment object containing data from curatedTCGAdata. Must include RNASeq2GeneNorm assay.
#' @param percent numeric(1) percentile of patients to compare 1-50
#'
#' @return Returns a MultiAssayExperiment object with column in colData identifying high, middle and low expression for one gene.
#'
#' @examples
#' gbm <- curatedTCGAData("gbm", "RNASeq2GeneNorm", FALSE)
#' rnasubset(gbm, "SOX10", 10)
#'
#' @export
rnasubset <- function(mae, gene, percent) {

  #Find which assay contains RNA data (assay names have suffixes and prefixes for each cancer)
  #To do: test whether this is really a MAE object, and whether it has the appropriate assay.
  assay_num <- grep("RNASeq2GeneNorm", names(mae))

  #Extract RNASeq data for desired gene, but keep it in MAE format.
  #To do: make sure that gene is a single character object (fail and exit if given multiple genes)
  pat.mae <- mae[gene, , assay_num]

  #Ensure that all samples in assay are present in colData
  pat.mae <- intersectColumns(pat.mae)

  #Get RNA data as a data.frame for calculations
  #To do: make sure that at least 3 RNA values are given.
  rna.mat <- longFormat(pat.mae)

  #Find RNA expression values corresponding to top and bottom percentile
  percents <- quantile(rna.mat$value, c(percent / 100, (1 - (percent / 100))))

  #Mark which samples have high / med / low expression, save as matrix
  level <- matrix(nrow = 1, ncol = nrow(rna.mat), dimnames = list("level", rna.mat$colname))
  level[rna.mat$value <= percents[1]] <- "low"
  level[rna.mat$value >= percents[2]] <- "high"
  level[rna.mat$value > percents[1] & rna.mat$value < percents[2]] <- "medium"

  #Append expression level matrix to original MAE object
  mae <- c(mae, ExpressionLevel = level, mapFrom = 1L)


  #return MAE object with new assay indicating expression level
  mae
}
