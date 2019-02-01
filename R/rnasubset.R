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

  #Merge replicates (now taking minimum, could be mean or maximum), and ensure that all samples in assay are present in colData
  pat.mae <- intersectColumns(pat.mae)
  pat.mae <- mergeReplicates(pat.mae, simplify = min)

  #Get RNA data as a data.frame for calculations
  #To do: make sure that at least 3 RNA values are given.
  rna.mat <- longFormat(pat.mae)

  #Find RNA expression values corresponding to top and bottom percentile
  percents <- quantile(rna.mat$value, c(percent / 100, (1 - (percent / 100))))

  #Find which samples have low and high values based on above percentiles
  low.samps <- rna.mat$primary[rna.mat$value <= percents[1]]
  high.samps <- rna.mat$primary[rna.mat$value >= percents[2]]
  med.samps <- rna.mat$primary[rna.mat$value > percents[1] & rna.mat$value < percents[2]]

  #Mark low, medium, and high expression values in original MAE colData. Mark samples NA if they have no RNA data.
  mae$level <- NA
  mae$level[mae$patientID %in% med.samps] <- "medium"
  mae$level[mae$patientID %in% low.samps] <- "low"
  mae$level[mae$patientID %in% high.samps] <- "high"

  #return MAE object with new colData column indicating expression level
  mae
}
