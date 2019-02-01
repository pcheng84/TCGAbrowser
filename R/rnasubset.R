#' rnasubset function
#'
#' Subsets MultiAssayExperiment object to one gene of interest, and add expression level column to colData.
#'
#' @param mae MultiAssayExperiment object containing data from curatedTCGAdata. Must include RNASeq2GeneNorm assay.
#' @param percent numeric(1) percentile of patients to compare 1-50
#'
#' @return Returns a MultiAssayExperiment object with column in colData identifying high, middle and low expression for one gene.
#'
#' @examples
#' data(skcm)
#' rnasubset(pat, rna, "SOX10", 10)
#'
#' @export
rnasubset <- function(mae, gene, percent) {

  #Find which assay contains RNA data (assay names have suffixes and prefixes for each cancer)
  #To do: test whether this is really a MAE object, and whether it has the appropriate assay.
  assay_num <- grep("RNASeq2GeneNorm", names(mae))

  #Extract RNASeq data for desired gene, but keep it in MAE format. This will be edited and returned.
  #To do: make sure that gene is a single character object (fail and exit if given multiple genes)
  pat.mae <- mae[gene, , assay_num]

  #Get RNA data as a matrix for calculations
  #To do: make sure that at least 3 RNA values are given.
  rna.mat <- assay(pat.mae)

  #Find RNA expression values corresponding to top and bottom percentile
  percents <- quantile(rna.mat, c(percent / 100, (1 - percent / 100)))

  #Mark samples with low, medium, and high values (in pat.mae colData) based on above percentiles
  pat.mae$level <- "middle"
  pat.mae$level[rna.mat <= percents[1]] <- "low"
  pat.mae$level[rna.mat >= percents[2]] <- "high"

  #return MAE object, now subsetted to gene of interest, with new colData column indicating expression level
  pat.mae
}
