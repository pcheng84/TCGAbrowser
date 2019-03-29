#' rnareact function
#'
#' Uses ReactomePA to find enriched pathways from the RNAseq differential gene analysis
#'
#' @param rnadeg data frame generated from rnadeg
#'
#' @import gage
#' @importFrom ReactomePA enrichPathway
#' @import DOSE
#'
#' @return data frame of enrichment for Reactome gene sets
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
#' egfr_deg <- rnadeg(lusc_t.egfr)
#' egfr_react <- rnareact(egfr_deg)
#'
#' @export
#'
rnareact <- function(deg) {
  data(lookup)
  limma.fc <- deg$logFC
  names(limma.fc) <- lookup$entrez[match(deg$genes, lookup$gene)]
  limma.names <- lookup$entrez[match(deg$genes, lookup$gene)]
  x <- enrichPathway(gene = limma.names, pvalueCutoff=1, readable=T)
}

