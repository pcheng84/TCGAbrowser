#' rnagage function
#'
#' Uses gage to find enriched pathways from the RNAseq differential gene analysis
#'
#' @param rnadeg data frame generated from rnadeg
#'
#' @import gage
#'
#' @return data frame of RNAseq values for only patients in the high and low gene expression groups
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
#' egfr_path <- rnagage(egfr_deg)
#'
#' @export
#'
rnagage <- function(deg) {
  data(kegg.gs, package = "gage")
  data(lookup)
  limma.fc <- deg$logFC
  names(limma.fc) <- lookup$entrez[match(deg$genes, lookup$gene)]
  gage.fc <- gage(limma.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
  gage.fc.g <- data.frame(pathway = rownames(gage.fc$greater), gage.fc$greater, enrichment = "greater")
  gage.fc.l <- data.frame(pathway = rownames(gage.fc$less), gage.fc$less, enrichment = "less")
  rbind(gage.fc.g, gage.fc.l)
}
