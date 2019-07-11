#' diffrna function
#'
#' Uses voom to find differentially expressed genes between the high and low groups
#'
#' @param mae MultiAssayExperiment object containing RNAseq assay and Cohort assay from rnasubset
#'
#' @import edgeR
#' @import limma
#'
#' @return data frame of significant differentially expressed genes between the two groups defined in rnasubset or mutsubset
#'
#' @examples
#' #using data from the cureatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_AllbyGene", "RPPAArray"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' egfr_deg <- diffrna(lusc_t.egfr)
#'
#' @export
#'
diffrna <- function(mae) {
  stopifnot(any(grepl("Cohort", names(mae))))
  #Get expression levels (which were appended by function rnasubset()), ID those with high / low expression
  level_assay_num <- grep("Cohort", names(mae))
  lvl <- mae[[level_assay_num]]
  good_levels <- names(lvl[1, lvl[1,] != "medium"])

  #Extract RNASeq data for all genes, subset out the "medium" expression group
  rna_assay_num <- grep("RNASeq", names(mae))
  counts <- assay(mae[[rna_assay_num]])
  counts <- counts[, good_levels]

  #Create DGEList
  deg <- DGEList(counts = counts, genes = rownames(counts), group = lvl[1, lvl[1,] != "medium"])

  #Only keeps genes with at least 1 count-per-million in at least half the samples
  isexpr <- rowSums(cpm(deg)>1) >= (ncol(deg)/2)
  deg <- deg[isexpr,]

  #Create model matrix
  design <- model.matrix(~factor(deg$samples$group, levels = c("low", "high")))

  #Limma - voom
  v2 <- voom(deg, design, plot = F)
  fit <- lmFit(v2, design)
  fit2 <- eBayes(fit)
  fit3 <- topTable(fit2, coef = 2, n = Inf, adjust.method = "BH", p.value = 0.05, lfc = 1,sort = "p")
  fit3
}
