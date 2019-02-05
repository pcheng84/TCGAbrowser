#' rnadeg function
#'
#' Uses voom to find differentially expressed genes between the high and low groups
#'
#' @param mae MultiAssayExperiment object containing RNASeq2GeneNorm assay, with gene expression level marked by rnasubset()
#'
#' @import edgeR
#' @import limma
#'
#' @return data frame of significant differentially expressed genes between the two groups defined in rnasubset or mutsubset
#'
#' @examples
#' gbm <- curatedTCGAData("gbm", "RNASeq2GeneNorm", FALSE)
#' gene <- "SOX10"
#' sox10.pat <- rnasubset(gbm, gene, 10)
#' sox10.deg <- rnadeg(sox10.pat)
#'
#' @export
#'
rnadeg <- function(mae) {

  #Get expression levels (which were appended by function rnasubset()), ID those with high / low expression
  level_assay_num <- grep("ExpressionLevel", names(mae))
  levels <- mae[[level_assay_num]]
  good_levels <- colnames(levels)[levels != "medium"]

  #Extract RNASeq data for all genes, subset out the "medium" expression group
  rna_assay_num <- grep("RNASeq2GeneNorm", names(mae))
  counts <- assay(mae[[rna_assay_num]])
  counts <- counts[, good_levels]

  #Create DGEList
  deg <- DGEList(counts = counts, genes = rownames(counts), group = levels[good_levels])

  #Only keeps genes with at least 1 count-per-million in at least half the samples
  isexpr <- rowSums(cpm(deg)>1) >= (ncol(deg)/2)
  deg <- deg[isexpr,]

  #Create model matrix
  design <- model.matrix(~factor(levels[, good_levels]))

  #Limma - voom
  v2 <- voom(deg, design, plot = F)
  fit <- lmFit(v2, design)
  fit2 <- eBayes(fit)
  fit3 <- topTable(fit2, coef = 2, n = Inf, adjust.method = "BH", p.value = 0.05, lfc = 1,sort = "p")
  fit3
}
