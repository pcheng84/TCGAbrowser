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

  #Find which assay contains RNA data (assay names have suffixes and prefixes for each cancer)
  #To do: test whether this is really a MAE object, and whether it has the appropriate assay.
  assay_num <- grep("RNASeq2GeneNorm", names(mae))

  #Extract RNASeq data for all genes, but keep it in MAE format.
  #To do: make sure that gene is a single character object (fail and exit if given multiple genes)
  pat.mae <- mae[, , assay_num]

  #Merge replicates (now taking minimum, could be mean or maximum), and ensure that all samples in assay are present in colData
  pat.mae <- intersectColumns(pat.mae)
  pat.mae <- mergeReplicates(pat.mae, simplify = min)

  #Convert counts to matrix, expression levels to vector
  counts <- assay(pat.mae[, pat.mae$level != "medium", ])
  level <- pat.mae$level[pat.mae$level != "medium"]

  #Create DGEList
  deg <- DGEList(counts = counts, genes = rownames(counts), group = level)

  #Only keeps genes with at least 1 count-per-million in at least half the samples
  isexpr <- rowSums(cpm(deg)>1) >= (ncol(deg)/2)
  deg <- deg[isexpr,]

  #Create model matrix
  design <- model.matrix(~factor(level))

  #Limma - voom
  v2 <- voom(deg, design, plot = F)
  fit <- lmFit(v2, design)
  fit2 <- eBayes(fit)
  fit3 <- topTable(fit2, coef = 2, n = Inf, adjust.method = "BH", p.value = 0.05, lfc = 1,sort = "p")
  fit3
}
