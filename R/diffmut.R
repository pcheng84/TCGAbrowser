#' diffmut function
#'
#' Computes differentially mutated genes between high and low group by Chi-squared test
#'
#' @param mae MultiAssayExperiment object. Must contain assays Mutation (from curatedTCGA) and ExpressionLevel (from function rnasubset())
#'
#' @import data.table
#' @import TCGAUtils
#' @return data frame with frequency of mutation in high and low group with p-value from Chi-squared test
#'
#' @examples
#' mae <- curatedTCGAData("SKCM", c("RNASeq2GeneNorm", "Mutation", "GISTIC", "RPPAArray") , FALSE)
#' gene <- "SOX10"
#' sox10 <- rnasubset(mae, gene, 10)
#' diffmut(sox10)
#'
#' @export
#'

diffmut <- function(mae) {

  # Find which assays contain mutation and expression level
  mut_assay <- grep("Mutation", names(mae))
  exp_assay <- grep("ExpressionLevel", names(mae))

  #Get names of samples with high and low expression, convert them to common format
  high <- TCGAbarcode(colnames(mae[[exp_assay]])[mae[[exp_assay]] == "high"])
  low  <- TCGAbarcode(colnames(mae[[exp_assay]])[mae[[exp_assay]] == "low"])

  #Find the patients in mutation assay corresponding to above
  overlap_high <- TCGAbarcode(colnames(mae[[mut_assay]])) %in% high
  overlap_high <- colnames(mae[[mut_assay]])[overlap_high]

  overlap_low <- TCGAbarcode(colnames(mae[[mut_assay]])) %in% low
  overlap_low <- colnames(mae[[mut_assay]])[overlap_low]
  all_overlaps <- c(overlap_high, overlap_low)

  #Get the mut data as a (large) sparse matrix, and extract the hugo symbols for use later. This is a huge assay.
  #Rows represent every possible site, and cells are either NA or Hugo Symbols instead of 1/0.
  #Next step reduces it to smaller matrix.
  mut <- mae[[mut_assay]]
  gene <- unique(as.data.frame(mut@assays)$Hugo_Symbol)
  mut <- mut[, c(overlap_high, overlap_low)]
  mut <- sparseAssay(mut)


  #Reduce mut assay into by-gene matrix instead of by-site
  mut_mat <- matrix(0, nrow= length(gene), ncol = length(all_overlaps),
                 dimnames = list(gene, all_overlaps))

  for(i in 1:length(all_overlaps)){ #Maybe data.table is faster? This takes 0.49 seconds
    sample <- all_overlaps[i]
    mut_genes <- mut[!is.na(mut[, sample]), sample]
    mut_mat[mut_genes, sample] <- 1
  }

  #calculates most mutated genes in the overlap patients
  m1 <- data.table(Gene = rownames(mut_mat),
                   N.high = rowSums(mut_mat[, overlap_high]),
                   N.hightotal = length(overlap_high),
                   N.low = rowSums(mut_mat[, overlap_low]),
                   N.lowtotal = length(overlap_low))

  #removes genes with no mutations or just 1 mutation in either group
  m2 <- m1[N.high + N.low != 0 & N.high + N.low != 1]

  #creates matrix with 4 rows
  xx = with(m2, matrix(c(N.hightotal - N.high, N.high, N.lowtotal - N.low, N.low), 4, byrow = TRUE))

  m2[, p.value := apply(xx, 2, function(x) {
    oopts <- options(warn = -1)
    on.exit(oopts)
    (chisq.test(matrix(x, 2))$p.value)
  })]

  m2[, FDR := p.adjust(p.value, method = "BH")]
  setkey(m2, p.value)
  m2
}
