#' plotlymut function
#'
#' Draws plotly plot of mutations
#'
#' @param mae MultiAssayExperiment object. Must contain assays Mutation (from curatedTCGA) and ExpressionLevel (from function rnasubset())
#' @param genemut data frame generated from diffmut function
#' @param gene character(1) gene of interest
#' @param num_genes numeric(1) number of genes to plot
#'
#' @import data.table
#' @import ggplot2
#' @import plotly
#'
#' @return interactive oncoprint-like plot for top differentially mutated genes
#'
#' @examples
#' mae <- curatedTCGAData("SKCM", c("RNASeq2GeneNorm", "Mutation", "GISTIC", "RPPAArray") , FALSE)
#' gene <- "SOX10"
#' sox10.pat <- rnasubset(mae, gene, 10)
#' sox10.mut <- diffmut(sox10.pat)
#' plotlymut(sox10.pat, mut, sox10.mut, gene)
#'
#' @export
#'
plotlymut <- function(mae, genemut, gene, num_genes) {

  #Find which assays contain expression levels and mutation data
  exp_assay <- grep("Expression", names(mae))
  mut_assay <- grep("Mutation", names(mae))

  #Get names of samples with high and low expression, convert them to common format
  high <- TCGAbarcode(colnames(mae[[exp_assay]])[mae[[exp_assay]] == "high"])
  low  <- TCGAbarcode(colnames(mae[[exp_assay]])[mae[[exp_assay]] == "low"])

  #Find the patients in mutation assay corresponding to above
  overlap_high <- TCGAbarcode(colnames(mae[[mut_assay]])) %in% high
  overlap_high <- colnames(mae[[mut_assay]])[overlap_high]

  overlap_low <- TCGAbarcode(colnames(mae[[mut_assay]])) %in% low
  overlap_low <- colnames(mae[[mut_assay]])[overlap_low]
  all_overlaps <- c(overlap_high, overlap_low)

  #Choose which genes to plot (num_genes)
  glist <- c(genemut$Gene[1:num_genes], gene)
  glist <- glist[!duplicated(glist)]

  #Get the mut data as a (large) sparse matrix, and extract the hugo symbols for use later. This is a huge assay.
  #Rows represent every possible site, and cells are either NA or Hugo Symbols instead of 1/0.
  #Next step reduces it to smaller matrix.
  mut <- mae[[mut_assay]]
  hugo <- unique(as.data.frame(mut@assays)$Hugo_Symbol)
  mut <- mut[, c(overlap_high, overlap_low)]
  mut <- sparseAssay(mut)

  #Reduce mut assay into by-gene matrix instead of by-site
  mut_mat <- matrix(0, nrow= length(hugo), ncol = length(all_overlaps),
                    dimnames = list(hugo, all_overlaps))

  for(i in 1:length(all_overlaps)){ #Maybe data.table is faster? This takes 0.49 seconds
    sample <- all_overlaps[i]
    mut_genes <- mut[!is.na(mut[, sample]), sample]
    mut_mat[mut_genes, sample] <- 1
  }

  #limit mut matrix to genes chosen above
  mut_mat <- mut_mat[glist,]

  #getting mutation table for high group
  mut1 <- mut_mat[glist, overlap_high]
  mut1 <- mut1[!(is.na(rowSums(mut1[, colnames(mut1)]))), ]

  #getting mutation table for low group
  mut2 <- mut_mat[glist, overlap_low]
  mut2 <- mut2[!(is.na(rowSums(mut2[, colnames(mut2)]))), ]

  #gets order of patients
  mut1.melt.high <- melt(mut1)
  gg1.high <- dcast(mut1.melt.high, Var2 ~ Var1)
  gene.order.high <- names(sort(apply(gg1.high[,colnames(gg1.high)[-1]], 2, sum), decreasing=T))
  mut1.melt.high$level <- "high"

  mut2.melt.low <- melt(mut2)
  gg1.low <- dcast(mut2.melt.low, Var2 ~ Var1)
  gene.order.low <- names(sort(apply(gg1.low[,colnames(gg1.low)[-1]], 2, sum), decreasing=T))
  mut2.melt.low$level <- "low"

  #combine mutation groups order patients according to sorting from individual dcasts
  mut1.melt <- rbind(mut1.melt.high, mut2.melt.low)
  names(mut1.melt) <- c("Gene", "barcode", "mutation", "level")

  mut1.melt$barcode <- factor(mut1.melt$barcode,
                                          levels = c(as.character(rev(gg1.high$Var2)),
                                                     as.character(rev(gg1.low$Var2))))
  mut1.melt$Gene <- factor(mut1.melt$Gene, levels = rev(gene.order.high))

  gg3 <- ggplot(mut1.melt, aes(x = barcode, y = Gene, fill = as.factor(mutation))) +
    geom_tile(colour = "white") +
    labs(x = "Patient", y = "Gene") +
    theme(title = element_text(size = 16),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    scale_fill_brewer(type = "qual",
                      name = "Legend",
                      palette = 6,
                      labels = c("Wild-type", "Mutated")) +
    ggtitle(paste(gene)) +
    facet_grid(. ~ level, scales = "free", space = "free")
  gg3
  #ggplotly(gg3)
}
