#' rnamethdeg function
#'
#' Uses limma to find differentially methylated islands between the high and low groups
#'
#' @param mae MultiAssayExperiment object containing Methylation assay and Cohort assay from rnasubset
#'
#' @import data.table
#' @import edgeR
#' @import limma
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#'
#' @return data frame of significant differentially methylated islands between the two groups defined in rnasubset or mutsubset
#'
#' @examples
#' #using data from the cureatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTICT", "RPPAArray", "Methylation_methyl450"), FALSE)
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' egfr_meth_deg <- rnamethdeg(lusc_t.egfr)
#'
#' @export
#'
rnamethdeg <- function(mae) {
  stopifnot(any(grepl("Cohort", names(mae))))
  stopifnot(any(grepl("Methylation", names(mae))))

  #Get Illumina methylation island data
  illumina <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  illumina <- data.table("REF" = rownames(illumina),
                         "Chr" = illumina@listData$chr,
                         "Loc" = illumina@listData$pos,
                         "Island" = illumina@listData$Relation_to_Island,
                         "Island_Name" = illumina@listData$Islands_Name,
                         "Gene" = illumina@listData$UCSC_RefGene_Name,
                         "UCSC_RefGene_Group" = illumina@listData$UCSC_RefGene_Group
  )
  setkey(illumina, "Island")
  illumina <- illumina["Island"]

  #Get expression levels (which were appended by function rnasubset()), ID those with high / low expression
  level_assay_num <- grep("Cohort", names(mae))
  lvl <- mae[[level_assay_num]]
  good_levels <- names(lvl[1, lvl[1,] != "medium"])

  #Extract methylation data for all genes, then
  #subset out the "medium" expression group, and
  #limit to loci with Island annotation from Illumina
  meth_assay_num <- grep("Methylation", names(mae))
  meth <- mae[[meth_assay_num]]
  meth <- meth[illumina$REF, TCGAbarcode(colnames(meth)) %in% TCGAbarcode(good_levels)]

  #Merge methylation data with annotation data
  merged <- as.data.table(assay(meth))
  merged[, "REF" := illumina$REF]
  merged <- merge(merged, illumina, by = "REF")

  #Calculate the mean methylation value per island
  methyltest6 <- merged[, lapply(.SD, function(x) median(x, na.rm=T)),
                        .SDcols = !c("REF", "Chr", "Loc", "Island", "Gene", "UCSC_RefGene_Group"), by = Island_Name]
  methyltest7 <- na.omit(methyltest6)

  #Associate expression levels with methylation names
  levels <- lvl
  colnames(levels) <- TCGAbarcode(colnames(levels))
  levels <- levels[, TCGAbarcode(colnames(meth))]

  #convert median methylation data to matrix and create model matrix
  meth.mat <- as.matrix(methyltest7[, 2:ncol(methyltest7)])
  row.names(meth.mat) <- methyltest7$Island_Name
  design <- model.matrix(~ factor(levels["level", ]))

  #Perform DGE analysis on asin transformed data and return result
  fit <- lmFit(asin(meth.mat), design)
  fit2 <- eBayes(fit, trend = TRUE)
  fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.05, lfc = 0.1, sort="p")
  fit3

}
