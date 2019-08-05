#' methylprep function
#'
#' Converts the methylation 450 SummarizedExperiment into a matrix of gene associated CpG Islands and samples
#' Uses CpG island definitions from `r Biocpkg("IlluminaHumanMethylation450kanno.ilmn12.hg19")`
#'
#' @param mae MultiAssayExperiment object containing Methylation assay and Cohort assay from rnasubset
#'
#' @import data.table
#' @import MultiAssayExperiment
#' @importFrom minfi getAnnotation
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#'
#' @return Returns the MultiAssayExperiment object with a new assay called SimpleMethyl contains a matrix of gene associated CpG islands and samples
#'
#'
#' @examples
#' #using data from the cureatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("RNASeq2GeneNorm", "Methylation_methyl450"), FALSE)
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' #prepare methylation data for simple differential methylation analysis
#' lusc_t <- methylprep(lusc_t)
#' #subset samples by top and bottom 10% of EGFR expression
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' #look for differentially methylated islands
#' egfr_meth_deg <- diffmeth(lusc_t.egfr)
#'
#' @export
#'
methylprep <- function(mae) {
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
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

  #Extract methylation data for all genes, then
  #limit to loci with Island annotation from Illumina
  meth_assay_num <- grep("Methylation", names(mae))
  meth <- mae[[meth_assay_num]]
  meth <- meth[illumina$REF,]

  #Merge methylation data with annotation data
  merged <- as.data.table(assay(meth))
  merged[, "REF" := illumina$REF]
  merged <- merge(merged, illumina, by = "REF")

  #Calculate the median methylation value per island
  methylisland <- merged[, lapply(.SD, function(x) median(x, na.rm=T)),
                        .SDcols = !c("REF", "Chr", "Loc", "Island", "Gene", "UCSC_RefGene_Group"), by = Island_Name]
  methylisland <- na.omit(methylisland)

  methylisland2 <- as.matrix(methylisland[, setdiff(colnames(methylisland), "Island_Name"), with = FALSE])
  rownames(methylisland2) <- methylisland$Island_Name
  #Append expression level matrix to original MAE object
  mae2 <- c(mae, SimpleMethyl = methylisland2, mapFrom = meth_assay_num)

  #return MAE object with new assay indicating expression level
  mae2

}
