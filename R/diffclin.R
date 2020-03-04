#' diffclin function
#'
#' Creates a table that summarizes the clinical features between the two cohorts
#'
#' @param mae MultiAssayExperiment object containing colData and Cohort assay from rnasubset
#' @param vars vector of columns from colData for summarize, default is years_to_birth, gender, histological_type and pathologic_stage
#' @import data.table
#' @import gtsummary
#'
#' @return gtsummary table of clinical features
#'
#' @examples
#' #using data from the cureatedTCGAdata set
#' library(curatedTCGAData)
#' library(TCGAutils)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm", "GISTIC_ThresholdedByGene", "RPPAArray", "Methylation_methyl450"), FALSE)
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' #subset LUSC data by expression of EGFR
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' #Look at clinical differences between the high and low group
#' diffclin(lusc_t.egfr)
#' @export
#'

diffclin <- function(mae, vars = c("years_to_birth", "gender", "histological_type", "pathologic_stage")) {
  #check if Cohort matrix has been made for MultiAssayExperiment
  stopifnot(any(grepl("Cohort", names(mae))))

  #only take clinical data for samples in Cohort
  mae <- intersectColumns(mae[, , "Cohort"])

  d1 <- longFormat(mae["level", , "Cohort"])
  d1$level <- as.character(d1$value)
  d2 <- longFormat(mae["counts", , "Cohort"])
  d2$counts <- as.numeric(levels(d2$value))
  if("vital_status" %in% colnames(colData(mae))) {
    cd1 <- merge(as.data.frame(d1), as.data.frame(colData(mae)[, c("patientID", "days_to_death", "vital_status", "days_to_last_followup", vars)]), by.x = "primary", by.y = "patientID")
    cd2 <- merge(as.data.frame(d2), cd1, by = c("primary", "assay", "colname"))


    #summarize clinical features using gtsummary
    setDT(cd2)
    setkey(cd2, level)

    tbl_summary(cd2[c("high", "low"), c("level", "counts", vars), with = FALSE], by = "level") %>%
      add_overall() %>%
      add_n() %>%
      add_p() %>%
      bold_labels()
  }
}

