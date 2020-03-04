#' plotclin function
#'
#' Plots box or bar plots of clinical features of the two cohorts
#'
#' @param mae MultiAssayExperiment object containing colData and Cohort assay from rnasubset
#' @param var vector of columns from colData for summarize, default is years_to_birth, gender, histological_type and pathologic stage
#' @import data.table
#' @import ggplot2
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
#' #Plot clinical differences between the high and low group
#' plotclin(lusc_t.egfr)
#' @export
#'

plotclin <- function(mae, var) {
  #check if Cohort matrix has been made for MultiAssayExperiment
  stopifnot(any(grepl("Cohort", names(mae))))
  if(var %in% colnames(colData(lusc.t3)) != TRUE)
  #only take clinical data for samples in Cohort
  mae <- intersectColumns(mae[, , "Cohort"])

  d1 <- longFormat(mae["level", , "Cohort"])
  d1$level <- as.character(d1$value)
  d2 <- longFormat(mae["counts", , "Cohort"])
  d2$counts <- as.numeric(levels(d2$value))
  cd1 <- merge(as.data.frame(d1), as.data.frame(colData(mae)[, c("patientID", var)]), by.x = "primary", by.y = "patientID")
  cd2 <- merge(as.data.frame(d2), cd1, by = c("primary", "assay", "colname"))
  setDT(cd2)
  setkey(cd2, level)

  if(class(cd2[, get(var)]) == "numeric" | class(cd2[, get(var)]) == "integer") {
    if(abs(median(cd2["low", get(var)], na.rm = TRUE) - median(cd2["high", get(var)], na.rm = TRUE)) < 1000) {
      ggplot(cd2[c("low", "high")], aes_string(x = "level", y = var, fill = "level")) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = "jitter") +
        theme_bw() } else {
          ggplot(cd2[c("low", "high")], aes_string(x = "level", y = var, fill = "level")) +
            geom_boxplot(outlier.shape = NA) +
            geom_point(position = "jitter") +
            scale_y_log10() +
            theme_bw()
        }
  } else if(class(cd2[, get(var)]) == "character") {
    ggplot(cd2[c("low", "high")], aes_string(x = "level", fill = var)) +
      geom_bar() +
      scale_fill_brewer(type = "qual", palette = 3) +
      theme_bw()
  }
}
