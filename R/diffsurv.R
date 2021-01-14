#' diffsurv function
#'
#' Returns a data table of median survival time and 95 confidence intervals of the two cohorts
#'
#' @param mae MultiAssayExperiment object with Cohort assay generated from rnasubset or mutsubset.
#' @param gene character(1) gene of interest
#'
#' @import survival
#' @importFrom survminer surv_pvalue
#'
#' @return data frame of survival statistics
#'
#' @examples
#' library(curatedTCGAData)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' diffsurv(lusc_t.egfr, "EGFR")
#'
#' @export
#'
diffsurv <- function(mae, gene) {
  #check if Cohort matrix has been made for MultiAssayExperiment
  stopifnot(any(grepl("Cohort", names(mae))))

  #only take clinical data for samples in Cohort
  mae <- intersectColumns(mae[, , "Cohort"])

  d1 <- longFormat(mae["level", , "Cohort"])
  if("vital_status" %in% colnames(colData(mae))) {
    cd1 <- merge(as.data.frame(d1), as.data.frame(colData(mae)[, c("patientID", "days_to_death", "vital_status", "days_to_last_followup")]), by.x = "primary", by.y = "patientID")
    cd1$years <- ifelse(cd1$vital_status == 1,
                        round(cd1$days_to_death/365.25,2),
                        round(cd1$days_to_last_followup/365.25, 2))
    cd2 <- cd1[cd1$value != "medium",]
    cd2$gene2 <- factor(cd2$value, levels = c("high", "low"))
    survplot <- survfit(Surv(years, vital_status) ~ gene2, data = cd2)

  } else {
    cd1 <- merge(as.data.frame(d1), as.data.frame(colData(mae)[, c("patientID", "days_to_death.x", "vital_status.x", "days_to_last_followup.x")]), by.x = "primary", by.y = "patientID")
    cd1$years <- ifelse(cd1$vital_status.x == 1,
                        round(cd1$days_to_death.x/365.25,2),
                        round(cd1$days_to_last_followup.x/365.25, 2))
    cd2 <- cd1[cd1$value != "medium",]
    cd2$gene2 <- factor(cd2$value, levels = c("high", "low"))
    survplot <- survfit(Surv(years, vital_status.x) ~ gene2, data = cd2)
  }

  t1 <- summary(survplot)$table
  data.table(Cohort = c(paste(gene, "high"),
                        paste(gene, "low")),
             n = c(t1[grep("high", rownames(t1)), "records"],
                   t1[grep("low", rownames(t1)), "records"]),
             events = c(t1[grep("high", rownames(t1)), "events"],
                        t1[grep("low", rownames(t1)), "events"]),
             median = c(t1[grep("high", rownames(t1)), "median"],
                        t1[grep("low", rownames(t1)), "median"]),
             LCI_95 = c(t1[grep("high", rownames(t1)), "0.95LCL"],
                        t1[grep("low", rownames(t1)), "0.95LCL"]),
             UCI_95  = c(t1[grep("high", rownames(t1)), "0.95UCL"],
                         t1[grep("low", rownames(t1)), "0.95UCL"]),
             pvalue = surv_pvalue(survplot)$pval)
  }

