#' plotsurv function
#'
#' Draws survminer survival plot for high and low gene populations
#'
#' @param mae MultiAssayExperiment object with Cohort assay generated from rnasubset or mutsubset.
#' @param gene character(1) gene of interest
#'
#' @import survival
#' @import survminer
#'
#' @return survminer styled Survival plot
#'
#' @examples
#' library(curatedTCGAData)
#' lusc <- curatedTCGAData("LUSC", c("Mutation", "RNASeq2GeneNorm"), FALSE)
#'
#' #split tumor and normal samples
#' lusc_tn <- splitAssays(lusc2, c("01", "11"))
#' #remake MultiAssayExperiment with only primary tumor samples
#' lusc_t <- lusc_tn[, , grep("^01", names(lusc_tn))]
#' lusc_t.egfr <- rnasubset(lusc_t, "EGFR", 10)
#' plotsurv(lusc_t.egfr, "EGFR")
#'
#' @export
#'
plotsurv <- function(mae, gene) {
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

  half <- summary(survplot)$table[,"median"]
  n <- summary(survplot)$table[,"records"]
  res <- ggsurvplot(survplot,
                    data = cd2,
                    pval=T,
                    pval.coord = c(1,0.1),
                    conf.int = T,
                    risk.table = T,
                    risk.table.y.text = F,
                    break.time.by = 1,
                    palette = c("#1a1a1a", "#d7191c"),
                    legend = c(0.8, 0.8),
                    legend.title = "Gene level",
                    legend.labs = c(sprintf("%s %s, n=%s\n Median OS %s years\n", gene, levels(cd2$gene2)[1], n[1], round(half[1],2)) ,
                                    sprintf("%s %s, n=%s\n Median OS %s years\n", gene, levels(cd2$gene2)[2], n[2], round(half[2],2))),
                    xlab = "Time Since Biopsy (Years)")
  res$table <- res$table + theme(axis.line = element_blank())
  print(res)
}
