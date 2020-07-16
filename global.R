#library(gplots)
library(plotly)
library(data.table)
library(shiny)
library(survival)
library(limma)
library(edgeR)
library(googleVis)
library(gage)
#library(plyr)
library(reshape2)
#library(STRINGdb)
library(ggplot2)
library(grid)
library(shinydashboard)
#library(rCharts)
#library(d3heatmap)
library(magrittr)
library(RColorBrewer)
library(DT)
#library(heatmap3)
library(TCGAbrowser)
library(ComplexHeatmap)
library(GSVAdata)

data(c2BroadSets)
data(kegg.gs)
cancers <- fread("nexus_tcga_cancer.txt")
setkey(cancers, Cancer_name)

load("./data/cancers/SKCM/multi.RData")
load("./data/cancers/SKCM/rppa.RData")

pat <- combi[[1]]
setkey(pat, bcr_patient_barcode, name)
pat$vitalstatus <- as.numeric(pat$vitalstatus)
pat$yearstobirth <- as.numeric(pat$yearstobirth)
pat$daystosubmittedspecimendx <- as.numeric(pat$daystosubmittedspecimendx)
pat <- pat[!(is.na(pat$name))]
pat <- pat[!(is.na(pat$days))]
pat[, years := round(days/365.25, 2)]
pat[, TCGA_day := days - daystosubmittedspecimendx]
pat[, TCGA_year := round(TCGA_day/365.25,2)]
pat <- pat[!is.na(pat$TCGA_day)]
pat$gender <- factor(pat$gender, levels=c("male", "female"))
setkey(pat, name)
pat$age_group <- cut(pat$yearstobirth + pat$daystosubmittedspecimendx/365.25, br=seq(from = 9, to = 99, by = 10), labels = c("10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99"))
pat$age_group <- factor(pat$age_group, levels=(sort(unique(pat$age_group))))


d1 <- combi[[2]]
setkey(d1, Gene)


lookup <- fread("140428_TCGA_Entrez_Gene.txt", sep="\t", header=T)
setkey(lookup, entrez)
gene.name <- d1$Gene

cp1 <- combi[[3]]
setkey(cp1, Gene)

m1 <- combi[[4]]
setkey(m1, Gene)

p1 <- data.table(rppa)
p1.names <- gsub("(TCGA-.*?-.*?-.*?)-.*", "\\1", colnames(p1))
setnames(p1, p1.names)

gene.name2 <- unique(m1$Gene)

rm(combi)

cellinfo <- function(x){
  if(is.null(x)) return(NULL)
  paste(x$callup)
}

cellinfo1 <- function(x){
  if(is.null(x)) return(NULL)
  paste(x$name, x$level)
}

