#loading tcga data list
load("cancers/SKCM/multi.RData")
load("cancers/SKCM/rppa.RData")
load("cancers/SKCM/cnasnp.RData")
library(devtools)
library(data.table)
library(TCGAbrowser)
gene.name <- fread("140428_TCGA_Entrez_Gene.txt")
cp.name <- fread("all_thresholded.by_genes.txt")
cp.name2 <- fread("ALANG_p_TCGA_180_SNP_1N_GenomeWideSNP_6_A01_895958.nocnv_hg19.seg.txt")
cp3 <- fread("c8a020ca-103f-4532-ac2b-4874b2e44b02/ALANG_p_TCGA_180_SNP_1N_GenomeWideSNP_6_A04_895912.nocnv_grch38.seg.txt")

#look for sysdata.rda
#writing R extensions
#create tempfile()
#dir.create(tempfile())

#cleaning up patient data
pat <- combi[[1]]
setkey(pat, bcr_patient_barcode, name)
#pat[, TNM := toupper(paste(pat$pathologyTstage, pat$pathologyNstage, pat$pathologyMstage))]
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
pat

#start test analysis
gene <- "SOX10"
percent <- 10
sox10 <- rnasubset(pat, combi[[2]], gene, percent)

#could use cut to define patient groups
cut(sox10$level, quantile(sox10$level, probs=c(0, 0.25, 0.75, 1)) ,labels = c("low", "medium", "high"))
.Machine$double.xmin

library(survival)
library(survminer)

#analysis for Tobias
INHBA
CYR61
ANGPTL4
FABP7

pat5 <- pat
#5 year
pat5$vitalstatus <- replace(pat5$vitalstatus, pat5$TCGA_year > 5, 0)
pat5$TCGA_year <- replace(pat5$TCGA_year, pat5$TCGA_year > 5, 5)

#3 year
pat3 <- pat
pat3$vitalstatus <- replace(pat3$vitalstatus, pat3$TCGA_year > 3, 0)
pat3$TCGA_year <- replace(pat3$TCGA_year, pat3$TCGA_year > 3, 3)

gene <- "INHBA"
inhba.new <- rnasubset(pat5, combi[[2]], gene, 11)
genesurv(inhba.new, gene)

gene <- "CYR61"
cyr61.new <- rnasubset(pat5, combi[[2]], gene, 50)
genesurv(cyr61.new, gene)

mcox <- coxph(Surv(TCGA_year, vitalstatus) ~(gene2 + yearstobirth + gender + pathologicstage), data = cyr61.new)
write.table(tidy(mcox), "CYR61_22_5year_mcox.txt", row.names=F, sep="\t")


gene <- "ANGPTL4"
ang.new <- rnasubset(pat5, combi[[2]], gene, 25)
genesurv(ang.new, gene)

gene <- "FABP7"
fabp7.new <- rnasubset(pat5, combi[[2]], gene, 25)
genesurv(fabp7.new, gene)

#survival plot from survminer
survplot <- survfit(Surv(TCGA_year, vitalstatus) ~ gene2, data=sox10, subset=gene2 != "middle")
half <- summary(survplot)$table[,"median"]
n <- summary(survplot)$table[,"records"]
res <- ggsurvplot(survplot,
                  pval=F,
                  pval.coord = c(1,0.1),
                  conf.int = F,
                  risk.table = T,
                  risk.table.y.text = F,
                  break.time.by = 1,
                  palette = brewer.pal(9, "Set1"),
                  legend = c(0.8, 0.8),
                  legend.title = "Gene level",
                  legend.labs = c(sprintf("%s high, n=%s\n Median OS %s years\n", gene, n[1], round(half[1],2)) ,
                                  sprintf("%s low, n=%s\n Median OS %s years\n", gene, n[2], round(half[2],2))),
                  xlab = "Time Since Biopsy (Years)")
res$table <- res$table + theme(axis.line = element_blank())
print(res)


#survival plot from highcharts
library(highcharter)
names(survplot$strata) <- c(sprintf("%s low", gene), sprintf("%s high", gene))
tooltip_table()

hchart(survplot, ranges=T) %>%
  hc_exporting(enabled=T) %>%
  hc_legend(align = "right", verticalAlign = "top", layout = "vertical", x = 0, y = 200) %>%
  hc_tooltip(shared=T, split=T, valueDecimals=2)

#survival plot from ggally and ggplot2
library(GGally)
library(plotly)
ggplotly(ggsurv(survplot, CI = TRUE) +
           geom_text(data=sox10, aes(label=name, x=TCGA_year, y=)) +
           theme_bw())

#expression plot with ggplot2
if( percent == 50) {
    cbPalette <- c( "#0072B2", "#D55E00")
  } else {
    cbPalette <- c( "#0072B2", "#999999","#D55E00")
  }
pat.d1.gene <- sox10
setkey(pat.d1.gene, level)
pat.d1.gene[, name := factor(name, levels=name)]
pat.d1.gene[, gene2 := factor(gene2, levels = c("high", "middle", "low"))]
g1 <- ggplot(data = pat.d1.gene, aes(x=name, y = level, colour = gene2)) +
  geom_point() +
  scale_colour_manual(values = cbPalette, labels = levels(pat.d1.gene$gene2), guide = guide_legend(title = NULL)) +
  xlab("Patient") +
  ylab("RSEM normalized read") +
  theme(axis.ticks=element_blank(), axis.text.x=element_blank(), panel.background=element_blank())
g1

#expression plot with high charts
setkey(sox10, level)
sox10[, name := factor(name, levels=name)]
hchart(sox10, "scatter", x = name, y = level, color = gene2)

#expression plot with plotly
library(plotly)
plot_ly(data = sox10, x = ~name, y = ~level, mode = "markers", color = ~gene2)

#mutation plot with plotly
#setkey(pat.gene, bcr_patient_barcode)
m1 <- combi[[4]]
setkey(pat.d1.gene, gene2, bcr_patient_barcode)
m1.high <- data.table(Gene = m1$Gene, N = rowSums(m1[, intersect(pat.d1.gene["high", bcr_patient_barcode], colnames(m1)), with =F]))
setkey(m1.high, N)
glist <- c(tail(m1.high[, Gene], 15), gene)
glist <- glist[!duplicated(glist)]
setkey(m1, Gene)
mut1 <- m1[glist,c("Gene", intersect(pat.d1.gene["high", bcr_patient_barcode], colnames(m1))), with=F]
mut1 <- mut1[!(is.na(rowSums(mut1[, colnames(mut1)[-1], with=F])))]
gg1 <- dcast.data.table((melt(mut1, variable.name="bcr_patient_barcode", id.vars="Gene")), bcr_patient_barcode ~ Gene)
setkeyv(gg1, names(sort(colSums(gg1[,colnames(gg1)[-1], with=F]), decreasing=T)))

#using order to reorder
colnames(mut1)[-1][with(mut1, order(colnames(mut1)[-1]))]
colnames(gg1)[-1][with(gg1, order(colnames(gg1)[-1]))]

#using do.call to order
gg1.order <- names(sort(colSums(gg1[, colnames(gg1)[-1], with=F]), decreasing=T))
setcolorder(gg1, c("bcr_patient_barcode", gg1.order))

mut1.order <- names(sort(colSums(mut1[, colnames(mut1)[-1], with=F]), decreasing=T))
setcolorder(mut1, c("Gene", mut1.order))

setkey(mut1, Gene)
mut1[gg1.order]

setkey(gg1, bcr_patient_barcode)
gg1[mut1.order]

gg1.order2 <- do.call(order, c(gg1[, colnames(gg1)[-1], with=F], decreasing = T))
mut1.order2 <-do.call(order, c(mut1[, colnames(mut1)[-1], with=F], decreasing = T))

nosort.gg1 <- gg1$bcr_patient_barcode[gg1.order2]
sort.gg1 <- gg1$bcr_patient_barcode[gg1.order2]

mut1$Gene[mut1.order2]


#using hclust to sort instead
hclust(dist(mut1[, colnames(mut1)[-1], with=F]))$order
hclust(dist(gg1[, colnames(gg1)[-1], with=F]))$order



mut2 <- melt(gg1)
mut2$bcr_patient_barcode <- factor(mut2$bcr_patient_barcode, levels=rev(gg1$bcr_patient_barcode))
mut2$variable <- factor(mut2$variable, levels=names(sort(colSums(gg1[,colnames(gg1)[-1], with=F]), decreasing=F)))
gg3 <- ggplot(mut2,aes(bcr_patient_barcode,variable,fill=as.factor(value))) + geom_tile(colour=c("white")) + labs(x = "Patient", y="Gene") +
  theme(title=element_text(size=16), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=12), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.text=element_text(size=12)) +
  scale_fill_brewer(type="qual",name="Legend", palette=6, labels=c("Wild-type", "Mutated")) +
  ggtitle(paste(gene, "high"))
gg3

ggplotly(gg3)

#differential expression with limma voom

setkey(pat2, gene2)
deg <- DGEList(counts=rna[,pat2[!("middle"), name], with=F], genes=rna$Gene, group=pat2[!("middle"),gene2])
isexpr <- rowSums(cpm(deg)>1) >= (ncol(deg)/2) #only keeps genes with at least 1 count-per-million in at least half the samples
deg <- deg[isexpr,]
design <- model.matrix(~factor(pat2[!("middle"),gene2], levels=c("low", "high")))
v2 <- voom(deg, design, plot=F)
fit <- lmFit(v2, design)
fit2 <- eBayes(fit)
fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.05, lfc=1,sort="p")
rnadeg <- fit3

#using gage to find KEGG pathway
data(kegg.gs)
limma.fc <- rnadeg$logFC
names(limma.fc) <- lookup$entrez[match(rnadeg$genes, lookup$gene)]
gage.fc <- gage(limma.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])

gage.fc.sig <- data.frame(rbind(gage.fc$greater[sel,],gage.fc$less[sel.l,]))
gage.fc.sig$enrichment <- rep(c("greater", "less"), c(table(sel)[[2]], table(sel.l)[[2]]))

path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

require(pathview)
out.suffix <- "limma"
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:3], function(pid) {
  pathview(gene.data = limma.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix)
})

#using GVSA to find more pathways

library(GSVA)
library(GSEABase)
library(GSVAdata)
data(c2BroadSets)

d1.gsva <- as.matrix(rna[,pat2[!("middle"), name], with=F])
rownames(d1.gsva) <- lookup$entrez[match(rna$Gene, lookup$gene)]
testgsva <- gsva(d1.gsva, c2BroadSets, min.sz=5, max.sz=500, mx.diff=TRUE, verbose=FALSE, rnaseq=TRUE)

design <- model.matrix(~factor(pat2[!("middle"),gene2], levels=c("low", "high")))

colnames(design) <- c(paste0(gene, "_low"), paste0(gene, "_high"))
fit <- lmFit(testgsva$es.obs, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=2, number=Inf)
DEgeneSets <- topTable(fit, coef=2, number=Inf, p.value=0.05, adjust="BH")
res <- decideTests(fit, p.value=0.05)
summary(res)

setkey(pat2, gene2)
h1 <- testgsva$es.obs[match(rownames(DEgeneSets)[1:100], rownames(testgsva$es.obs)),]
h1.t <- t(apply(h1, 1, scale))
colnames(h1.t) <- colnames(h1)
df <- data.frame(as.character(pat2[!("middle"),gene2]))
colnames(df) <- gene
col1 <- list(Cell = c("high" = "#ca0020",
                      "low" = "#0571b0"))
names(col1) <- gene
top_ha <- HeatmapAnnotation(df = df, col = col1)

h1.heat <- Heatmap(h1.t, top_annotation = top_ha, name = "color scale",
                   show_column_names = T,
                   row_names_gp = gpar(fontsize=8),
                   column_names_gp = gpar(fontsize=6),
                   row_names_max_width = unit(4, "cm"))

draw(h1.heat, heatmap_legend_side = "left", annotation_legend_side = "left")
for(an in colnames(df)) {
  decorate_annotation(an, {
    # annotation names on the right
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
    # annotation names on the left
    #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
  })
}

#import RPPA data
p1 <- data.table(rppa)
p1.names <- gsub("(TCGA-.*?-.*?-.*?)-.*", "\\1", colnames(p1))
setnames(p1, p1.names)

setkey(pat, name)
p1.pat <- pat[intersect(p1.names, pat$name)]

setkey(sox10, name)
rppa.sox10 <- intersect(colnames(p1), sox10[sox10$gene2 != "middle", name])

p1.sox10 <- sox10[rppa.sox10]


p1.gene <- as.matrix(p1[,p1.sox10$name, with=F])
rownames(p1.gene) <- p1$CE.REF

setkey(p1.sox10, gene2)
df <- data.frame(as.character(p1.sox10[!("middle"),gene2]),
                 p1.sox10[!("middle"), exprs_rank])
colnames(df) <- c(paste0(gene, "_group"), paste0(gene, "_expression"))
cellcol <- c("#ca0020", "#0571b0")
names(cellcol) <- c(levels(p1.sox10$gene2)[1], levels(p1.sox10$gene2)[2])
col1 <- list(Cell = cellcol,
             expression = colorRamp2(range(p1.sox10$exprs_rank), c("white", "purple")))
names(col1) <- c(paste0(gene, "_group"), paste0(gene, "_expression"))
top_ha <- HeatmapAnnotation(df = df, col = col1)

x <- Heatmap(p1.gene, top_annotation = top_ha, name = "color scale",
             show_column_names = T,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 6),
             column_dend_reorder = as.numeric(df[, 1]))
draw(x)


rna.log <- data.table(Gene = rna$Gene, cpm(rna[,pat2[!("middle"), name], with = FALSE], log = TRUE, normalized.lib.sizes = FALSE, prior.count = 8))
setkey(rna.log, Gene)
dgenes <- deg$genes[1:100][!is.na(deg$genes[1:100])]
h1 <- as.matrix(rna.log[dgenes, colnames(rna.log)[-1], with = FALSE])
rownames(h1) <- dgenes
h1.t <- t(apply(h1, 1, scale))
colnames(h1.t) <- colnames(h1)
setkey(rna, Gene)
df <- data.frame(as.character(pat2[!("middle"),gene2]),
                 pat2[!("middle"), exprs_rank])
colnames(df) <- c(paste0(gene, "_group"), paste0(gene, "_expression"))
cellcol <- c("#ca0020", "#0571b0")
names(cellcol) <- c(levels(pat2$gene2)[1], levels(pat2$gene2)[2])
col1 <- list(Cell = cellcol,
             expression = colorRamp2(c(1, nrow(pat2)), c("white", "purple")))
names(col1) <- c(paste0(gene, "_group"), paste0(gene, "_expression"))
top_ha <- HeatmapAnnotation(df = df, col = col1)

x <- Heatmap(h1.t, top_annotation = top_ha, name = "color scale",
             show_column_names = T,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 6),
             column_dend_reorder = as.numeric(df[, 1]))
draw(x)
#differential RPPA
design <- model.matrix(~p1.pat$gene)
fit <- lmFit(p1.gene, design)
fit2 <- eBayes(fit)
summary(decideTests(fit))
fit3 <- topTable(fit2, coef=2)
fit3$Antibody <- pro.name[as.numeric(rownames(fit3))]
fit3 <- fit3[,c(7,1:6)]
head(fit3)


#playing with RPPA data
rppadeg <- reactive({
  p1.pat <- rppapat()
  p1.gene <- rppagene()
  design <- model.matrix(~p1.pat$gene)
  fit <- lmFit(p1.gene, design)
  fit2 <- eBayes(fit)
  fit3 <- topTable(fit2, coef=2, n=Inf, adjust.method="BH", p.value=0.5, sort="p")
  fit3$Antibody <- pro.name[as.numeric(rownames(fit3))]
  fit3 <- fit3[,c(7,1:6)]
  fit3
})

#GoSeq package


#remove warnings
function() {
oopts <- options(warn = -1)
on.exit(options = oopts)
foo()
}

#for estimating copynumber
v <- c(1,0,1)

tabulate(v + 2, nbins=3)

v1 <- c(-1, 1, 1)

tabulate(v + 2, nbins=3)

cpmelt <- melt(cp, id.vars="Gene", variable.name = "name")
?xtabs
xtabs( ~ value + Gene,data=cpmelt)



#copy numner analysis
disjoin( withrevmap = T)
# revmap will make hits object to map back to original ranges

#test copy number analysis using disjoin granges from segmented files.
v1 <- "C:/Users/Phil/Documents/Phil/Post-doc/NGS/melarray/MelArray_Michael/ptarget1/10B_10A.vcf.gz"

#TCGAbiolinks copy number

library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-SKCM", data.category = "Copy Number Variation",
                    data.type = "Copy Number Segment")

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Clinical")
GDCdownload(query)

clinical.patient <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical.followup <- GDCprepare_clinic(query, clinical.info = "follow_up")
clinical.nte <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")

load("cancers/BRCA/multi.RData")

brca.rna <- combi[[2]]
brca.names <- gsub("(TCGA-.*?-.*?)-.*", "\\1", colnames(brca.rna)[-1])
brca.pat <- data.table(clinical.patient)

brca1 <- data.table(bcr_patient_barcode = brca.names, name = colnames(brca.rna)[-1])

brca.pat2 <- merge(brca.pat, brca1, by = "bcr_patient_barcode")

brca.pat2[ , c("ER", "PR", "HER2") := list(replace(breast_carcinoma_estrogen_receptor_status, breast_carcinoma_estrogen_receptor_status == "" | breast_carcinoma_estrogen_receptor_status == "Indeterminate", "Negative"),
                                           replace(breast_carcinoma_progesterone_receptor_status, breast_carcinoma_progesterone_receptor_status == "" | breast_carcinoma_progesterone_receptor_status == "Indeterminate", "Negative"),
                                           replace(lab_proc_her2_neu_immunohistochemistry_receptor_status, lab_proc_her2_neu_immunohistochemistry_receptor_status == "" | lab_proc_her2_neu_immunohistochemistry_receptor_status == "Equivocal" | lab_proc_her2_neu_immunohistochemistry_receptor_status == "Indeterminate", "Negative"))]
brca.pat2[, subtype := interaction(brca.pat2[, .(ER, PR, HER2)], sep= "_")]
brca.pat2[, subtype := factor(subtype, levels = unique(subtype))]
brca.pat2[, .N, by=subtype]

write.table(brca.pat2[, .N, by=subtype], "BRCA_subtype.txt", sep="\t", row.names = F)

gene <- "RNASE2"
brca.rna <- rnasubset(brca.pat2, combi[[2]], gene, 50)
brca.deg <- rnadeg(brca.rna, combi[[2]])
rnaheat(brca.rna, combi[[2]], brca.deg, gene)

pat2 <- brca.rna
rna <- combi[[2]]
deg <- brca.deg
setkey(pat2, gene2)
rna.log <- data.table(Gene = rna$Gene, cpm(rna[,pat2[!("middle"), name], with = FALSE], log = TRUE, normalized.lib.sizes = FALSE, prior.count = 8))
setkey(rna.log, Gene)
dgenes <- deg$genes[1:100][!is.na(deg$genes[1:100])]
h1 <- as.matrix(rna.log[dgenes, colnames(rna.log)[-1], with = FALSE])
rownames(h1) <- dgenes
h1.t <- t(apply(h1, 1, scale))
colnames(h1.t) <- colnames(h1)
setkey(rna, Gene)
df <- data.frame(as.character(pat2[!("middle"),gene2]),
                 pat2[!("middle"), exprs_rank],
                 pat2[!("middle"), subtype])
colnames(df) <- c(paste0(gene, "_group"), paste0(gene, "_expression"), "Subtype")
scols <- brewer.pal(8, "Set1")
names(scols) <- unique(brca.pat2$subtype)

col1 <- list(Cell = c("high" = "#ca0020",
                      "low" = "#0571b0"),
             expression = colorRamp2(c(1, nrow(pat2)), c("white", "purple")),
             Subtype = scols)
names(col1) <- c(paste0(gene, "_group"), paste0(gene, "_expression"), "Subtype")
top_ha <- HeatmapAnnotation(df = df, col = col1)

x <- Heatmap(h1.t, top_annotation = top_ha, name = "color scale",
             show_column_names = T,
             row_names_gp = gpar(fontsize = 8),
             column_names_gp = gpar(fontsize = 6),
             column_dend_reorder = as.numeric(df[, 1]))
draw(x)
for(an in colnames(df)) {
  decorate_annotation(an, {
    # annotation names on the right
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
    # annotation names on the left
    #grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right")
  })
}


options(shiny.reactlog=TRUE)

mdm2.rna <- rnasubset(pat, combi[[2]], "MDM2", 25)
mdm2.mut <- diffmut(mdm2.rna, combi[[4]])
mdm2.gg <- plotlymut(mdm2.rna, combi[[4]], mdm2.mut, "MDM2")

ggsave(filename="SKCM_MDM2_25_mutation.pdf", mdm2.gg, device=cairo_pdf, width = 12, height = 8, units = "in", dpi = 300)
genes <- c(mdm2.mut$Gene[1:20], "MDM2")
setkey(mdm2.mut, Gene)
write.table(mdm2.mut, "SKCM_MDM2_25_mutation.txt", sep="\t", row.names=F)

mdm4.rna <- rnasubset(pat, combi[[2]], "MDM4", 25)
mdm4.mut <- diffmut(mdm4.rna, combi[[4]])
mdm4.gg <- plotlymut(mdm4.rna, combi[[4]], mdm4.mut, "MDM4")

ggsave(filename="SKCM_MDM4_25_mutation.pdf", mdm4.gg, device=cairo_pdf, width = 12, height = 8, units = "in", dpi = 300)
genes <- c(mdm4.mut$Gene[1:20], "MDM4")
setkey(mdm4.mut, Gene)
write.table(mdm4.mut, "SKCM_MDM4_25_mutation.txt", sep="\t", row.names=F)

tp53.rna <- rnasubset(pat, combi[[2]], "TP53", 25)
tp53.mut <- diffmut(tp53.rna, combi[[4]])
tp53.gg <- plotlymut(tp53.rna, combi[[4]], tp53.mut, "TP53")

ggsave(filename="SKCM_TP53_25_mutation.pdf", tp53.gg, device=cairo_pdf, width = 12, height = 8, units = "in", dpi = 300)
genes <- c(tp53.mut$Gene[1:20], "TP53")
setkey(tp53.mut, Gene)
write.table(tp53.mut, "SKCM_TP53_25_mutation.txt", sep="\t", row.names=F)

##methylation test

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

load("cancers/SKCM/met.RData")

met1 <- data.table(REF = rownames(met), met)

z2 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

methylannot <- data.table("REF" = rownames(z2), "Island" = z2@listData$Relation_to_Island, "Island_Name" = z2@listData$Islands_Name, "Gene" = z2@listData$UCSC_RefGene_Name)
setkey(methylannot, "Island")

methylannot1 <- methylannot["Island"]
methyltest4 <- merge(met1, methylannot1, by="REF")
methyltest4[, "Hugo_Symbol" := tstrsplit(methyltest4$Gene, ";")[1]]
methyltest5a <- na.omit(methyltest4)
methyltest5b <- na.omit(methyltest4, invert=T)
setkey(methyltest5b, REF)

methyltest6 <- methyltest4[, lapply(.SD, function(x)median(x, na.rm=T)), .SDcols = !c("REF", "Gene", "Hugo_Symbol", "Island"), by="Island_Name"]
methyltest7 <- na.omit(methyltest6)
methyltest8 <- merge(methyltest7, unique(methyltest4[, .(Island_Name, Hugo_Symbol)]), by="Island_Name")
setnames(methyltest8, colnames(methyltest8), sub("(TCGA-.*?-.*?-.*?)-.*", "\\1", colnames(methyltest8)))
methyltest8

test1 <- as.numeric(methyltest8[450, setdiff(colnames(methyltest8), c("REF", "Hugo_Symbol", "Island", "Island_Name")), with = FALSE])
plot(density(test1))
plot(density(asin(sqrt(test1))))


pat <- fread("/data/Phil/TCGA_SKCM_141203/shinyapp/megatev2/141203_TCGA_patient_shiny_SKCM.txt")
methylshiny <- methyltest8[, c("GeneSymbol", "Island_Name", pat$name), with=F]
setkey(methylshiny, GeneSymbol)
d1 <- fread("/data/Phil/TCGA_SKCM_141203/shinyapp/megatev2/141203_TCGA_RNAseq_norm_shiny_SKCM.txt")
setkey(d1, Gene)

gene <- "SOX9"
apply(methylshiny[gene, colnames(methylshiny)[-1:-2], with=F], 1, function(x) {
  plot(x ~ as.numeric(d1[gene, colnames(d1)[-1], with=F]))
})
methylshiny2 <- methylshiny[!(is.na(methylshiny$GeneSymbol))]

library(TCGAbiolinks)

