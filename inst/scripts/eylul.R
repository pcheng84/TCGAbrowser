skcm <- readRDS("data/skcmMAEO.rds")
skcm

aa <- colnames(skcm)
aa[[1L]] <- colnames(skcm)[["RNASeq2GeneNorm"]][1:10]
skcm["SOX10", aa, "RNASeq2GeneNorm"]

m <- matrix(rpois(20 * 30, 0.2), 20)
tidy <- melt(m) %>%  tbl_df
tidy %>% group_by(Var2) %>%  mutate(order=order(value))


#analysis for Eylul
library(edgeR)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)

rna <- combi[[2]]
rna.log <- data.table(Gene = rna$Gene, cpm(rna[,colnames(rna)[-1], with=F], log=T, normalized.lib.sizes=F))
setkey(rna.log, Gene)
e1 <- c("MITF", "ZEB1", "ZEB2", "AXL", "WNT5A", "SOX10", "SOX9", "TGFB2", "BMP7")
h1 <- as.matrix(rna.log[e1, colnames(rna.log)[-1], with = FALSE])
rownames(h1) <- e1
h1.t <- t(apply(h1, 1, scale))
colnames(h1.t) <- colnames(h1)

pdf(file = "161108_eylul_genes.pdf", width = 14, height = 8)
dendrow <- as.dendrogram(hclust(dist(h1.t)))
dendrow <- set(dendrow, "branches_lwd", 1.5)

dendcol <- as.dendrogram(hclust(dist(t(h1.t))))
dendcol <- set(dendcol, "branches_lwd", 1.5)
Heatmap(h1.t, name = "color scale", cluster_rows = dendrow, cluster_columns = dendcol, column_dend_height = unit(2, "cm"), show_column_names = FALSE)
dev.off()

cor.test(as.numeric(rna.log["SMAD7", colnames(rna.log)[-1], with=F]),
    as.numeric(rna.log["ZEB1", colnames(rna.log)[-1], with=F]))


ey.cor <- melt(rna.log[c("MITF", "SOX10", "ZEB1", "WNT5A", "SOX9", "AXL", "SMAD7", "ZEB2", "BMP7", "TGFB2")])
ey.cor2 <- dcast.data.table(ey.cor, variable ~ Gene)

cor.test(ey.cor2$SMAD7, ey.cor2$ZEB1)
ggplot(data=ey.cor2, aes(x=SMAD7, y=ZEB1)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() +
  annotate("text", x = 1, y = 7.5, label = "p < 0.001") +
  annotate("text", x = 1, y = 8.5, label = "r = 0.53") +
  xlab("SMAD7 (logCPM)") +
  ylab("ZEB1 (logCPM)")
  #geom_text_repel(aes(label = variable))

ggsave("Eylul_SMAD7_ZEB1_correlation.pdf", width=6, height=6)

cor.test(ey.cor2$SMAD7, ey.cor2$CCND1)
ggplot(data=ey.cor2, aes(x=SMAD7, y=CCND1)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() +
  annotate("text", x = 1, y = 13, label = "p = 0.23") +
  annotate("text", x = 1, y = 13.5, label = "r = -0.05") +
  xlab("SMAD7 (logCPM)") +
  ylab("CCND1 (logCPM)")
#geom_text_repel(aes(label = variable))

ggsave("Eylul_SMAD7_CCND1_correlation.pdf", width=6, height=6)

cor.test(ey.cor2$SMAD7, ey.cor2$ZEB2)
ggplot(data=ey.cor2, aes(x=SMAD7, y=ZEB2)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() +
  annotate("text", x = 1, y = 9, label = "p < 0.001") +
  annotate("text", x = 1, y = 9.5, label = "r = 0.26") +
  xlab("SMAD7 (logCPM)") +
  ylab("ZEB2 (logCPM)")
#geom_text_repel(aes(label = variable))

ggsave("Eylul_SMAD7_MITF_correlation.pdf", width=6, height=6)

cor.test(ey.cor2$SMAD7, ey.cor2$MITF)
ggplot(data=ey.cor2, aes(x=SMAD7, y=MITF)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() +
  annotate("text", x = 6.5, y = 10, label = "p < 0.001") +
  annotate("text", x = 6.5, y = 10.5, label = "r = -0.19") +
  xlab("SMAD7 (logCPM)") +
  ylab("MITF (logCPM)")
#geom_text_repel(aes(label = variable))

ggsave("Eylul_SMAD7_MITF_correlation.pdf", width=6, height=6)

cor.test(ey.cor2$SMAD7, ey.cor2$SOX10)
ggplot(data=ey.cor2, aes(x=SMAD7, y=SOX10)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() +
  annotate("text", x = 1, y = 8, label = "p < 0.001") +
  annotate("text", x = 1, y = 8.5, label = "r = -0.24") +
  xlab("SMAD7 (logCPM)") +
  ylab("SOX10 (logCPM)")
#geom_text_repel(aes(label = variable))

ggsave("Eylul_SMAD7_SOX10_correlation.pdf", width=6, height=6)

cor.test(ey.cor2$SMAD7, ey.cor2$AXL)
ggplot(data=ey.cor2, aes(x=SMAD7, y=AXL)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
  theme_bw() +
  annotate("text", x = 1, y = 9, label = "p < 0.001") +
  annotate("text", x = 1, y = 9.5, label = "r = 0.28") +
  xlab("SMAD7 (logCPM)") +
  ylab("AXL (logCPM)")
#geom_text_repel(aes(label = variable))

ggsave("Eylul_SMAD7_AXL_correlation.pdf", width=6, height=6)



x1 <- "SMAD7"
y1 <- "AXL"

ecor <- function(x1, y1) {
  eycor <- with(ey.cor2, cor.test(get(x1), get(y1)))
  gg1 <- ggplot(data=ey.cor2, aes(x=get(x1), y=get(y1))) +
    geom_point() +
    geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95) +
    theme_bw() +
    annotate("text",
             x = with(ey.cor2, max(get(x1))) - 0.5,
             y = with(ey.cor2, max(get(y1))) - 1,
             label = sprintf("p <  %s", signif(eycor$p.value, 3))) +
    annotate("text",
             x = with(ey.cor2, max(get(x1))) - 0.5,
             y = with(ey.cor2, max(get(y1))) - 0.5,
             label = sprintf("r =  %s", signif(eycor$estimate, 3))) +
    xlab(sprintf("%s (logCPM)", x1)) +
    ylab(sprintf("%s (logCPM)", y1))
#geom_text_repel(aes(label = variable))
  ggsave(sprintf("Eylul_%s_%s_correlation.pdf", x1, y1), plot = gg1, width=6, height=6)
}

ecor("BMP7", "AXL")
ecor("BMP7", "MITF")
ecor("BMP7", "SMAD7")
ecor("BMP7", "SOX10")
ecor("BMP7", "SOX9")
ecor("BMP7", "TGFB2")
ecor("BMP7", "WNT5A")
ecor("BMP7", "ZEB1")
ecor("BMP7", "ZEB2")

ecor("TGFB2", "AXL")
ecor("TGFB2", "MITF")
ecor("TGFB2", "SMAD7")
ecor("TGFB2", "SOX10")
ecor("TGFB2", "SOX9")
ecor("TGFB2", "BMP7")
ecor("TGFB2", "WNT5A")
ecor("TGFB2", "ZEB1")
ecor("TGFB2", "ZEB2")

library(corrplot)
eycor2 <- cor(ey.cor2[, .(AXL, BMP7, MITF, SMAD7, SOX10, SOX9, TGFB2, WNT5A, ZEB1, ZEB2)])
corrplot(eycor2, order = "hclust")
