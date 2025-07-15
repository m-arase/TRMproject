library(Seurat)
library(Signac)
library(hdf5r)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)

# FigS4D
result <- read.csv(file = "RUNX2OEvolcano.csv")
result$log2FC <- ifelse(result$FC > 0, log2(result$FC), -log2(abs(result$FC)))
result$log2FC[is.nan(result$log2FC)] <- 0
result$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result$diffexpressed[result$log2FC > 1 & result$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result$diffexpressed[result$log2FC < -1 & result$pvalue < 0.05] <- "DOWN"

result$delabel <- NA
result$delabel[result$diffexpressed != "NO"] <- result$gene[result$diffexpressed != "NO"]
result$volcano <- -log10(result$pvalue)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)

pdf("RUNX2OE_volcano.pdf", width = 10, height = 10)
ggplot(data=result, aes(x=log2FC, y=volcano, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  scale_y_continuous(limits=c(0,8),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

result <- read.csv(file = "BHLHE40OEvolcano.csv")
result$log2FC <- ifelse(result$FC > 0, log2(result$FC), -log2(abs(result$FC)))
result$log2FC[is.nan(result$log2FC)] <- 0
result$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result$diffexpressed[result$log2FC > 1 & result$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result$diffexpressed[result$log2FC < -1 & result$pvalue < 0.05] <- "DOWN"

result$delabel <- NA
result$delabel[result$diffexpressed != "NO"] <- result$gene[result$diffexpressed != "NO"]
result$volcano <- -log10(result$pvalue)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)

pdf("BHLHE40OE_volcano.pdf", width = 10, height = 10)
ggplot(data=result, aes(x=log2FC, y=volcano, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  scale_y_continuous(limits=c(0,7),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

result <- read.csv(file = "BothOEvolcano.csv")
result$log2FC <- ifelse(result$FC > 0, log2(result$FC), -log2(abs(result$FC)))
result$log2FC[is.nan(result$log2FC)] <- 0
result$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result$diffexpressed[result$log2FC > 1 & result$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result$diffexpressed[result$log2FC < -1 & result$pvalue < 0.05] <- "DOWN"

result$delabel <- NA
result$delabel[result$diffexpressed != "NO"] <- result$gene[result$diffexpressed != "NO"]
result$volcano <- -log10(result$pvalue)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)

pdf("BothOE_volcano.pdf", width = 10, height = 10)
ggplot(data=result, aes(x=log2FC, y=volcano, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  scale_y_continuous(limits=c(0,9),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

# Multiome data from Figure5
# FigS4F
pdf("Coverageplot_IFNG.pdf", width = 6, height = 3)
CoveragePlot(
  object = Multiome_data,
  region = "IFNG",
  idents = idents.plot,
  extend.upstream = 40000,
  extend.downstream = 40000
)
dev.off()

pdf("Coverageplot_GZMB.pdf", width = 6, height = 3)
CoveragePlot(
  object = Multiome_data,
  region = "GZMB",
  idents = idents.plot,
  extend.upstream = 50000,
  extend.downstream = 50000
)
dev.off()
