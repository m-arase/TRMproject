library(Seurat)
library(Signac)
library(hdf5r)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
set.seed(1234)

# obj : CITE-seq data from Figure1
# FigS1A
pdf("CITE_UMAP_sample_2.pdf", width = 6, height = 6)
DimPlot(obj, reduction = 'umap.rpca', group.by = "sample", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()
write.csv(table(obj$sample, obj$annotationnew), file = "CITE_prop_sample.csv")

# FigS1B
# RNA
DefaultAssay(obj) <- "RNA"
result <- FindMarkers(obj, ident.1 = "TRM_2", group.by = "annotationnew")
result$gene <- rownames(result)

result$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result$diffexpressed[result$avg_log2FC > 1 & result$p_val_adj < 10^-50] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result$diffexpressed[result$avg_log2FC < -1 & result$p_val_adj < 10^-50] <- "DOWN"

result$delabel <- ifelse(result$diffexpressed != "NO", result$gene, NA)
result$volcano <- -log10(result$p_val_adj)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)

pdf("CITE_CD_TRM_2_DEG_volcano.pdf", width = 10, height = 10)
ggplot(data=result, aes(x=avg_log2FC, y=volcano, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel(size = 5) +
  scale_color_manual(values=c("blue", "black", "red")) +
  scale_y_continuous(limits=c(0,320),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

result2 <- result[result$volcano > 290, ]
pdf("CITE_CD_TRM_2_DEG_volcano2.pdf", width = 30, height = 10)
ggplot(data=result2, aes(x=avg_log2FC, y=volcano, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel(size = 5, max.overlaps = 30) +
  scale_color_manual(values=c("blue", "black", "red")) +
  ylab("-log10_p_val")
dev.off()

# ADT
DefaultAssay(obj) <- "ADT"
result <- FindMarkers(obj, ident.1 = "TRM_2", group.by = "annotationnew")
result$gene <- rownames(result)

result$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result$diffexpressed[result$avg_log2FC > 0.5 & result$p_val_adj < 10^-50] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result$diffexpressed[result$avg_log2FC < -0.5 & result$p_val_adj < 10^-50] <- "DOWN"

result$delabel <- ifelse(result$diffexpressed != "NO", result$gene, NA)
result$volcano <- -log10(result$p_val_adj)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)

pdf("CITE_CD_TRM_2_DEADT_volcano.pdf", width = 10, height = 10)
ggplot(data=result, aes(x=avg_log2FC, y=volcano, col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel(size = 5) +
  scale_color_manual(values=c("blue", "black", "red")) +
  scale_y_continuous(limits=c(0,320),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

# FigS1C
library(Seurat)
library(tidyverse)
library(escape)
set.seed(1234)

TRMgenelist <- read.csv(file = "/Volumes/BUFFALOHDD/CITEseqAnalysis/TRMpublic.csv", skip = 1)
gene.sets1 <- list(TRMgenelist[1:100,1])
names(gene.sets1) <- "gsea1"
trial.ssGSEA1 <- escape.matrix(obj@assays$RNA$counts, 
                               method = "ssGSEA", 
                               gene.sets = gene.sets1, 
                               groups = 5000,
                               min.size = 0)

obj <- AddMetaData(obj, metadata = trial.ssGSEA1)

jpeg("CITE_natcom_TRM_GSEA.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(obj, features = "gsea1", min.cutoff = 4000, reduction = "umap.rpca") + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# FigS1D
# data from PNAS 2023
sc_y <- readRDS('../yokoi_rap_scanorama_simul12sample.rds')

clean_gene_names <- function(gene_names) {
  cleaned <- sapply(gene_names, function(gene) {
    if (grepl("pAbO$", gene)) {
      return(gene) 
    } else {
      return(strsplit(gene, "\\.")[[1]][1]) 
    }
  })
  
  cleaned_unique <- make.unique(cleaned)
  
  return(cleaned_unique)
}

cleaned_genes <- clean_gene_names(rownames(sc_y@assays$RNA@counts))

head(cleaned_genes)

length(cleaned_genes) == nrow(sc_y@assays$RNA@counts)

rownames(sc_y@assays$RNA@counts) <- cleaned_genes
rownames(sc_y@assays$RNA@data) <- cleaned_genes
rownames(sc_y@assays$panorama@counts) <- cleaned_genes
rownames(sc_y@assays$panorama@data) <- cleaned_genes

reference <- readRDS("../CITEwithTCR.rds")
DefaultAssay(reference) <- "RNA"

sc_y[["RNA"]]@meta.features <- data.frame(row.names = rownames(sc_y[["RNA"]]))
DefaultAssay(sc_y) <- "RNA"

sc_y <- NormalizeData(sc_y, assay = "RNA", verbose = TRUE)
sc_y <- FindVariableFeatures(sc_y, assay = "RNA", selection.method = "vst", nfeatures = 2000)
sc_y <- ScaleData(sc_y, assay = "RNA", features = VariableFeatures(sc_y), verbose = TRUE)
sc_y <- RunPCA(
  object = sc_y,
  assay = "RNA",
  features = VariableFeatures(sc_y),
  verbose = TRUE
)

sc_y$panorama_RefMap <- as.character(sc_y$panorama_snn_res.0.8)
sc_y$panorama_RefMap[
  sc_y$panorama_snn_res.0.8 == 12 & sc_y$surgroups == "CD"
] <- "12-CD"


sc_y$panorama_RefMap[
  sc_y$panorama_snn_res.0.8 == 12 & sc_y$surgroups == "UC"
] <- "12-UC"

sc_y$panorama_RefMap[
  sc_y$panorama_snn_res.0.8 == 12 & sc_y$surgroups == "CK"
] <- "12-CK"

# extract CD4 T cells
Idents(sc_y) <- "panorama_snn_res.0.8"
sc_y_subset <- subset(sc_y, subset = panorama_snn_res.0.8 %in% c("4", "5", "2", "8", "0", "6", "1", "21", "12"))

LP5.anchors <- FindTransferAnchors(
  reference = sc_y_subset,
  query = reference,
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize", 
  reference.reduction = "pca",
  dims = 1:30
)

sc_y_subset <- RunUMAP(sc_y_subset, reduction = "pca", dims = 1:30, reduction.name = "umap", return.model = T)
obj <- MapQuery(
  anchorset = LP5.anchors,
  query = obj,
  reference = sc_y_subset,
  refdata = list(
    celltype.annotation = "panorama_RefMap"
  ),
  reference.reduction = "pca",
  reduction.model = "umap" 
)

jpeg("CITE_refmap_CD.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(
  obj, 
  features = c("12-CD"), 
  reduction = "umap.rpca"
)
dev.off()

jpeg("CITE_refmap_CK.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(
  obj, 
  features = c("12-CK"), 
  reduction = "umap.rpca"
)
dev.off()

jpeg("CITE_refmap_UC.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(
  obj, 
  features = c("12-UC"), 
  reduction = "umap.rpca"
)
dev.off()





