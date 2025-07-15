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
library(escape)

set.seed(1234)
# CITE-seq data
obj <- readRDS(file = "CITE_rpca_annotation.rds")
# Multiome data
Multiome_data <- readRDS(file = "/Volumes/BUFFALOHDD/CITEseqAnalysis/20240611_referencemapping/Multiome20240701_annotation.rds")

# Fig7G
# Both OE genes
gene.sets1 <- read.csv("/Volumes/BUFFALOHDD/lab_results/bulkRNAseq/Lenti_Runx2_Bhlhe40_overexpression/RUBH_fc2_pval0.05.csv")
gene.sets1 <- list(gene.sets1[1:100,1])
names(gene.sets1) <- "gseaRUBH"

ssGSEA1 <- escape.matrix(obj@assays$RNA$counts, 
                              method = "ssGSEA", 
                              gene.sets = gene.sets1, 
                              groups = 5000,
                              min.size = 0)
obj <- AddMetaData(obj, metadata = ssGSEA1)

# RUNX2 OE genes
gene.sets3 <- read.csv("/Volumes/BUFFALOHDD/lab_results/bulkRNAseq/Lenti_Runx2_Bhlhe40_overexpression/RU_fc2_pval0.05.csv")
gene.sets3 <- list(gene.sets3[1:100,1])
names(gene.sets3) <- "gseaRU"
ssGSEA2 <- escape.matrix(obj@assays$RNA$counts, 
                               method = "ssGSEA", 
                               gene.sets = gene.sets3, 
                               groups = 5000,
                               min.size = 0)

obj <- AddMetaData(obj, metadata = ssGSEA2)

# BHLHE40 OE genes
gene.sets5 <- read.csv("/Volumes/BUFFALOHDD/lab_results/bulkRNAseq/Lenti_Runx2_Bhlhe40_overexpression/BH_fc2_pval0.05.csv")

gene.sets5 <- list(gene.sets5[1:100,1])
names(gene.sets5) <- "gseaBH"
ssGSEA3 <- escape.matrix(obj@assays$RNA$counts, 
                               method = "ssGSEA", 
                               gene.sets = gene.sets5, 
                               groups = 5000,
                               min.size = 0)
obj <- AddMetaData(obj, metadata = ssGSEA3)

jpeg("CITE_RUNX2BHLHE40_scGSEA.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(obj, features = "gseaRUBH", reduction = "umap.rpca", min.cutoff = 0, max.cutoff = 1000) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()
jpeg("CITE_RUNX2_scGSEA.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(obj, features = "gseaRU", reduction = "umap.rpca", min.cutoff = 0, max.cutoff = 1000) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()
jpeg("CITE_BHLHE40_scGSEA.jpg", width = 2048, height = 2048, res = 400)
FeaturePlot(obj, features = "gseaBH", reduction = "umap.rpca", min.cutoff = 0, max.cutoff = 1000) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig7J
new_levels <- c("Naive T","TEM","TRM_1","TRM_2","Tfh","Treg","TEMRA","Mixed_1","Mixed_2","Unannotated")
new_levels <- rev(new_levels)
Multiome_data$celltype <- factor(x = Multiome_data$celltype, levels = new_levels)

DefaultAssay(Multiome_data) <- "ATAC"
Multiome_data <- RegionStats(Multiome_data, genome = BSgenome.Hsapiens.UCSC.hg38)
Multiome_data <- LinkPeaks(
  object = Multiome_data,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = c("BHLHE40")
)
CoveragePlot(
  object = Multiome_data,
  region = "BHLHE40",
  features = "BHLHE40",
  idents = idents.plot,
  extend.upstream = 20000,
  extend.downstream = 30000
)
