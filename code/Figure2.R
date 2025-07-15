library(Seurat)
library(hdf5r)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)
library(harmony)
library(SingleCellExperiment)
set.seed(1234)

LPhash1 <- readRDS(file = '/Volumes/BUFFALOHDD/scRNAseq/20220215_CD103_hashtag/20220217_Analysis_Demultiplexing/LPhashbeforeQC.rds')
LPhash2 <- readRDS(file = '/Volumes/BUFFALOHDD/scRNAseq/20220324_CD103_hashtag/20220324_Analysis/LPhashbeforeQC.rds')

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
LPhash1 <- NormalizeData(LPhash1, assay = "HTO", normalization.method = "CLR")
LPhash2 <- NormalizeData(LPhash2, assay = "HTO", normalization.method = "CLR")
# Demultiplex cells based on HTO enrichment
LPhash1 <- HTODemux(LPhash1, assay = "HTO", positive.quantile = 0.99)
LPhash2 <- HTODemux(LPhash2, assay = "HTO", positive.quantile = 0.99)

# Extract the singlets
Idents(LPhash1) <- "HTO_classification.global"
LPhash1 <- subset(LPhash1, idents = "Singlet")
Idents(LPhash2) <- "HTO_classification.global"
LPhash2 <- subset(LPhash2, idents = "Singlet")

head(LPhash1)
# rename samples
LPhash1$sample <- 'CK206'
LPhash1$sample[LPhash1$HTO_classification == 'Hashtag2'] <- 'CK231'
LPhash1$sample[LPhash1$HTO_classification == 'Hashtag3'] <- 'CK238'
LPhash1$sample[LPhash1$HTO_classification == 'Hashtag4'] <- 'CD253'
LPhash1$sample[LPhash1$HTO_classification == 'Hashtag5'] <- 'CD259'
LPhash1$sample[LPhash1$HTO_classification == 'Hashtag6'] <- 'CD268'

LPhash2$sample <- 'CK221'
LPhash2$sample[LPhash2$HTO_classification == 'Hashtag2'] <- 'CK228'
LPhash2$sample[LPhash2$HTO_classification == 'Hashtag3'] <- 'CK242'
LPhash2$sample[LPhash2$HTO_classification == 'Hashtag4'] <- 'CD24'
LPhash2$sample[LPhash2$HTO_classification == 'Hashtag5'] <- 'CD27'
LPhash2$sample[LPhash2$HTO_classification == 'Hashtag6'] <- 'CD189'

LPhash$batch <- 'batch1'
LPhash$batch[LPhash$sample == 'CK221'] <- 'batch2'
LPhash$batch[LPhash$sample == 'CK228'] <- 'batch2'
LPhash$batch[LPhash$sample == 'CK242'] <- 'batch2'
LPhash$batch[LPhash$sample == 'CD24'] <- 'batch2'
LPhash$batch[LPhash$sample == 'CD27'] <- 'batch2'
LPhash$batch[LPhash$sample == 'CD189'] <- 'batch2'

LPhash$disease <- 'CK'
LPhash$disease[LPhash$sample == 'CD189'] <- 'CD'
LPhash$disease[LPhash$sample == 'CD24'] <- 'CD'
LPhash$disease[LPhash$sample == 'CD253'] <- 'CD'
LPhash$disease[LPhash$sample == 'CD259'] <- 'CD'
LPhash$disease[LPhash$sample == 'CD27'] <- 'CD'
LPhash$disease[LPhash$sample == 'CD268'] <- 'CD'

# -----------------------just merge
LPhash <- merge(LPhash1, y = LPhash2, add.cell.ids = c("hash1", "hash2"), project = "LPhash")
head(LPhash)
LPhash <- subset(LPhash, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
table(LPhash$sample)

LPhash <- NormalizeData(LPhash, normalization.method = "LogNormalize", scale.factor = 10000)
LPhash <- FindVariableFeatures(LPhash, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(LPhash)
LPhash <- ScaleData(LPhash, features = all.genes)

LPhash <- RunPCA(LPhash, features = VariableFeatures(object = LPhash))
ElbowPlot(LPhash)

LPhash <- RunHarmony(LPhash, "batch")
LPhash <- RunUMAP(LPhash, reduction = "harmony", reduction.name = "harmony.umap", dims = 1:20)
LPhash <- FindNeighbors(LPhash, reduction = "harmony", dims = 1:20)
LPhash <- FindClusters(LPhash, resolution = 0.5)

# LP.markers <- FindAllMarkers(LPhash, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(LP.markers, "clustermarkers.csv")

# Add TCR for FigS2 
tcr1 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220215_CD103_hashtag/F3804_220215_152314_cellranger/LPMC_TCR/outs/filtered_contig_annotations.csv")
tcr2 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220324_CD103_hashtag/F3970_220323_162237_cellranger/CD103_TCR/outs/filtered_contig_annotations.csv")
head(tcr1)
# Clonotype-centric info.
clono1 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220215_CD103_hashtag/F3804_220215_152314_cellranger/LPMC_TCR/outs/clonotypes.csv")
clono2 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220324_CD103_hashtag/F3970_220323_162237_cellranger/CD103_TCR/outs/clonotypes.csv")
head(clono1)

tcr1 <- tcr1[!duplicated(tcr1$barcode), ]
tcr2 <- tcr2[!duplicated(tcr2$barcode), ]

tcr1 <- tcr1[,c("barcode", "raw_clonotype_id")]
head(tcr1)
names(tcr1)[names(tcr1) == "raw_clonotype_id"] <- "clonotype_id"
tcr1[["barcode"]] <- paste0("hash1_", tcr1[["barcode"]])
tcr2 <- tcr2[,c("barcode", "raw_clonotype_id")]
head(tcr2)
names(tcr2)[names(tcr2) == "raw_clonotype_id"] <- "clonotype_id"
tcr2[["barcode"]] <- paste0("hash2_", tcr2[["barcode"]])

# Slap the AA sequences onto our original table by clonotype_id.
tcr1 <- merge(tcr1, clono1[, c("clonotype_id", "cdr3s_aa")])
tcr2 <- merge(tcr2, clono2[, c("clonotype_id", "cdr3s_aa")])
# Reorder so barcodes are first column and set them as rownames.
tcr1 <- tcr1[, c(2,1,3)]
rownames(tcr1) <- tcr1[,1]
tcr1[,1] <- NULL
tcr2 <- tcr2[, c(2,1,3)]
rownames(tcr2) <- tcr2[,1]
tcr2[,1] <- NULL
head(tcr2)
tcr1[["clonotype_id"]] <- paste0(tcr1[['clonotype_id']], '_1')
tcr2[["clonotype_id"]] <- paste0(tcr2[['clonotype_id']], '_2')
head(tcr1)

head(LPhash)
tcr <- rbind(tcr1, tcr2)
tail(tcr)

LPhash2 <- AddMetaData(object=LPhash, metadata=tcr)

CD103LP2 <- subset(x = LPhash, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))

levels(Idents(CD103LP2))
Idents(CD103LP2) <- factor(Idents(CD103LP2), levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
DimPlot(CD103LP2, reduction = "harmony.umap", label = T, group.by = "RNA_snn_res.0.5")
new_names <- c("TXNIP high", "HSP high_1", "FOSB high", "CCL20 high",
               "Cytotoxic", "MT high", "HSP high_2", "Treg", "NK_like", "naive")
names(new_names) <- levels(CD103LP2)
CD103LP2 <- RenameIdents(CD103LP2, new_names)
CD103LP2$annotationnew <- Idents(CD103LP2)
new_levels <- c("TXNIP high","CCL20 high","Cytotoxic", "FOSB high", "HSP high_1","HSP high_2","Treg", "NK_like", "naive","MT high")
CD103LP2$annotationnew <- factor(x = CD103LP2$annotationnew, levels = new_levels)

# Fig2B
pdf("CD103_UMAP.pdf", width = 6, height = 6)
DimPlot(CD103LP2, reduction = "harmony.umap", label = F, group.by = "annotationnew") + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig2C
dotfeaturesRNA <- c("TXNIP", "CXCR4", "HSPA1A", "HSPD1", "HSPA6", "FOSB", "JUN", "CCL20", "TRAF1",
                    "HLA-DRB5", "IFNG", "GZMB", "GZMK", "NKG7", "TIGIT", "FOXP3", "CD27", "IKZF2", "LEF1", "SELL", "CCR7")

CD103LP2$annotation <- factor(x = CD103LP2$annotation, levels = rev(new_names))

pdf("CD103_RNA_Dot.pdf", width = 11, height = 5)
DotPlot(CD103LP2, features = dotfeaturesRNA, dot.scale=10, group.by = "annotationnew", cols="RdBu") + 
  RotatedAxis()
dev.off()

# Fig2D
pdf("CD103_UMAP_disease.pdf", width = 6, height = 6)
DimPlot(CD103LP2, reduction = "harmony.umap", label = F, group.by = "disease") + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig2E
write.csv(prop.table(table(CD103LP2$disease, CD103LP2$annotation), margin = 1), file = "CD103_celltype_per_disease_prop.csv")

# Fig2F
DefaultAssay(CD103LP2) <- "RNA"
pdf("CD103_IFNG.pdf", width = 3, height = 3)
FeaturePlot(CD103LP2, 
            features = c("IFNG"), 
            reduction = "harmony.umap") + NoAxes() 
dev.off()

pdf("CD103_GZMB.pdf", width = 3, height = 3)
FeaturePlot(CD103LP2, 
            features = c("GZMB"), 
            reduction = "harmony.umap") + NoAxes() 
dev.off()

pdf("CD103_GZMH.pdf", width = 3, height = 3)
FeaturePlot(CD103LP2, 
            features = c("GZMH"), 
            reduction = "harmony.umap") + NoAxes() 
dev.off()

pdf("CD103_PRF1.pdf", width = 3, height = 3)
FeaturePlot(CD103LP2, 
            features = c("PRF1"), 
            reduction = "harmony.umap") + NoAxes() 
dev.off()

pdf("CD103_HLA-DRB1.pdf", width = 3, height = 3)
FeaturePlot(CD103LP2, 
            features = c("HLA-DRB1"), 
            reduction = "harmony.umap") + NoAxes() 
dev.off()

pdf("CD103_ENTPD1.pdf", width = 3, height = 3)
FeaturePlot(CD103LP2, 
            features = c("ENTPD1"), 
            reduction = "harmony.umap") + NoAxes() 
dev.off()





