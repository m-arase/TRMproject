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

# CD31 : CD1
# Hy350 : CD2
# Hy392 : CD3
# CK265 : Ctrl1
# CK278 : Ctrl2
# CK282 : Ctrl3

counts <- Read10X_h5("/Volumes/BUFFALOHDD/MultiomeAnalysis/20230703_Multiome6sample/CD3CK3/outs/filtered_feature_bc_matrix.h5")
fragpath <- "/Volumes/BUFFALOHDD/MultiomeAnalysis/20230703_Multiome6sample/CD3CK3/outs/atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
LP <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
head(LP)
# head(str_sub(Cells(LP), start = -1, end = -1))
LP$disease <- "CK"
LP$disease[str_sub(Cells(LP), start = -1, end = -1) == 1] <- "CD"
LP$disease[str_sub(Cells(LP), start = -1, end = -1) == 5] <- "CD"
LP$disease[str_sub(Cells(LP), start = -1, end = -1) == 6] <- "CD"
LP$sample <- "CD31"
LP$sample[str_sub(Cells(LP), start = -1, end = -1) == 2] <- "CK265"
LP$sample[str_sub(Cells(LP), start = -1, end = -1) == 3] <- "CK278"
LP$sample[str_sub(Cells(LP), start = -1, end = -1) == 4] <- "CK282"
LP$sample[str_sub(Cells(LP), start = -1, end = -1) == 5] <- "Hy350"
LP$sample[str_sub(Cells(LP), start = -1, end = -1) == 6] <- "Hy392"

# create ATAC assay and add it to the object
LP[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
LP

DefaultAssay(LP) <- "ATAC"

LP <- NucleosomeSignal(LP)
LP <- TSSEnrichment(LP)

saveRDS(LP, file = "LPbeforeQC.rds")
LP <- readRDS(file = "LPbeforeQC.rds")

VlnPlot(
  object = LP,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells
LP <- subset(
  x = LP,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 30000 &
    nCount_ATAC > 500 &
    nCount_RNA > 500 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

# call peaks using MACS2
peaks <- CallPeaks(LP, macs2.path = "/opt/anaconda3/envs/r41python27/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(LP),
  features = peaks,
  cells = colnames(LP)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
LP[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(LP) <- "RNA"
LP[["percent.mt"]] <- PercentageFeatureSet(LP, pattern = "^MT-")
VlnPlot(LP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
LP <- subset(LP, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
LP.list <- SplitObject(LP, split.by = "sample")

# normalize and identify variable features for each dataset independently
LP.list <- lapply(X = LP.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = LP.list)

immune.anchors <- FindIntegrationAnchors(object.list = LP.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

DefaultAssay(immune.combined) <- "ATAC"

# process the combined dataset
immune.combined <- FindTopFeatures(immune.combined, min.cutoff = 10)
immune.combined <- RunTFIDF(immune.combined)
immune.combined <- RunSVD(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "lsi", reduction.name = "umap.atac.merge", dims = 2:30)

LP.list <- SplitObject(immune.combined, split.by = "sample")
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = LP.list,
  anchor.features = rownames(LP.list[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = immune.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", reduction.name = "umap.atac", dims = 2:30)
DimPlot(integrated, group.by = 'sample', reduction = "umap.atac", pt.size = 0.1)

integrated[["integrated.rna"]] <- immune.combined[["integrated"]]
integrated[["pca.integrated.rna"]] <- immune.combined[["pca"]]
integrated[["umap.integrated.rna"]] <- immune.combined[["umap"]]

integrated <- FindMultiModalNeighbors(integrated, reduction.list = list("pca.integrated.rna", "integrated_lsi"), dims.list = list(1:30, 2:30))
integrated <- RunUMAP(integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
integrated <- FindClusters(integrated, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.8)

saveRDS(integrated, file = "RNAandATACintegrate.rds")
Multiome_data <- readRDS(file = "/Volumes/BUFFALOHDD/MultiomeAnalysis/20230721_Multiome6sample/RNAandATACintegrate.rds")

# obj : CITE-seq data from Figure1
LP.anchors <- FindTransferAnchors(reference = obj, query = Multiome_data, dims = 1:50,
                                  reference.reduction = "integrated.rpca")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", return.model = T)
Multiome_data <- MapQuery(
  anchorset = LP.anchors,
  query = Multiome_data,
  reference = obj,
  refdata = list(
    celltype.annotation = "annotationnew",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "integrated.rpca", 
  reduction.model = "umap.rpca"
)

Idents(Multiome_data) <- "seurat_clusters"
newnames <- c("TEM_mixed", "TRM_1", "eTreg_2_1", "TRM_2", "Tfh_1_1", "Naive T_1", "eTreg_2_2", "Tfh_1_2", "Naive T_2",
              "TEM_2", "eTreg_2_3", "Mixed_1", "TEMRA_1", "eTreg_2_4", "eTreg_2_5", "Mixed_2", "Tfh_1_3", "Tfh_2", "Unannotatad_2_1",
              "Unannotatad_1_1", "Unannotatad_2_2", "eTreg_2_6", "Unannotatad_1_2")
names(newnames) <- levels(Multiome_data)
Multiome_data <- RenameIdents(Multiome_data, newnames)

Multiome_data$annotation <- Idents(Multiome_data)
head(Multiome_data)
DimPlot(Multiome_data, group.by = "annotation", reduction = "wnn.umap")

new_levels <- c("Naive T_1", "Naive T_2", "TEM_2", "TEM_mixed", "TRM_1", "TRM_2", "Tfh_1_1", "Tfh_1_2", "Tfh_1_3","Tfh_2","eTreg_2_1",
                "eTreg_2_2", "eTreg_2_3", "eTreg_2_4", "eTreg_2_5", "eTreg_2_6", "TEMRA_1", "Mixed_1", "Mixed_2", "Unannotatad_1_1",
                "Unannotatad_1_2", "Unannotatad_2_1", "Unannotatad_2_2")
Multiome_data$annotation <- factor(x = Multiome_data$annotation, levels = new_levels)

Multiome_data$celltype <- "Unannotated"
Multiome_data$celltype[Multiome_data$annotation == "Naive T_1"] <- "Naive T"
Multiome_data$celltype[Multiome_data$annotation == "Naive T_2"] <- "Naive T"

Multiome_data$celltype[Multiome_data$annotation == "TEM_2"] <- "TEM"
Multiome_data$celltype[Multiome_data$annotation == "TEM_mixed"] <- "TEM"

Multiome_data$celltype[Multiome_data$annotation == "TRM_1"] <- "TRM_1"
Multiome_data$celltype[Multiome_data$annotation == "TRM_2"] <- "TRM_2"

Multiome_data$celltype[Multiome_data$annotation == "Tfh_1_1"] <- "Tfh"
Multiome_data$celltype[Multiome_data$annotation == "Tfh_1_2"] <- "Tfh"
Multiome_data$celltype[Multiome_data$annotation == "Tfh_1_3"] <- "Tfh"
Multiome_data$celltype[Multiome_data$annotation == "Tfh_2"] <- "Tfh"

Multiome_data$celltype[Multiome_data$annotation == "eTreg_2_1"] <- "Treg"
Multiome_data$celltype[Multiome_data$annotation == "eTreg_2_2"] <- "Treg"
Multiome_data$celltype[Multiome_data$annotation == "eTreg_2_3"] <- "Treg"
Multiome_data$celltype[Multiome_data$annotation == "eTreg_2_4"] <- "Treg"
Multiome_data$celltype[Multiome_data$annotation == "eTreg_2_5"] <- "Treg"
Multiome_data$celltype[Multiome_data$annotation == "eTreg_2_6"] <- "Treg"

Multiome_data$celltype[Multiome_data$annotation == "TEMRA_1"] <- "TEMRA"

Multiome_data$celltype[Multiome_data$annotation == "Mixed_1"] <- "Mixed_1"
Multiome_data$celltype[Multiome_data$annotation == "Mixed_2"] <- "Mixed_2"

new_levels <- c("Naive T","TEM","TRM_1","TRM_2","Tfh","Treg","TEMRA","Mixed_1","Mixed_2","Unannotated")
Multiome_data$celltype <- factor(x = Multiome_data$celltype, levels = new_levels)

# Fig5A
cluster_colors4 <- c(
  "Naive T" = "#F8766D",  
  "TEM_2" = "#B79F00",  
  "TEM_mixed" = "#E58539", 
  "TRM_1" = "#0000FF",  
  "TRM_2" = "#FF0000",  
  "Tfh_1" = "#00BF7D", 
  "Tfh_2" = "#6BD76B", 
  "eTreg_2" = "#619CFF",  
  "TEMRA_1" = "#F564E3",  
  "Mixed_1" = "#A3A500",
  "Mixed_2" = "#00B0F6",
  "Unannotated" = "#A9A9A9"  
)

jpeg("Multiome_UMAP6.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Multiome_data, reduction = "wnn.umap", group.by = "celltype", cols = cluster_colors4, label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig5B
pdf("Multiome_UMAP_disease2.pdf", width = 6, height = 6)
DimPlot(Multiome_data, reduction = "wnn.umap", group.by = "disease", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig5C
pdf("Multiome_UMAP_CD103_ADT.pdf", width = 5, height = 5)
DefaultAssay(Multiome_data) <- "predicted_ADT"
FeaturePlot(Multiome_data, 
            features = c("CD103"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_HLADR_ADT.pdf", width = 5, height = 5)
DefaultAssay(Multiome_data) <- "predicted_ADT"
FeaturePlot(Multiome_data, 
            features = c("HLA-DR"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_ADT_Vln_CD103.pdf", width = 4, height = 2.5)
DefaultAssay(Multiome_data) <- "predicted_ADT"
VlnPlot(Multiome_data, features = "CD103", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "eTreg", "Tfh", "Naive T"),
        pt.size = 0)+ NoLegend()
dev.off()

pdf("Multiome_ADT_Vln_HLADR.pdf", width = 4, height = 2.5)
DefaultAssay(Multiome_data) <- "predicted_ADT"
VlnPlot(Multiome_data, features = "HLA-DR", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "eTreg", "Tfh", "Naive T"),
        pt.size = 0)+ NoLegend()
dev.off()

# Fig5D
pdf("Multiome_UMAP_IFNG.pdf", width = 5, height = 5)
FeaturePlot(Multiome_data, 
            features = c("IFNG"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_GZMB.pdf", width = 5, height = 5)
FeaturePlot(Multiome_data, 
            features = c("GZMB"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_RNA_Vln_IFNG.pdf", width = 4, height = 2.5)
DefaultAssay(Multiome_data) <- "RNA"
VlnPlot(Multiome_data, features = "IFNG", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "eTreg", "Tfh", "Naive T"),
        pt.size = 0)+ NoLegend()
dev.off()

pdf("Multiome_RNA_Vln_GZMB.pdf", width = 4, height = 2.5)
DefaultAssay(Multiome_data) <- "RNA"
VlnPlot(Multiome_data, features = "GZMB", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "eTreg", "Tfh", "Naive T"),
        pt.size = 0)+ NoLegend()
dev.off()

# Fig5E
idents.plot <- c("TRM_2", "TRM_1", "naive T_2", "TEM_mixed")
DefaultAssay(Multiome_data) <- "peaks"

pdf("Coverageplot_IFNG.pdf", width = 6, height = 3)
CoveragePlot(
  object = Multiome_data,
  region = "IFNG",
  idents = idents.plot,
  extend.upstream = 5000,
  extend.downstream = 5000
)
dev.off()

pdf("Coverageplot_GZMB.pdf", width = 6, height = 3)
CoveragePlot(
  object = Multiome_data,
  region = "GZMB",
  idents = idents.plot,
  extend.upstream = 5000,
  extend.downstream = 5000
)
dev.off()

# Fig5F
DefaultAssay(Multiome_data) <- "peaks"
da_peaks <- FindMarkers(
  object = Multiome_data,
  ident.1 = 'TRM_2',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

write.csv(da_peaks, file = "TRM_2_peaks.csv")

top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])
enriched.motifs <- FindMotifs(
  object = Multiome_data,
  features = top.da.peak
)

enriched.motifs <- read.csv(file = "TRM_2_motifs.csv", row.names = 1)
result <- enriched.motifs

result$plot <- "NO"
result$plot[result$motif.name == "RUNX2"] <- "YES"
result$plot[result$motif.name == "BHLHE40"] <- "YES"

result$delabel <- NA
result$delabel[result$plot != "NO"] <- result$motif.name[result$plot != "NO"]
result$volcano <- -log10(result$pvalue)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)

pdf("TRM_2_motif_volcano_RUNX2BHLHE40.pdf", width = 6, height = 6)
ggplot(data=result, aes(x=fold.enrichment, y=volcano, col=plot, label=delabel)) +
  geom_point(show.legend = FALSE) + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values=c("gray", "red")) +
  scale_y_continuous(limits=c(0,60),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

# ChIP-Atlas
chipatlas <- read.table(file = "TRM_2_top500_CHIPAtlas_delete_common.txt",  sep="\t")
row.names(chipatlas) <- paste0(chipatlas$V3, "_",chipatlas$V1)

chipatlas$plot <- "NO"
chipatlas$plot[chipatlas$V3 == "RUNX2"] <- "YES"
chipatlas$plot[chipatlas$V3 == "BHLHE40"] <- "YES"

chipatlas$delabel[chipatlas$plot != "NO"] <- chipatlas$V3[chipatlas$plot != "NO"]
chipatlas$volcano <- -(chipatlas$V9)
chipatlas["volcano"] <- lapply(chipatlas["volcano"], gsub, pattern="Inf", replacement = 300)
chipatlas["V11"] <- lapply(chipatlas["V11"], gsub, pattern="Inf", replacement = 363)
chipatlas["volcano"] <- lapply(chipatlas["volcano"], as.double)
chipatlas["V11"] <- lapply(chipatlas["V11"], as.double)


pdf("TRM_2_chipaltas_volcano_RUNX2BHLHE40.pdf", width = 6, height = 6)
ggplot(data=chipatlas, aes(x=V11, y=volcano, col=plot, label=delabel)) +
  geom_point(show.legend = FALSE) + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values=c("gray", "red")) +
  scale_y_continuous(limits=c(0,170),expand=c(0,0)) +
  xlab("log2_fold_change") +
  ylab("-log10_p_val")
dev.off()

# TFonly_DEG
Multiome_data <- ScaleData(Multiome_data)
result <- FindMarkers(Multiome_data, ident.1 = "TRM_2", group.by = "annotation", 
                      min.pct = 0, logfc.threshold = 0)
write.csv(result, file = "TRM_2_markers.csv")

# extract TF
result <- read.csv(file = "TRM_2_markers_TFonly.csv")
result$plot <- "NO"
result$plot[result$X == "RUNX2"] <- "YES"
result$plot[result$X == "BHLHE40"] <- "YES"
result$delabel <- NA
result$delabel[result$plot != "NO"] <- result$X[result$plot != "NO"]
result$volcano <- -log10(result$p_val_adj)
result["volcano"] <- lapply(result["volcano"], gsub, pattern="Inf", replacement = 300)
result["volcano"] <- lapply(result["volcano"], as.double)
pdf("TRM_2_TFonly_DEG_volcano_RUNX2BHLHE40.pdf", width = 3, height = 3)
ggplot(data=result, aes(x=avg_log2FC, y=volcano, col=plot, label=delabel)) +
  geom_point(show.legend = FALSE) + 
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  theme(axis.line = element_line(colour="black")) +
  geom_text_repel(show.legend = FALSE) +
  scale_color_manual(values=c("gray", "red")) +
  scale_y_continuous(limits=c(0,310),expand=c(0,0)) +
  ylab("-log10_p_val")
dev.off()

# Fig5J
p1 <- FeaturePlot(Multiome_data, 
                  features = c("RBPJ"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p2 <- FeaturePlot(Multiome_data, 
                  features = c("PPARG"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p3 <- FeaturePlot(Multiome_data, 
                  features = c("BHLHE40"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p4 <- FeaturePlot(Multiome_data, 
                  features = c("RUNX2"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p5 <- FeaturePlot(Multiome_data, 
                  features = c("CEBPD"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p6 <- FeaturePlot(Multiome_data, 
                  features = c("MYB"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p7 <- FeaturePlot(Multiome_data, 
                  features = c("STAT1"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p8 <- FeaturePlot(Multiome_data, 
                  features = c("IRF4"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p9 <- FeaturePlot(Multiome_data, 
                  features = c("VDR"), 
                  reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p10 <- FeaturePlot(Multiome_data, 
                   features = c("TFEB"), 
                   reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p11 <- FeaturePlot(Multiome_data, 
                   features = c("TBX21"), 
                   reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p12 <- FeaturePlot(Multiome_data, 
                   features = c("BATF3"), 
                   reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p13 <- FeaturePlot(Multiome_data, 
                   features = c("GFI1"), 
                   reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p14 <- FeaturePlot(Multiome_data, 
                   features = c("TFDP1"), 
                   reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
p15 <- FeaturePlot(Multiome_data, 
                   features = c("TWIST1"), 
                   reduction = "wnn.umap") + NoLegend() + NoAxes() + ggtitle(NULL)
jpeg("Multiome_features_top.jpg", width = 4096, height = 512, res = 400)
p1 | p2 | p3 | p4 | p5 | p6 | p7 | p8
dev.off()
jpeg("Multiome_features_bot.jpg", width = 3584, height = 512, res = 400)
p9 | p10 | p11 | p12 | p13 | p14 | p15
dev.off()

# Fig5K
RNAlist <- c("RUNX2", "BHLHE40", "RBPJ", "PPARG",  "CEBPD", "MYB", "STAT1",
             "IRF4", "VDR", "TFEB", "TBX21", "BATF3", "GFI1", "TFDP1", "TWIST1")
DefaultAssay(Multiome_data) <- "RNA"
new_levels <- c("Naive T","TEM","TRM_1","TRM_2","Tfh","Treg","TEMRA","Mixed_1","Mixed_2","Unannotated")
new_levels <- rev(new_levels)
Multiome_data$celltype <- factor(x = Multiome_data$celltype, levels = new_levels)
pdf("Multiome_RNA_TF_Dot.pdf", width = 12, height = 6)
DotPlot(subset(Multiome_data, idents = c("Naive T","TEM","TRM_1","TRM_2","Tfh","Treg","TEMRA","Mixed_1","Mixed_2")), 
        features = RNAlist, dot.scale=10, cols="RdBu",dot.min = 0,group.by = "celltype") + 
  RotatedAxis()
dev.off()

