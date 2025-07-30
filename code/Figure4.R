library(Seurat)
library(hdf5r)
library(ggrepel)
library(patchwork)
library(cowplot)
library(tidyverse)
library(scRepertoire)
library(circlize)
library(scales)
library(ggraph)
library(UpSetR)
library(harmony)
set.seed(1234)

#------------------------sample1

LP.umis1 <- Read10X(data.dir = "/Users/mitsuru/Documents/lab/20250115_blood_lymph_colon_CITE/20241220/outs/per_sample_outs/20241220/count/sample_filtered_feature_bc_matrix")

rownames(LP.umis1[["Antibody Capture"]])[
  rownames(LP.umis1[["Antibody Capture"]]) == "Hashtag_1"
] <- "Hashtag1"
rownames(LP.umis1[["Antibody Capture"]])[
  rownames(LP.umis1[["Antibody Capture"]]) == "Hashtag_2"
] <- "Hashtag2"
rownames(LP.umis1[["Antibody Capture"]])[
  rownames(LP.umis1[["Antibody Capture"]]) == "Hashtag_3"
] <- "Hashtag3"
rownames(LP.umis1[["Antibody Capture"]])[
  rownames(LP.umis1[["Antibody Capture"]]) == "Hashtag_4"
] <- "Hashtag4"
rownames(LP.umis1[["Antibody Capture"]])

# Setup Seurat object
LPhash1 <- CreateSeuratObject(counts = LP.umis1$`Gene Expression`)

DefaultAssay(LPhash1) <- "RNA"
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
LPhash1[["percent.mt"]] <- PercentageFeatureSet(LPhash1, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(LPhash1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Add HTO data as a new assay independent from RNA
LPhash1[["HTO"]] <- CreateAssayObject(counts = LP.umis1$`Antibody Capture`[1:4,])
LPhash1[["ADT"]] <- CreateAssayObject(counts = LP.umis1$`Antibody Capture`[5:nrow(LP.umis1$`Antibody Capture`), ])

DefaultAssay(LPhash1) <- "HTO"
LPhash1 <- NormalizeData(LPhash1, assay = "HTO", normalization.method = "CLR")
# Demultiplex cells based on HTO enrichment
LPhash1 <- HTODemux(LPhash1, assay = "HTO", positive.quantile = 0.99)
#head(LPhash1)
table(LPhash1$HTO_classification)
table(LPhash1$HTO_classification.global)

# add TCR
TCR1 <- read.csv("/Users/mitsuru/Documents/lab/20250115_blood_lymph_colon_CITE/20241220/outs/per_sample_outs/20241220/vdj_t/filtered_contig_annotations.csv")
contig.list1 <- createHTOContigList(TCR1, 
                                    LPhash1, 
                                    group.by = "HTO_classification")

combined.TCR1 <- combineTCR(contig.list1,
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = FALSE)

scRep1 <- combineExpression(combined.TCR1, 
                            LPhash1, 
                            cloneCall="gene", 
                            proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

addTCR <- function(so, tcrpath, clonopath){
  TCR <- read.csv(tcrpath)
  clono <- read.csv(clonopath)
  TCR <- TCR[!duplicated(TCR$barcode), ]
  TCR <- TCR[,c("barcode", "raw_clonotype_id")]
  names(TCR)[names(TCR) == "raw_clonotype_id"] <- "clonotype_id"
  TCR <- merge(TCR, clono[, c("clonotype_id", "cdr3s_aa")])
  TCR <- TCR[, c(2,1,3)]
  rownames(TCR) <- TCR[,1]
  TCR[,1] <- NULL
  TCR[["clonotype_id"]] <- paste0(TCR[['clonotype_id']], '_', as.character(eval(substitute(alist(so)))))
  so <- AddMetaData(object=so, metadata=TCR)
  return(so)
}

tcr_1 <- "/Users/mitsuru/Documents/lab/20250115_blood_lymph_colon_CITE/20241220/outs/per_sample_outs/20241220/vdj_t/filtered_contig_annotations.csv"
clono1 <- "/Users/mitsuru/Documents/lab/20250115_blood_lymph_colon_CITE/20241220/outs/per_sample_outs/20241220/vdj_t/clonotypes.csv"

scRep1 <- addTCR(scRep1, tcr_1, clono1)

#------------------------sample2

LP.umis2 <- Read10X(data.dir = "/Users/mitsuru/Documents/lab/20250321_CITE/20250220/outs/per_sample_outs/20250220/count/sample_filtered_feature_bc_matrix")
rownames(LP.umis2[["Antibody Capture"]])
rownames(LP.umis2[["Antibody Capture"]])[
  rownames(LP.umis2[["Antibody Capture"]]) == "CD103"
] <- "anti_human_CD103"
# Setup Seurat object
LPhash2 <- CreateSeuratObject(counts = LP.umis2$`Gene Expression`)
DefaultAssay(LPhash2) <- "RNA"
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
LPhash2[["percent.mt"]] <- PercentageFeatureSet(LPhash2, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(LPhash2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Add HTO data as a new assay independent from RNA
LPhash2[["HTO"]] <- CreateAssayObject(counts = LP.umis2$`Antibody Capture`[1:4,])
LPhash2[["ADT"]] <- CreateAssayObject(counts = LP.umis2$`Antibody Capture`[1:5,])

DefaultAssay(LPhash2) <- "HTO"
LPhash2 <- NormalizeData(LPhash2, assay = "HTO", normalization.method = "CLR")
# Demultiplex cells based on HTO enrichment
LPhash2 <- HTODemux(LPhash2, assay = "HTO", positive.quantile = 0.99)
#head(LPhash2)
table(LPhash2$HTO_classification)
table(LPhash2$HTO_classification.global)

# add TCR
TCR2 <- read.csv("/Users/mitsuru/Documents/lab/20250321_CITE/20250220/outs/per_sample_outs/20250220/vdj_t/filtered_contig_annotations.csv")
contig.list2 <- createHTOContigList(TCR2, 
                                    LPhash2, 
                                    group.by = "HTO_classification")

combined.TCR2 <- combineTCR(contig.list2,
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = FALSE)

scRep2 <- combineExpression(combined.TCR2, 
                            LPhash2, 
                            cloneCall="gene", 
                            proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))


tcr_2 <- "/Users/mitsuru/Documents/lab/20250321_CITE/20250220/outs/per_sample_outs/20250220/vdj_t/filtered_contig_annotations.csv"
clono2 <- "/Users/mitsuru/Documents/lab/20250321_CITE/20250220/outs/per_sample_outs/20250220/vdj_t/clonotypes.csv"

scRep2 <- addTCR(scRep2, tcr_2, clono2)

#------------------------sample3

LP.umis3 <- Read10X(data.dir = "/Users/mitsuru/Documents/lab/20250422_CITE/20250328/outs/per_sample_outs/20250328/count/sample_filtered_feature_bc_matrix")
rownames(LP.umis3[["Antibody Capture"]])[
  rownames(LP.umis3[["Antibody Capture"]]) == "CD103"
] <- "anti_human_CD103"
# Setup Seurat object
LPhash3 <- CreateSeuratObject(counts = LP.umis3$`Gene Expression`)
DefaultAssay(LPhash3) <- "RNA"
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
LPhash3[["percent.mt"]] <- PercentageFeatureSet(LPhash3, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(LPhash3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# Add HTO data as a new assay independent from RNA
LPhash3[["HTO"]] <- CreateAssayObject(counts = LP.umis3$`Antibody Capture`[1:4,])
LPhash3[["ADT"]] <- CreateAssayObject(counts = LP.umis3$`Antibody Capture`[1:5,])

DefaultAssay(LPhash3) <- "HTO"
LPhash3 <- NormalizeData(LPhash3, assay = "HTO", normalization.method = "CLR")
# Demultiplex cells based on HTO enrichment
LPhash3 <- HTODemux(LPhash3, assay = "HTO", positive.quantile = 0.99)
#head(LPhash3)
table(LPhash3$HTO_classification)
table(LPhash3$HTO_classification.global)

# add TCR
TCR3 <- read.csv("/Users/mitsuru/Documents/lab/20250422_CITE/20250328/outs/per_sample_outs/20250328/vdj_t/filtered_contig_annotations.csv")
contig.list3 <- createHTOContigList(TCR3, 
                                    LPhash3, 
                                    group.by = "HTO_classification")

combined.TCR3 <- combineTCR(contig.list3,
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = FALSE)

scRep3 <- combineExpression(combined.TCR3, 
                            LPhash3, 
                            cloneCall="gene", 
                            proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))


tcr_3 <- "/Users/mitsuru/Documents/lab/20250422_CITE/20250328/outs/per_sample_outs/20250328/vdj_t/filtered_contig_annotations.csv"
clono3 <- "/Users/mitsuru/Documents/lab/20250422_CITE/20250328/outs/per_sample_outs/20250328/vdj_t/clonotypes.csv"

scRep3 <- addTCR(scRep3, tcr_3, clono3)

scRep1$batch <- "batch1"
scRep2$batch <- "batch2"
scRep3$batch <- "batch3"

scRep1$HTO_batch <- paste0(scRep1$HTO_classification, "_", scRep1$batch)
scRep2$HTO_batch <- paste0(scRep2$HTO_classification, "_", scRep2$batch)
scRep3$HTO_batch <- paste0(scRep3$HTO_classification, "_", scRep3$batch)

scRep1 <- subset(scRep1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
scRep2 <- subset(scRep2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)
scRep3 <- subset(scRep3, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15)

combined <- merge(scRep1, y = list(scRep2, scRep3), add.cell.ids = c("batch1", "batch2", "batch3"))

DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

combined <- subset(combined, subset = HTO_classification.global == "Singlet")
combined <- RunHarmony(combined, "HTO_batch")

combined <- combined %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.8)

combined <- RunUMAP(combined, reduction = "harmony", reduction.name = "harmony.umap", dims = 1:30)

DefaultAssay(combined) <- "ADT"
combined <- NormalizeData(combined, assay = "ADT", normalization.method = "CLR", margin = 2)

combined$tissue <- "PBMC"
combined$tissue[combined$HTO_classification == "Hashtag1"] <- "uninflamed"
combined$tissue[combined$HTO_classification == "Hashtag2"] <- "inflamed"
combined$tissue[combined$HTO_classification == "Hashtag3"] <- "Lymphnode"


reference <- readRDS(file = "/Volumes/BUFFALOHDD/CITEseqAnalysis/20240724_Figure/CITEwithTCR.rds")
DefaultAssay(combined) <- "RNA"
DefaultAssay(reference) <- "RNA"
LP.anchors <- FindTransferAnchors(reference = reference, query = combined, dims = 1:50, reference.reduction = "integrated.rpca")
head(reference)
reference <- RunUMAP(reference, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", return.model = T)

new.cluster.ids <- c("TEM_1", "eTreg_1", "Tfh_1", "TRM_1", "naïve T", "TRM_2", 
                     "TCM", "Tfh_2", "TEM_2", "naïve Treg", "eTreg_2",
                     "TEMRA_1", "TEMRA_2", "unannotated_1", "unannotated_2")

Idents(reference)
names(new.cluster.ids) <- levels(reference)
reference <- RenameIdents(reference, new.cluster.ids)
reference$annotationnew <- Idents(reference)

new_levels <- c("naïve T","TEM_1","TEM_2", "TCM", "TRM_1","TRM_2","Tfh_1","Tfh_2","naïve Treg","eTreg_1","eTreg_2",
                "TEMRA_1","TEMRA_2","unannotated_1", "unannotated_2")

reference$annotationnew <- factor(x = reference$annotationnew, levels = new_levels)

Mapped_data <- MapQuery(
  anchorset = LP.anchors,
  query = combined,
  reference = reference,
  refdata = list(
    celltype.annotation = "annotationnew",
    arase_ADT = "ADT"
  ),
  reference.reduction = "integrated.rpca", 
  reduction.model = "umap.rpca"
)

saveRDS(Mapped_data, file = "CITE_TCR_harmony_mapped.rds")

Mapped_data <- readRDS(file = "CITE_TCR_harmony_mapped.rds")

Mapped_data <- subset(Mapped_data, subset = RNA_snn_res.0.8 != "13")
plot <- DimPlot(Mapped_data, reduction = "harmony.umap") + NoLegend()
cells_to_remove <- CellSelector(plot)
Mapped_data_filtered <- subset(Mapped_data, cells = setdiff(Cells(Mapped_data), cells_to_remove))
DimPlot(Mapped_data_filtered, reduction = "harmony.umap")

uninflamed <- subset(Mapped_data_filtered, subset = tissue == "uninflamed")
inflamed <- subset(Mapped_data_filtered, subset = tissue == "inflamed")
Lymphnode <- subset(Mapped_data_filtered, subset = tissue == "Lymphnode")
PBMC <- subset(Mapped_data_filtered, subset = tissue == "PBMC")

# uninflamed
DefaultAssay(uninflamed) <- "RNA"
uninflamed <- NormalizeData(uninflamed)
uninflamed <- FindVariableFeatures(uninflamed)
uninflamed <- ScaleData(uninflamed)
uninflamed <- RunPCA(uninflamed, reduction.name = "pca_uninf")
#head(uninflamed)
uninflamed <- RunHarmony(uninflamed, "batch", reduction = "pca_uninf", reduction.save = "uninf_harmony")
uninflamed <- FindNeighbors(uninflamed, reduction = "uninf_harmony")
uninflamed <- FindClusters(uninflamed, resolution = 1.0)

uninflamed <- RunUMAP(uninflamed, reduction = "uninf_harmony", reduction.name = "uninf_harmony.umap", dims = 1:30)
DimPlot(uninflamed, reduction = "uninf_harmony.umap", label = TRUE)

#inflamed
DefaultAssay(inflamed) <- "RNA"
inflamed <- NormalizeData(inflamed)
inflamed <- FindVariableFeatures(inflamed)
inflamed <- ScaleData(inflamed)
inflamed <- RunPCA(inflamed, reduction.name = "pca_inf")
#head(inflamed)
inflamed <- RunHarmony(inflamed, "batch", reduction = "pca_inf", reduction.save = "inf_harmony")
inflamed <- FindNeighbors(inflamed, reduction = "inf_harmony")
inflamed <- FindClusters(inflamed, resolution = 1.2)

inflamed <- RunUMAP(inflamed, reduction = "inf_harmony", reduction.name = "inf_harmony.umap", dims = 1:30)
DimPlot(inflamed, reduction = "inf_harmony.umap", label = TRUE)

#Lymphnode
DefaultAssay(Lymphnode) <- "RNA"
Lymphnode <- NormalizeData(Lymphnode)
Lymphnode <- FindVariableFeatures(Lymphnode)
Lymphnode <- ScaleData(Lymphnode)
Lymphnode <- RunPCA(Lymphnode, reduction.name = "pca_lymph")
#head(Lymphnode)
Lymphnode <- RunHarmony(Lymphnode, "batch", reduction = "pca_lymph", reduction.save = "lymph_harmony")
Lymphnode <- FindNeighbors(Lymphnode, reduction = "lymph_harmony")
Lymphnode <- FindClusters(Lymphnode, resolution = 0.6)

Lymphnode <- RunUMAP(Lymphnode, reduction = "lymph_harmony", reduction.name = "lymph_harmony.umap", dims = 1:30)
DimPlot(Lymphnode, reduction = "lymph_harmony.umap", label = TRUE)

# PBMC
DefaultAssay(PBMC) <- "RNA"
PBMC <- NormalizeData(PBMC)
PBMC <- FindVariableFeatures(PBMC)
PBMC <- ScaleData(PBMC)
PBMC <- RunPCA(PBMC, reduction.name = "pca_pbmc")
#head(PBMC)
PBMC <- RunHarmony(PBMC, "batch", reduction = "pca_pbmc", reduction.save = "pbmc_harmony")
PBMC <- FindNeighbors(PBMC, reduction = "pbmc_harmony")
PBMC <- FindClusters(PBMC, resolution = 0.5)

PBMC <- RunUMAP(PBMC, reduction = "pbmc_harmony", reduction.name = "pbmc_harmony.umap", dims = 1:30)

saveRDS(uninflamed, file = "uninflamed_only.rds")
saveRDS(inflamed, file = "inflamed_only.rds")
saveRDS(Lymphnode, file = "lymphnode_only.rds")
saveRDS(PBMC, file = "PBMC_only.rds")


markers_uninf <- FindAllMarkers(uninflamed, only.pos = F)
write.csv(markers_uninf, file = "uninf_RNA_markers.csv")

markers_inf <- FindAllMarkers(inflamed, only.pos = F)
write.csv(markers_inf, file = "inf_RNA_markers.csv")

markers_lymph <- FindAllMarkers(Lymphnode, only.pos = F)
write.csv(markers_lymph, file = "lymph_RNA_markers.csv")

markers_pbmc <- FindAllMarkers(PBMC, only.pos = F)
write.csv(markers_pbmc, file = "pbmc_RNA_markers.csv")

get_clones <- function(obj) na.omit(obj$clonotype_id)
uninflamed_clones <- get_clones(uninflamed)
inflamed_clones   <- get_clones(inflamed)
lymph_clones      <- get_clones(Lymphnode)
pbmc_clones       <- get_clones(PBMC)

clone_lists <- list(
  uninflamed = unique(uninflamed_clones),
  inflamed   = unique(inflamed_clones),
  Lymphnode  = unique(lymph_clones),
  PBMC       = unique(pbmc_clones)
)

all_clones <- unlist(clone_lists)
clone_table <- table(all_clones)

shared_clones <- names(clone_table[clone_table >= 2])

count_clones <- function(seurat_obj, clone_list) {
  tbl <- table(seurat_obj$clonotype_id)
  df <- data.frame(clonotype_id = names(tbl), count = as.integer(tbl), stringsAsFactors = FALSE)
  df <- subset(df, clonotype_id %in% clone_list)
  return(df)
}

uninflamed_counts <- count_clones(uninflamed, shared_clones)
inflamed_counts   <- count_clones(inflamed,   shared_clones)
lymph_counts      <- count_clones(Lymphnode,  shared_clones)
pbmc_counts       <- count_clones(PBMC,       shared_clones)

df <- Reduce(function(x, y) merge(x, y, by = "clonotype_id", all = TRUE),
             list(uninflamed_counts, inflamed_counts, lymph_counts, pbmc_counts))

colnames(df) <- c("clonotype_id", "uninflamed", "inflamed", "Lymphnode", "PBMC")

df[is.na(df)] <- 0

write.csv(df, file = "shared_clone.csv")
head(uninflamed)
uninflamed$annotation <- "Prolifrating T"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "0"] <- "TRM_2"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "1"] <- "TEM_1"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "2"] <- "eTreg2"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "3"] <- "Mixed_1"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "4"] <- "eTreg1"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "5"] <- "TRM_1"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "6"] <- "Tfh_1_1"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "7"] <- "Tfh_1_2"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "8"] <- "Mixed"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "9"] <- "naive T"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "10"] <- "Tfh_2"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "11"] <- "TEM_2"
uninflamed$annotation[uninflamed$RNA_snn_res.1 == "12"] <- "naive Treg"

head(inflamed)
inflamed$annotation <- "Prolifrating T"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "0"] <- "Tfh_1_1"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "1"] <- "TRM_1"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "2"] <- "eTreg_1"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "3"] <- "TEM_2"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "4"] <- "TRM_2"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "5"] <- "eTreg_1"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "6"] <- "naive Treg"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "7"] <- "naive T"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "8"] <- "eTreg_2_1"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "9"] <- "Tfh_1_2"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "10"] <- "Tfh_2_1"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "11"] <- "Tfh_2_2"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "12"] <- "Unannotated"
inflamed$annotation[inflamed$RNA_snn_res.1.2 == "13"] <- "eTreg2_2"
DimPlot(inflamed, reduction = "inf_harmony.umap", group.by = "annotation", label = T)

head(Lymphnode)
Lymphnode$annotation <- "Prolifrating T"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "0"] <- "Tfh_1"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "1"] <- "naive Treg"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "2"] <- "Tfh_2"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "3"] <- "naive T"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "4"] <- "eTreg"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "5"] <- "TCM"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "6"] <- "Teff_Th1"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "7"] <- "Unannotated"
Lymphnode$annotation[Lymphnode$RNA_snn_res.0.6 == "8"] <- "Teff_Th17"

PBMC$annotation <- "Prolifrating T"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "0"] <- "naive T"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "1"] <- "TEM_1"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "2"] <- "Unannotated_1"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "3"] <- "TCM_1"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "4"] <- "TEM_2"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "5"] <- "TCM_2"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "6"] <- "Treg"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "7"] <- "Teff_Th1"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "8"] <- "Teff_Th17"
PBMC$annotation[PBMC$RNA_snn_res.0.5 == "9"] <- "Unannotated_2"

saveRDS(uninflamed, file = "uninflamed_annotated.rds")
saveRDS(inflamed, file = "inflamed_annotated.rds")
saveRDS(Lymphnode, file = "Lymphnode_annotated.rds")
saveRDS(PBMC, file = "PBMC_annotated.rds")

uninflamed <- readRDS(file = "uninflamed_annotated.rds")
inflamed <- readRDS(file = "inflamed_annotated.rds")
Lymphnode <- readRDS(file = "Lymphnode_annotated.rds")
PBMC <- readRDS(file = "PBMC_annotated.rds")

# Fig4A

jpeg("uninflamed_UMAP.jpg", width = 2048, height = 2048, res = 400)
DimPlot(uninflamed, reduction = 'uninf_harmony.umap', group.by = "seurat_clusters", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

jpeg("inflamed_UMAP.jpg", width = 2048, height = 2048, res = 400)
DimPlot(inflamed, reduction = 'inf_harmony.umap', group.by = "seurat_clusters", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

jpeg("Lymphnode_UMAP.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Lymphnode, reduction = 'lymph_harmony.umap', group.by = "seurat_clusters", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

jpeg("PBMC_UMAP.jpg", width = 2048, height = 2048, res = 400)
DimPlot(PBMC, reduction = 'pbmc_harmony.umap', group.by = "seurat_clusters", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("sample0_clonotype6_scRep1",
                             "sample0_clonotype3_scRep1",
                             "sample0_clonotype4_scRep1",
                             "sample0_clonotype25_scRep1",
                             "sample0_clonotype2_scRep1",
                             "sample0_clonotype18_scRep1"
)

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = uninflamed, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "blue", "green", "magenta", "purple", "orange")

jpeg("uninf_TCR_1.jpg", width = 2048, height = 2048, res = 400)
DimPlot(uninflamed, reduction = "uninf_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("sample0_clonotype6_scRep1",
                             "sample0_clonotype25_scRep1",
                             "sample0_clonotype18_scRep1"
)

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = inflamed, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "magenta", "orange")

jpeg("inf_TCR_1.jpg", width = 2048, height = 2048, res = 400)
DimPlot(inflamed, reduction = "inf_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("sample0_clonotype6_scRep1",
                             "sample0_clonotype3_scRep1",
                             "sample0_clonotype4_scRep1",
                             "sample0_clonotype25_scRep1",
                             "sample0_clonotype2_scRep1",
                             "sample0_clonotype18_scRep1"
)

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = Lymphnode, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "blue", "green", "magenta", "purple", "orange")

jpeg("Lymph_TCR_1.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Lymphnode, reduction = "lymph_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("sample0_clonotype6_scRep1",
                             "sample0_clonotype3_scRep1",
                             "sample0_clonotype4_scRep1",
                             "sample0_clonotype2_scRep1"
)

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = PBMC, subset = clonotype_id == clonotype))
})


highlight_colors <- c("red", "blue", "green", "purple")

jpeg("pbmc_TCR_1.jpg", width = 2048, height = 2048, res = 400)
DimPlot(PBMC, reduction = "pbmc_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()


clonotypes_to_highlight <- c("20250220_clonotype17_scRep2")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = uninflamed, subset = clonotype_id == clonotype))
})

highlight_colors <- c("green")

jpeg("uninf_TCR_2.jpg", width = 2048, height = 2048, res = 400)
DimPlot(uninflamed, reduction = "uninf_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("20250220_clonotype5_scRep2",
                             "20250220_clonotype2_scRep2",
                             "20250220_clonotype4_scRep2",
                             "20250220_clonotype20_scRep2",
                             "20250220_clonotype13_scRep2")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = inflamed, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red","blue" , "magenta", "purple", "orange")

jpeg("inf_TCR_2.jpg", width = 2048, height = 2048, res = 400)
DimPlot(inflamed, reduction = "inf_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("20250220_clonotype5_scRep2",
                             "20250220_clonotype20_scRep2")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = Lymphnode, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "purple")

jpeg("Lymph_TCR_2.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Lymphnode, reduction = "lymph_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("20250220_clonotype5_scRep2",
                             "20250220_clonotype2_scRep2",
                             "20250220_clonotype17_scRep2",
                             "20250220_clonotype4_scRep2",
                             "20250220_clonotype20_scRep2",
                             "20250220_clonotype13_scRep2")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = PBMC, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red","blue" ,"green", "magenta", "purple", "orange") 

jpeg("pbmc_TCR_2.jpg", width = 2048, height = 2048, res = 400)
DimPlot(PBMC, reduction = "pbmc_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

clonotypes_to_highlight <- c("20250328_clonotype1_scRep3",
                             "20250328_clonotype2_scRep3",
                             "20250328_clonotype3_scRep3",
                             "20250328_clonotype4_scRep3",
                             "20250328_clonotype5_scRep3",
                             "20250328_clonotype7_scRep3")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = uninflamed, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "blue", "green", "magenta", "purple", "orange")

jpeg("uninf_TCR_3.jpg", width = 2048, height = 2048, res = 400)
DimPlot(uninflamed, reduction = "uninf_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()


clonotypes_to_highlight <- c("20250328_clonotype1_scRep3",
                             "20250328_clonotype2_scRep3",
                             "20250328_clonotype3_scRep3",
                             "20250328_clonotype4_scRep3",
                             "20250328_clonotype7_scRep3")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = inflamed, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "blue", "green", "magenta", "orange")


jpeg("inf_TCR_3.jpg", width = 2048, height = 2048, res = 400)
DimPlot(inflamed, reduction = "inf_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()


clonotypes_to_highlight <- c("20250328_clonotype1_scRep3",
                             "20250328_clonotype2_scRep3",
                             "20250328_clonotype5_scRep3")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = Lymphnode, subset = clonotype_id == clonotype))
})

highlight_colors <- c("red", "blue", "purple")

jpeg("Lymph_TCR_3.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Lymphnode, reduction = "lymph_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()


clonotypes_to_highlight <- c("20250328_clonotype3_scRep3",
                             "20250328_clonotype4_scRep3",
                             "20250328_clonotype5_scRep3",
                             "20250328_clonotype7_scRep3")

highlight_cell_list <- lapply(clonotypes_to_highlight, function(clonotype) {
  Cells(subset(x = PBMC, subset = clonotype_id == clonotype))
})

highlight_colors <- c("green", "magenta", "purple", "orange")

jpeg("pbmc_TCR_3.jpg", width = 2048, height = 2048, res = 400)
DimPlot(PBMC, reduction = "pbmc_harmony.umap",
        cells.highlight = highlight_cell_list,
        cols.highlight = highlight_colors,
        sizes.highlight = 2) +
  NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()


# Fig4B
get_clones <- function(obj) unique(na.omit(obj$clonotype_id))

clone_sets <- list(
  uninflamed = get_clones(uninflamed),
  inflamed   = get_clones(inflamed),
  Lymphnode  = get_clones(Lymphnode),
  PBMC       = get_clones(PBMC)
)

library(ggvenn)
ggvenn(clone_sets, 
       fill_color = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
       stroke_size = 0.5,
       set_name_size = 4)
ggsave("venn_clonotypes.pdf", width = 6, height = 6, dpi = 300)

# Fig4C
#which cluster share clone
count_clone_clusters <- function(seurat_obj, clone_list) {
  meta <- seurat_obj@meta.data
  df <- meta %>%
    filter(!is.na(clonotype_id), clonotype_id %in% clone_list) %>%
    select(clonotype_id, annotation) %>%
    group_by(clonotype_id, annotation) %>%
    summarise(cell_count = n(), .groups = "drop")
  return(df)
}

uninflamed_df <- count_clone_clusters(uninflamed, df$clonotype_id) %>% mutate(dataset = "uninflamed")
inflamed_df   <- count_clone_clusters(inflamed, df$clonotype_id) %>% mutate(dataset = "inflamed")
lymph_df      <- count_clone_clusters(Lymphnode, df$clonotype_id) %>% mutate(dataset = "Lymphnode")
pbmc_df       <- count_clone_clusters(PBMC, df$clonotype_id) %>% mutate(dataset = "PBMC")

combined_df <- bind_rows(uninflamed_df, inflamed_df, lymph_df, pbmc_df)

trm2_clones <- combined_df %>% 
  filter(annotation == "TRM_2") %>% 
  filter(dataset == "inflamed") %>% 
  pull(clonotype_id) %>% 
  unique()

shared_clones <- combined_df %>% 
  filter(clonotype_id %in% trm2_clones) %>% 
  group_by(dataset, annotation) %>% 
  summarise(
    unique_clones = n_distinct(clonotype_id),
    total_cells = sum(cell_count),
    .groups = "drop"
  ) %>% 
  arrange(dataset, annotation)

head(shared_clones)
write.csv(shared_clones, file = "TCR_shared_with_infTRM2_eachcluster.csv")

trm2_clones <- combined_df %>% 
  filter(annotation == "TRM_2") %>% 
  filter(dataset == "uninflamed") %>% 
  pull(clonotype_id) %>% 
  unique()

shared_clones <- combined_df %>% 
  filter(clonotype_id %in% trm2_clones) %>% 
  group_by(dataset, annotation) %>% 
  summarise(
    unique_clones = n_distinct(clonotype_id),
    total_cells = sum(cell_count),
    .groups = "drop"
  ) %>% 
  arrange(dataset, annotation)

head(shared_clones)
write.csv(shared_clones, file = "TCR_shared_with_uninfTRM2_eachcluster.csv")


