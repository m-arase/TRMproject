library(Seurat)
library(Signac)
library(hdf5r)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
set.seed(1234)

CD31 <- Read10X(data.dir ="/Users/sumitaninaoki/bioinformatics/CD/data/cellranger_v6/CD31/filtered_feature_bc_matrix")
Hy350 <- Read10X(data.dir ="/Users/sumitaninaoki/bioinformatics/CD/data/cellranger_v6/Hy350/filtered_feature_bc_matrix")
CK265 <- Read10X(data.dir ="/Users/sumitaninaoki/bioinformatics/CD/data/cellranger_v6/CK265/filtered_feature_bc_matrix")
CK278 <- Read10X(data.dir ="/Users/sumitaninaoki/bioinformatics/CD/data/cellranger_v6/CK278/filtered_feature_bc_matrix")
h5_file <- "/Users/sumitaninaoki/bioinformatics/CD/data/cellranger_v6/CK282/F5738_v783_loupe_and_matrix/filtered_feature_bc_matrix_CK282.h5"
CK282 <- Read10X_h5(filename = h5_file)
h5_file_2 <- "/Users/sumitaninaoki/bioinformatics/CD/data/cellranger_v6/Hy392/F6047_v819_loupe_and_matrix/filtered_feature_bc_matrix_Hy392.h5"
Hy392 <- Read10X_h5(filename = h5_file_2)

makeaseuratobject <- function(v1, v2, v3){
  tmp <- CreateSeuratObject(counts = v1[1])
  adt_assay <- CreateAssay5Object(counts = v1[2])
  tmp[["ADT"]] <- adt_assay
  tmp <- AddMetaData(tmp, metadata = v2, col.name="disease")
  tmp <- AddMetaData(tmp, metadata = v3, col.name="sample")#オブジェクト名をサンプル名としてメタデータに入力
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT.")
  tmp <- subset(tmp, subset = nFeature_RNA > 200  & percent.mt < 15)
  return(tmp)}

lst1 <- list(CD31, Hy350, CK265, CK278, CK282, Hy392)
lst2 <- list("CD","CD","CT","CT","CT","CD")
lst3 <- list("CD31","Hy350","CK265","CK278","CK282","Hy392")

CD31_obj <- makeaseuratobject(lst1[[1]], lst2[[1]], lst3[[1]])
Hy350_obj <- makeaseuratobject(lst1[[2]], lst2[[2]], lst3[[2]])
CK265_obj <- makeaseuratobject(lst1[[3]], lst2[[3]], lst3[[3]])
CK278_obj <- makeaseuratobject(lst1[[4]], lst2[[4]], lst3[[4]])
CK282_obj <- makeaseuratobject(lst1[[5]], lst2[[5]], lst3[[5]])
Hy392_obj <- makeaseuratobject(lst1[[6]], lst2[[6]], lst3[[6]])


obj <- merge(x = CD31_obj, y = c(Hy350_obj, CK265_obj, CK278_obj, CK282_obj, Hy392_obj), add.cell.ids = c("CD31", "Hy350", "CK265", "CK278", "CK282", "Hy392"))
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DefaultAssay(obj) <- "RNA"
obj <- IntegrateLayers(
  object = obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca"
)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "rpca_clusters0.5")
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

new.cluster.ids <- c("TEM_1", "eTreg_1", "Tfh_1", "TRM_1", "naïve T", "TRM_2", 
                     "TCM", "Tfh_2", "TEM_2", "naïve Treg", "eTreg_2",
                     "TEMRA_1", "TEMRA_2", "unannotated_1", "unannotated_2")

names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
obj$annotationnew <- Idents(obj)

new_levels <- c("naïve T","TEM_1","TEM_2", "TCM", "TRM_1","TRM_2","Tfh_1","Tfh_2","naïve Treg","eTreg_1","eTreg_2",
                "TEMRA_1","TEMRA_2","unannotated_1", "unannotated_2")
obj$annotationnew <- factor(x = obj$annotationnew, levels = new_levels)

cluster_colors4 <- c(
  "naïve T" = "#F8766D",  
  "TEM_1" = "#E58539", 
  "TEM_2" = "#B79F00",  
  "TCM" = "#FFFF00",  
  "TRM_1" = "#0000FF",  
  "TRM_2" = "#FF0000",  
  "Tfh_1" = "#00BF7D", 
  "Tfh_2" = "#6BD76B", 
  "naïve Treg" = "#00BFC4",  
  "eTreg_1" = "#00B0F6",
  "eTreg_2" = "#619CFF",  
  "TEMRA_1" = "#F564E3",  
  "TEMRA_2" = "#FF61C3",  
  "unannotated_1" = "#696969", 
  "unannotated_2" = "#A9A9A9"  
)

saveRDS(obj, "/Users/sumitaninaoki/bioinformatics/CD/data/seurat5/obj.rds")

# Fig1B
jpeg("CITE_UMAP_3_legend.jpg", width = 2048, height = 2048, res = 400)
DimPlot(obj, reduction = 'umap.rpca', group.by = "annotationnew", label = F, cols = cluster_colors4) + NoAxes() + ggtitle(NULL)
dev.off()

# Fig1C
# DotPlot
objsub <- subset(x = obj, idents=c("naïve T","TEM_1","TEM_2", "TCM", "TRM_1","TRM_2","Tfh_1","Tfh_2","naïve Treg","eTreg_1","eTreg_2",
                                   "TEMRA_1","TEMRA_2"))
new_levels <- c("naïve T","TEM_1","TEM_2", "TCM", "TRM_1","TRM_2","Tfh_1","Tfh_2","naïve Treg","eTreg_1","eTreg_2",
                "TEMRA_1","TEMRA_2")
new_levels <- rev(new_levels)
objsub$annotationnew <- factor(x = objsub$annotationnew, levels = new_levels)
DefaultAssay(objsub) <- "RNA"
dotfeaturesRNA <- c("CCR7", "TCF7", "SELL", "CD27", "CCR5", "KLRB1", "CXCR5", "CXCL13", "ICOS", "PDCD1", 
                    "BCL6", "TBX21", "GATA3", "RORC", "IL7R", "CD69", "ITGAE","S1PR1", "KLF2",
                    "FOXP3", "IKZF4", "CTLA4", "NKG7", "IFNG", "GZMA", "GZMB", "PRF1")
pdf("CITE_RNA_Dot2.pdf", width = 15, height = 6)
DotPlot(objsub, features = dotfeaturesRNA, dot.scale=10, group.by = "annotationnew", cols="RdBu") + 
  RotatedAxis()
dev.off()

# Fig1D
# DotPlot
dotfeaturesADT <- c("CD45RA", "CD45RO", "CD103", "CD69.1", "CD185", "CD62L", "CD25")
pdf("CITE_ADT_Dot.pdf", width = 6, height = 6)
DotPlot(objsub, features = dotfeaturesADT, dot.scale=10, group.by = "annotationnew", cols="RdBu") + 
  RotatedAxis()
dev.off()

# Fig1E
write.csv(prop.table(table(obj$disease, obj$annotationnew), margin = 1), file = "CITE_celltype_per_disease_prop_2.csv")

# Fig1F
pdf("CITE_UMAP_disease_nolegend.pdf", width = 10, height = 10)
DimPlot(obj, reduction = 'umap.rpca', group.by = "disease", label = F) + NoLegend()
dev.off()

# Fig1G
library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(ggplot2)

obj <- readRDS(file = "/Users/sumitaninaoki/bioinformatics/CD/data/seurat5/obj2_rename.rds")
obj <- JoinLayers(obj)
DefaultAssay(obj) <- "RNA"

obj[["RNA2"]] <- as(object = obj[["RNA"]], Class = "Assay")
DefaultAssay(obj) <- "RNA2"
obj_milo <- as.SingleCellExperiment(obj, assay = "RNA2")
reducedDim(obj_milo, "UMAP", withDimnames=TRUE) <- obj[['umap.cca']]@cell.embeddings

obj_milo <- Milo(obj_milo)
obj_milo <- buildGraph(obj_milo, k = 30, d = 30)
obj_milo <- makeNhoods(obj_milo, prop = 0.1, k = 30, d=30, refined = TRUE)
plotNhoodSizeHist(obj_milo)
obj_milo <- countCells(obj_milo, meta.data = data.frame(colData(obj_milo)), samples="sample")
traj_design <- data.frame(colData(obj_milo))[,c("sample", "disease")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
traj_design <- traj_design[colnames(nhoodCounts(obj_milo)), , drop=FALSE]
traj_design

obj_milo <- calcNhoodDistance(obj_milo, d=30)
class(obj_milo)
rownames(traj_design) <- traj_design$sample
da_results <- testNhoods(obj_milo, design = ~ disease, design.df = traj_design)
#da_results <- testNhoods(obj_milo, design = ~ sample, design.df = traj_design)
da_results %>% arrange(- SpatialFDR) %>% head()
da_results <- annotateNhoods(obj_milo, da_results, coldata_col = "ident")
head(da_results)
pdf("/Users/sumitaninaoki/bioinformatics/CD/data/obj_milo_beeswarm.pdf", width = 10, height = 7)
plotDAbeeswarm(da_results, group.by = "ident")
dev.off()

obj_milo <- buildNhoodGraph(obj_milo)
pdf("/Users/sumitaninaoki/bioinformatics/CD/data/obj2_milo.pdf", width = 10, height = 7)
plotUMAP(obj_milo) + plotNhoodGraphDA(obj_milo, da_results, alpha=0.1) +
  plot_layout(guides="collect")
dev.off()

nh_graph_pl <- plotNhoodGraphDA(obj_milo, da_results, layout="UMAP",alpha=0.1)
ggsave(file = "/Users/sumitaninaoki/bioinformatics/CD/data/obj2_milo.png", plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 300)

saveRDS(obj_milo, file = "/Users/sumitaninaoki/bioinformatics/CD/data/seurat5/obj2_milo.rds")

# Fig1H
objsub <- subset(x = obj, idents=c("TEM_1", "Treg_1", "Tfh_1", "TRM_1", "naïve T", "TRM_2", 
                                   "TEM_2", "Tfh_2", "TEM_3", "Treg_2", "Treg_3",
                                   "TEMRA_1", "TEMRA_2"))
objsub$celltype <- "TEMRA"
objsub$celltype[objsub$annotation == "TEM_1"] <- "TEM"
objsub$celltype[objsub$annotation == "TEM_2"] <- "TCM"
objsub$celltype[objsub$annotation == "TEM_3"] <- "TEM"

objsub$celltype[objsub$annotation == "Treg_1"] <- "Treg"
objsub$celltype[objsub$annotation == "Treg_2"] <- "Treg"
objsub$celltype[objsub$annotation == "Treg_3"] <- "Treg"

objsub$celltype[objsub$annotation == "Tfh_1"] <- "Tfh"
objsub$celltype[objsub$annotation == "Tfh_2"] <- "Tfh"

objsub$celltype[objsub$annotation == "TRM_1"] <- "TRM_1"
objsub$celltype[objsub$annotation == "TRM_2"] <- "TRM_2"

objsub$celltype[objsub$annotation == "naïve T"] <- "naïve T"

pdf("CITE_RNA_Vln_ITGAE.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "RNA"
VlnPlot(objsub, features = "ITGAE", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_RNA_Vln_S1PR1.pdf", width = 4, height = 2.5)
DefaultAssay(obj) <- "RNA"
VlnPlot(obj, features = "S1PR1", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_RNA_Vln_IFNG.pdf", width = 4, height = 2.5)
DefaultAssay(obj) <- "RNA"
VlnPlot(obj, features = "IFNG", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_RNA_Vln_CD69.pdf", width = 4, height = 2.5)
DefaultAssay(obj) <- "RNA"
VlnPlot(obj, features = "CD69", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_RNA_Vln_GZMB.pdf", width = 4, height = 2.5)
DefaultAssay(obj) <- "RNA"
VlnPlot(obj, features = "GZMB", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_RNA_Vln_PRF1.pdf", width = 4, height = 2.5)
DefaultAssay(obj) <- "RNA"
VlnPlot(obj, features = "PRF1", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_RNA_Vln_GZMH.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "RNA"
VlnPlot(objsub, features = "GZMH", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_IFNG_2.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("IFNG"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_GZMB.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("GZMB"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_PRF1.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("PRF1"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_ITGAE.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("ITGAE"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_KLF2.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("KLF2"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_S1PR1.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("S1PR1"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_GZMA.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("GZMA"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

pdf("CITE_GZMH.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("GZMH"), 
            reduction = "umap.rpca") + NoAxes()
dev.off()

# Fig1I
pdf("CITE_ADT_Vln_CD103.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "ADT"
VlnPlot(objsub, features = "CD103", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_ADT_Vln_CD38.1.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "ADT"
VlnPlot(objsub, features = "CD38.1", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_ADT_Vln_HLA-DR.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "ADT"
VlnPlot(objsub, features = "HLA-DR", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_ADT_Vln_CD39.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "ADT"
VlnPlot(objsub, features = "CD39", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_ADT_Vln_CD29.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "ADT"
VlnPlot(objsub, features = "CD29", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

pdf("CITE_ADT_Vln_CD101.pdf", width = 4, height = 2.5)
DefaultAssay(objsub) <- "ADT"
VlnPlot(objsub, features = "CD101", 
        group.by = "celltype", idents = c("TRM_1", "TRM_2", "TEM", "Treg", "Tfh", "naïve T"),
        pt.size = 0)
dev.off()

DefaultAssay(obj) <- "ADT"
pdf("CITE_CD103_adt_2.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("CD103"), 
            cols = c("lightgrey", "darkgreen"),
            reduction = "umap.rpca",
            max.cutoff = 3) + NoAxes()
dev.off()

DefaultAssay(obj) <- "ADT"
pdf("CITE_CD38_adt.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("CD38.1"), 
            cols = c("lightgrey", "darkgreen"),
            reduction = "umap.rpca") + NoAxes()
dev.off()

DefaultAssay(obj) <- "ADT"
pdf("CITE_CD39_adt.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("CD39"), 
            cols = c("lightgrey", "darkgreen"),
            reduction = "umap.rpca") + NoAxes()
dev.off()

DefaultAssay(obj) <- "ADT"
pdf("CITE_HLADR_adt.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("HLA-DR"), 
            cols = c("lightgrey", "darkgreen"),
            reduction = "umap.rpca") + NoAxes()
dev.off()

DefaultAssay(obj) <- "ADT"
pdf("CITE_CD29_adt.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("CD29"), 
            cols = c("lightgrey", "darkgreen"),
            reduction = "umap.rpca"
) + NoAxes()
dev.off()

DefaultAssay(obj) <- "ADT"
pdf("CITE_CD101_adt.pdf", width = 3, height = 3)
FeaturePlot(obj, 
            features = c("CD101"), 
            cols = c("lightgrey", "darkgreen"),
            reduction = "umap.rpca"
) + NoAxes()
dev.off()






