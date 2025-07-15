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

# Multiome data from Figure5
# FigS3B
DefaultAssay(Multiome_data) <- "prediction.score.celltype.annotation"

pdf("Multiome_UMAP_predictionscore_naive.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("naïve T"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TEM_1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TEM-1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TCM.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TCM"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TEM_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TEM-2"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TRM_1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TRM-1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TRM_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TRM-2"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_Tfh_1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("Tfh-1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_Tfh_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("Tfh-2"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_eTreg_1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("eTreg-1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_naive_Treg_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("naïve Treg"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_eTreg_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("eTreg-2"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TEMRA_1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TEMRA-1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_TEMRA_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("TEMRA-2"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_unannotated_1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("unannotated-1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_predictionscore_unannotated_2.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("unannotated-2"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

# FigS3C
jpeg("Multiome_UMAP_sample.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Multiome_data, reduction = "wnn.umap", group.by = "sample", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

write.csv(table(Multiome_data$disease, Multiome_data$celltype), file = "Multiome_proportion_disease.csv")
write.csv(table(Multiome_data$sample, Multiome_data$celltype), file = "Multiome_proportion_sample.csv")

# FigS3D
jpeg("Multiome_UMAP6.jpg", width = 2048, height = 2048, res = 400)
DimPlot(Multiome_data, reduction = "wnn.umap", group.by = "annotation", label = F) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# FigS3E
DefaultAssay(Multiome_data) <- "RNA"
dotfeaturesRNA <- c("CCR7", "TCF7", "SELL", "CD27", "CCR5", "KLRB1", "CXCR5", "CXCL13", "ICOS", "PDCD1", 
                    "BCL6", "TBX21", "GATA3", "RORC", "IL7R", "CD69", "ITGAE","S1PR1", "KLF2",
                    "FOXP3", "IKZF4", "CTLA4", "NKG7", "IFNG", "GZMA", "GZMB", "PRF1")
pdf("Multiome_RNA_Dot.pdf", width = 12, height = 6)
DotPlot(Multiome_data, features = dotfeaturesRNA, dot.scale=10, group.by = "celltype", cols="RdBu") + 
  RotatedAxis()
dev.off()

# FigS3H
pdf("Multiome_UMAP_RUNX1.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("RUNX1"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

pdf("Multiome_UMAP_RUNX3.pdf", width = 3, height = 3)
FeaturePlot(Multiome_data, 
            features = c("RUNX3"), 
            reduction = "wnn.umap") + NoAxes()
dev.off()

# FigS3I
library(decoupleR)

net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
# Extract the normalized log-transformed counts
mat <- as.matrix(CD103LP@assays$RNA@data)

# Run wmean
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)

CD103LP[['tfswmean']] <- acts %>%
  dplyr::filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

DefaultAssay(object = CD103LP) <- "tfswmean"
CD103LP <- ScaleData(CD103LP)
CD103LP@assays$tfswmean@data <- CD103LP@assays$tfswmean@scale.data

CD103LP3 <- subset(x = CD103LP, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
new_names <- c("T_CD_TXNIP", "T_CK_HSP_1", "T_CK_FOSB", "T_CD_CCL20",
               "T_CD_Cytotoxic", "MT_high", "T_CK_HSP_2", "Treg", "NK_like", "naive")
names(new_names) <- levels(CD103LP3)
CD103LP3 <- RenameIdents(CD103LP3, new_names)
CD103LP3$annotation <- Idents(CD103LP3)

pdf("CD103_decoupleR_RUNX2.pdf", width = 3, height = 3)
FeaturePlot(CD103LP3, features = c("RUNX2"), reduction = "harmony.umap") + NoAxes() & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
dev.off()

pdf("CD103_decoupleR_BHLHE40.pdf", width = 3, height = 3)
FeaturePlot(CD103LP3, features = c("BHLHE40"), reduction = "harmony.umap") + NoAxes() & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
dev.off()








