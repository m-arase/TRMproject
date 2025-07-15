library(Seurat)
library(SeuratObject)
library(tidyverse)
library(scRepertoire)
library(circlize)
library(scales)
library(ggraph)
set.seed(1234)

# FigS2A
obj2 <- readRDS(file = "/Volumes/BUFFALOHDD/CITEseqAnalysis/20240724_Figure/CITE20240731_repertoire.rds")
new.cluster.ids <- c("TEM_1", "eTreg_1", "Tfh_1", "TRM_1", "naïve T", "TRM_2", 
                     "TCM", "Tfh_2", "TEM_2", "naïve Treg", "eTreg_2",
                     "TEMRA_1", "TEMRA_2")

names(new.cluster.ids) <- levels(obj2)
obj2 <- RenameIdents(obj2, new.cluster.ids)
obj2$annotationnew <- Idents(obj2)

new_levels <- c("naïve T","TEM_1","TEM_2", "TCM", "TRM_1","TRM_2","Tfh_1","Tfh_2","naïve Treg","eTreg_1","eTreg_2",
                "TEMRA_1","TEMRA_2","unannotated_1", "unannotated_2")

obj2$annotationnew <- factor(x = obj2$annotationnew, levels = new_levels)

CD31_obj <- subset(x = obj2, subset = sample == "CD31")
Hy350_obj <- subset(x = obj2, subset = sample == "Hy350")
Hy392_obj <- subset(x = obj2, subset = sample == "Hy392")

CD31_allTCR <- CD31_obj[[c("annotationnew", "clonotype_id")]]
CD31_allTCR <- CD31_allTCR %>% filter(!is.na(clonotype_id))

write.csv(table(CD31_allTCR$clonotype_id), file = "CD31_clone.csv")

Hy350_allTCR <- Hy350_obj[[c("annotationnew", "clonotype_id")]]
Hy350_allTCR <- Hy350_allTCR %>% filter(!is.na(clonotype_id))

write.csv(table(Hy350_allTCR$clonotype_id), file = "Hy350_clone.csv")

Hy392_allTCR <- Hy392_obj[[c("annotationnew", "clonotype_id")]]
Hy392_allTCR <- Hy392_allTCR %>% filter(!is.na(clonotype_id))

write.csv(table(Hy392_allTCR$clonotype_id), file = "Hy392_clone.csv")

CD31_trm2_data <- CD31_allTCR[CD31_allTCR$annotationnew == "TRM_2", ]
Hy350_trm2_data <- Hy350_allTCR[Hy350_allTCR$annotationnew == "TRM_2", ]
Hy392_trm2_data <- Hy392_allTCR[Hy392_allTCR$annotationnew == "TRM_2", ]

write.csv(table(CD31_trm2_data $clonotype_id), file = "CD31_TRM2_clone.csv")
write.csv(table(Hy350_trm2_data$clonotype_id), file = "Hy350_TRM2_clone.csv")
write.csv(table(Hy392_trm2_data$clonotype_id), file = "Hy392_TRM2_clone.csv")

# FigS2B
p1 <- vizGenes(obj_TRM_2_CD31, 
               x.axis = "TRAV",
               y.axis = NULL,
               plot = "barplot",  
               scale = FALSE) + 
  theme(axis.text = element_text(size = 20)
  )
data1 <- p1$data
data1$y.axis <- NULL
data1$sd <- NULL

p2 <- vizGenes(obj_TRM_2_Hy350, 
               x.axis = "TRAV",
               y.axis = NULL,
               plot = "barplot",  
               scale = FALSE) + 
  theme(axis.text = element_text(size = 20)
  )
data2 <- p2$data
data2$y.axis <- NULL
data2$sd <- NULL

p3 <- vizGenes(obj_TRM_2_Hy392, 
               x.axis = "TRAV",
               y.axis = NULL,
               plot = "barplot",  
               scale = FALSE) + 
  theme(axis.text = element_text(size = 20)
  )
data3 <- p3$data
data3$y.axis <- NULL
data3$sd <- NULL

head(data1)
TRAVdata <- merge(data1, data2, by = "x.axis", all = TRUE)
TRAVdata <- merge(TRAVdata, data3, by = "x.axis", all = TRUE)
colnames(TRAVdata) <- c("TRAV", "CD31", "Hy350", "Hy392")
write.csv(TRAVdata, file = "CD_TRM_TRAV.csv")

p4 <- vizGenes(obj_TRM_2_CD31, 
               x.axis = "TRBV",
               y.axis = NULL,
               plot = "barplot",  
               scale = FALSE) + 
  theme(axis.text = element_text(size = 20)
  )
data4 <- p4$data
data4$y.axis <- NULL
data4$sd <- NULL

p5 <- vizGenes(obj_TRM_2_Hy350, 
               x.axis = "TRBV",
               y.axis = NULL,
               plot = "barplot",  
               scale = FALSE) + 
  theme(axis.text = element_text(size = 20)
  )
data5 <- p5$data
data5$y.axis <- NULL
data5$sd <- NULL

p6 <- vizGenes(obj_TRM_2_Hy392, 
               x.axis = "TRBV",
               y.axis = NULL,
               plot = "barplot",  
               scale = FALSE) + 
  theme(axis.text = element_text(size = 20)
  )
data6 <- p6$data
data6$y.axis <- NULL
data6$sd <- NULL

TRBVdata <- merge(data1, data2, by = "x.axis", all = TRUE)
TRBVdata <- merge(TRBVdata, data3, by = "x.axis", all = TRUE)
colnames(TRBVdata) <- c("TRBV", "CD31", "Hy350", "Hy392")
write.csv(TRBVdata, file = "CD_TRM_TRBV.csv")

# FigS2C
# CD103LP2 from Figure2
CD103LP2$clonotype_id_2 <- paste0(CD103LP2$sample, "_", CD103LP2$clonotype_id)
head(CD103LP2)

library(scRepertoire)
batch1TCR <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220215_CD103_hashtag/F3804_220215_152314_cellranger/LPMC_TCR/outs/filtered_contig_annotations.csv")
batch2TCR <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220324_CD103_hashtag/F3970_220323_162237_cellranger/CD103_TCR/outs/filtered_contig_annotations.csv")

batch1 <- subset(CD103LP2, subset = batch == "batch1")
batch2 <- subset(CD103LP2, subset = batch == "batch2")
head(batch1)
head(batch1TCR)
batch1TCR$barcode <- paste0("hash1_", batch1TCR$barcode)
batch2TCR$barcode <- paste0("hash2_", batch2TCR$barcode)

contig.list1 <- createHTOContigList(batch1TCR, 
                                    batch1,
                                    group.by = "sample")
contig.list2 <- createHTOContigList(batch2TCR,
                                    batch2,
                                    group.by = "sample")

combined.TCR1 <- combineTCR(contig.list1,
                            samples = c("CK206", "CD259", "CD253", "CK231", "CK238", "CD268"),
                            removeNA = FALSE,
                            removeMulti = FALSE,
                            filterMulti = FALSE)

combined.TCR2 <- combineTCR(contig.list2,
                            samples = c("CK228", "CK221", "CD24", "CD189", "CK242", "CD27"),
                            removeNA = FALSE,
                            removeMulti = FALSE,
                            filterMulti = FALSE)
for(i in 1:6){
  combined.TCR1[[i]]$barcode <- sub(".*hash1", "hash1", combined.TCR1[[i]]$barcode)
}

for(i in 1:6){
  combined.TCR2[[i]]$barcode <- sub(".*hash2", "hash2", combined.TCR2[[i]]$barcode)
}

batch1 <- combineExpression(combined.TCR1, 
                            batch1, 
                            cloneCall="gene", 
                            group.by = "sample", 
                            proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
head(batch1)

batch2 <- combineExpression(combined.TCR2, 
                            batch2, 
                            cloneCall="gene", 
                            group.by = "sample", 
                            proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

CD103LPTCR <- merge(batch1, y = batch2)
CD103LPTCR
DimPlot(CD103LPTCR, reduction = "harmony.umap", label = T, group.by = "annotation")
head(CD103LPTCR)

CD103LP2@meta.data <- CD103LPTCR@meta.data

saveRDS(CD103LP2, file = "CD103_scRepertiore.rds")

plot <- DimPlot(CD103LP2, reduction = "harmony.umap")
select.cells <- CellSelector(plot = plot)
CD103LP2 <- CD103LP2[,!colnames(CD103LP2) %in% select.cells]

DimPlot(CD103LP2, reduction = "harmony.umap")

CD103LP2_CD <- subset(x = CD103LP2, subset = disease == "CD")
CD103LP2_CK <- subset(x = CD103LP2, subset = disease == "CK")

pdf("CD103_CD_TCR_clone.pdf", width = 6, height = 6)
Seurat::DimPlot(CD103LP2_CD, group.by = "cloneSize", reduction = "harmony.umap", pt.size = 1.5) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

pdf("CD103_CK_TCR_clone.pdf", width = 6, height = 6)
Seurat::DimPlot(CD103LP2_CK, group.by = "cloneSize", reduction = "harmony.umap", pt.size = 1.5) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# FigS2D
CD103LP2_CD <- subset(x = CD103LP2, subset = disease == "CD")
CD103LP2_CK <- subset(x = CD103LP2, subset = disease == "CK")

write.csv(table(CD103LP2_CD$clonotype_id), file = "CD103_CD_allclonotype.csv")
write.csv(table(CD103LP2_CK$clonotype_id), file = "CD103_CK_allclonotype.csv")

# FigS2E
# TCR circle
pdf("CD103_CD_TCR_circle_new.pdf", width = 5, height = 5)
makeTCRcircle(CD103LP2_CD)
dev.off()

pdf("CD103_CK_TCR_circle_new.pdf", width = 5, height = 5)
makeTCRcircle(CD103LP2_CK)
dev.off()

# CD UpsetR
CD103_CD_allTCR <- CD103LP2_CD[[c("annotationnew", "clonotype_id")]]
CD103_CD_allTCR <- CD103_CD_allTCR %>% filter(!is.na(clonotype_id))
resultCD <- CD103_CD_allTCR %>%
  group_by(annotationnew) %>%
  summarize(clonotype_ids = list(clonotype_id))

result_listCD <- setNames(as.list(resultCD$clonotype_ids), resultCD$annotationnew)

all_elements <- sort(unique(unlist(result_listCD)))

binary_matrix <- matrix(0, nrow = length(result_listCD), ncol = length(all_elements))
rownames(binary_matrix) <- names(result_listCD)
colnames(binary_matrix) <- all_elements

for (i in 1:length(result_listCD)) {
  binary_matrix[i, all_elements %in% result_listCD[[i]]] <- 1
}
binary_matrix <- t(binary_matrix)



sets <- c("TXNIP high", "HSP high_1", "FOSB high", "CCL20 high",
          "Cytotoxic", "MT high", "HSP high_2", "Treg", "NK_like", "naive")
reduced_data <- binary_matrix[rowSums(binary_matrix[, sets]) > 1, ]

reduced_data <- data.frame(reduced_data)
sets <- c("TXNIP.high", "HSP.high_1", "FOSB.high", "CCL20.high",
          "Cytotoxic", "MT.high", "HSP.high_2", "Treg", "NK_like", "naive")
sets <- rev(sets)

upset(reduced_data , sets = sets, order.by = "freq", nintersects = NA, keep.order = TRUE)

pdf("CD103_CD_TCR_upset_new.pdf", width = 9, height = 9)
upset(reduced_data ,sets = sets, order.by = "freq", nintersects = NA, keep.order = TRUE)
dev.off()

pdf("CD103_CD_TCR_upset_new_legend.pdf", width = 9, height = 9)
sets <- c("TXNIP high", "HSP high_1", "FOSB high", "CCL20 high",
          "Cytotoxic", "MT high", "HSP high_2", "Treg", "NK_like", "naive")
sets <- rev(sets)
upset(fromList(result_listCD), sets = sets, nintersects = NA, order.by = "freq", keep.order = TRUE)
dev.off()

# CK UpsetR
CD103_CK_allTCR <- CD103LP2_CK[[c("annotationnew", "clonotype_id")]]
CD103_CK_allTCR <- CD103_CK_allTCR %>% filter(!is.na(clonotype_id))
resultCK <- CD103_CK_allTCR %>%
  group_by(annotationnew) %>%
  summarize(clonotype_ids = list(clonotype_id))

result_listCK <- setNames(as.list(resultCK$clonotype_ids), resultCK$annotationnew)

all_elements <- sort(unique(unlist(result_listCK)))

binary_matrix <- matrix(0, nrow = length(result_listCK), ncol = length(all_elements))
rownames(binary_matrix) <- names(result_listCK)
colnames(binary_matrix) <- all_elements

for (i in 1:length(result_listCK)) {
  binary_matrix[i, all_elements %in% result_listCK[[i]]] <- 1
}
binary_matrix <- t(binary_matrix)


sets <- c("TXNIP high", "HSP high_1", "FOSB high", "CCL20 high",
          "Cytotoxic", "MT high", "HSP high_2", "Treg", "NK_like", "naive")
reduced_data <- binary_matrix[rowSums(binary_matrix[, sets]) > 1, ]
reduced_data <- data.frame(reduced_data)
sets <- c("TXNIP.high", "HSP.high_1", "FOSB.high", "CCL20.high",
          "Cytotoxic", "MT.high", "HSP.high_2", "Treg", "NK_like", "naive")
sets <- rev(sets)


pdf("CD103_CK_TCR_upset_new.pdf", width = 9, height = 9)
upset(reduced_data ,sets = sets, order.by = "freq", nintersects = NA, keep.order = TRUE)
dev.off()

pdf("CD103_CK_TCR_upset_new_legend.pdf", width = 9, height = 9)
sets <- c("TXNIP high", "HSP high_1", "FOSB high", "CCL20 high",
          "Cytotoxic", "MT high", "HSP high_2", "Treg", "NK_like", "naive")
sets <- rev(sets)
upset(fromList(result_listCK), sets = sets, nintersects = 15, order.by = "freq", keep.order = TRUE)
dev.off()





