library(Seurat)
library(SeuratObject)
library(tidyverse)
library(scRepertoire)
library(circlize)
library(scales)
library(ggraph)
library(UpSetR)
set.seed(1234)

# CD31 : CD1
# Hy350 : CD2
# Hy392 : CD3
# CK265 : Ctrl1
# CK278 : Ctrl2
# CK282 : Ctrl3

obj <- readRDS(file = "/Volumes/BUFFALOHDD/CITEseqAnalysis/20240611_referencemapping/CITE_rpca_annotation.rds")

makeTCRdf <- function(tcrpath, clonopath, samplename){
  TCR <- read.csv(tcrpath)
  clono <- read.csv(clonopath)
  TCR <- TCR[!duplicated(TCR$barcode), ]
  TCR <- TCR[,c("barcode", "raw_clonotype_id")]
  names(TCR)[names(TCR) == "raw_clonotype_id"] <- "clonotype_id"
  TCR <- merge(TCR, clono[, c("clonotype_id", "cdr3s_aa")])
  TCR <- TCR[, c(2,1,3)]
  TCR[,1] <- paste0(samplename, "_", TCR[,1])
  TCR[["clonotype_id"]] <- paste0(TCR[['clonotype_id']], '_', samplename)
  return(TCR)
}

tcrCD31 <- "/Volumes/BUFFALOHDD/scRNAseq/20220526_CITE-seq/F4241_220525_152545_cellranger/CD31_CITE_TCR/outs/filtered_contig_annotations.csv"
clonoCD31 <- "/Volumes/BUFFALOHDD/scRNAseq/20220526_CITE-seq/F4241_220525_152545_cellranger/CD31_CITE_TCR/outs/clonotypes.csv"
CD31TCR <- makeTCRdf(tcrCD31, clonoCD31, "CD31")
tcrCK265 <- "/Volumes/BUFFALOHDD/scRNAseq/20221027_CK265_CITE/F4852_221026_122440_cellranger/CK265_TCR/outs/filtered_contig_annotations.csv"
clonoCK265 <- "/Volumes/BUFFALOHDD/scRNAseq/20221027_CK265_CITE/F4852_221026_122440_cellranger/CK265_TCR/outs/clonotypes.csv"
CK265TCR <- makeTCRdf(tcrCK265, clonoCK265, "CK265")
tcrCK278 <- "/Volumes/BUFFALOHDD/scRNAseq/20230220_CK278_CITE/F5438_230220_114708_cellranger/CK278_TCR/outs/filtered_contig_annotations.csv"
clonoCK278 <- "/Volumes/BUFFALOHDD/scRNAseq/20230220_CK278_CITE/F5438_230220_114708_cellranger/CK278_TCR/outs/clonotypes.csv"
CK278TCR <- makeTCRdf(tcrCK278, clonoCK278, "CK278")
tcrHy350 <- "/Volumes/BUFFALOHDD/scRNAseq/20230306_Hy350_CITE/F5489_230306_124037_cellranger/Hy350_TCR/outs/filtered_contig_annotations.csv"
clonoHy350 <- "/Volumes/BUFFALOHDD/scRNAseq/20230306_Hy350_CITE/F5489_230306_124037_cellranger/Hy350_TCR/outs/clonotypes.csv"
Hy350TCR <- makeTCRdf(tcrHy350, clonoHy350, "Hy350")
tcrCK282 <- "/Volumes/BUFFALOHDD/scRNAseq/20230414_CK282_CITE/CK282/outs/per_sample_outs/CK282/vdj_t/filtered_contig_annotations.csv"
clonoCK282 <- "/Volumes/BUFFALOHDD/scRNAseq/20230414_CK282_CITE/CK282/outs/per_sample_outs/CK282/vdj_t/clonotypes.csv"
CK282TCR <- makeTCRdf(tcrCK282, clonoCK282, "CK282")
tcrHy392 <- "/Volumes/BUFFALOHDD/scRNAseq/20230620_Hy392_CITE/Hy392/outs/per_sample_outs/Hy392/vdj_t/filtered_contig_annotations.csv"
clonoHy392 <- "/Volumes/BUFFALOHDD/scRNAseq/20230620_Hy392_CITE/Hy392/outs/per_sample_outs/Hy392/vdj_t/clonotypes.csv"
Hy392TCR <- makeTCRdf(tcrHy392, clonoHy392, "Hy392")

head(CD31TCR)
All_TCR <- dplyr::bind_rows(CD31TCR, CK265TCR, CK278TCR, Hy350TCR, CK282TCR, Hy392TCR)
head(All_TCR)
TCRlist <- left_join(celllist, All_TCR, by = c('celllist'='barcode'))
head(TCRlist)
rownames(TCRlist) <- TCRlist[,1]
TCRlist[,1] <- NULL

obj <- AddMetaData(object=obj, metadata=TCRlist)

library(scRepertoire)
obj <- readRDS(file = "/Volumes/BUFFALOHDD/CITEseqAnalysis/20240724_Figure/CITEwithTCR.rds")
CD31 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20220526_CITE-seq/F4241_220525_152545_cellranger/CD31_CITE_TCR/outs/filtered_contig_annotations.csv")
Hy350 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20230306_Hy350_CITE/F5489_230306_124037_cellranger/Hy350_TCR/outs/filtered_contig_annotations.csv")
Hy392 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20230620_Hy392_CITE/Hy392/outs/per_sample_outs/Hy392/vdj_t/filtered_contig_annotations.csv")
CK265 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20221027_CK265_CITE/F4852_221026_122440_cellranger/CK265_TCR/outs/filtered_contig_annotations.csv")
CK278 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20230220_CK278_CITE/F5438_230220_114708_cellranger/CK278_TCR/outs/filtered_contig_annotations.csv")
CK282 <- read.csv("/Volumes/BUFFALOHDD/scRNAseq/20230414_CK282_CITE/CK282/outs/per_sample_outs/CK282/vdj_t/filtered_contig_annotations.csv")

contig_list <- list(CD31, Hy350, Hy392, CK265, CK278, CK282)
combined.TCR <- combineTCR(contig_list, 
                           samples = c("CD31", "Hy350", "Hy392", 
                                       "CK265","CK278", "CK282"),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

scRep <- combineExpression(combined.TCR, 
                           obj, 
                           cloneCall="gene", 
                           group.by = "sample", 
                           proportion = FALSE, 
                           cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

obj <- scRep
CDobj <- subset(x = obj, subset = disease == "CD")
CKobj <- subset(x = obj, subset = disease == "CT")

# Fig3A
pdf("CITE_CD_TCR_clone.pdf", width = 6, height = 6)
Seurat::DimPlot(CDobj, group.by = "cloneSize", reduction = "umap.rpca", pt.size = 1.5) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

pdf("CITE_CK_TCR_clone.pdf", width = 4.5, height = 6)
Seurat::DimPlot(CKobj, group.by = "cloneSize", reduction = "umap.rpca", pt.size = 1.5) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) + NoLegend() + NoAxes() + ggtitle(NULL)
dev.off()

# Fig3B
write.csv(table(CDobj$clonotype_id), file = "CITE_CD_allclonotype.csv")
write.csv(table(CKobj$clonotype_id), file = "CITE_CK_allclonotype.csv")

# Fig3C
makeTCRcircle <- function(so){
  circles <- getCirclize(so, group.by = "annotationnew")
  grid.cols <- c(
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
    "TEMRA_2" = "#FF61C3")
  
  p <- chordDiagram(circles, self.link = 1, grid.col = grid.cols)
  return(p)
}

pdf("CITE_CD_TCR_circle.pdf", width = 5, height = 5)
makeTCRcircle(CDobj)
dev.off()

pdf("CITE_CK_TCR_circle.pdf", width = 5, height = 5)
makeTCRcircle(CKobj)
dev.off()

#Fig3D
CD_allTCR <- CDobj[[c("annotationnew", "clonotype_id")]]
CD_allTCR <- CD_allTCR %>% filter(!is.na(clonotype_id))
trm2_data <- CD_allTCR[CD_allTCR$annotationnew == "TRM_2", ]
write.csv(table(trm2_data), file = "CD_TRM_2_clone.csv")

#Fig3E
resultCD <- CD_allTCR %>%
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

sets <- c("TEM_1", "eTreg_1", "Tfh_1", "TRM_1", "naïve T", "TRM_2", 
          "TCM", "Tfh_2", "TEM_2", "naïve Treg", "eTreg_2",
          "TEMRA_1", "TEMRA_2")

reduced_data <- binary_matrix[rowSums(binary_matrix[, sets]) > 1 , ]

reduced_data <- data.frame(reduced_data)
reduced_data_filtered <- reduced_data[reduced_data$TRM_2 >= 1, ]

trm2_data <- CD_allTCR[CD_allTCR$annotationnew == "TRM_2", ]

duplicate_clonotypes <- trm2_data$clonotype_id[duplicated(trm2_data$clonotype_id) | 
                                                 duplicated(trm2_data$clonotype_id, fromLast = TRUE)]

TRM_2_dupli_clone <- unique(duplicate_clonotypes)
noncommon_from_TRM_2_dupli_clone <- setdiff(TRM_2_dupli_clone, rownames(reduced_data_filtered))
TRM_2_dupli_clone_num <- length(noncommon_from_TRM_2_dupli_clone)

new_rows <- data.frame(
  TEM_1 = rep(0, TRM_2_dupli_clone_num),
  eTreg_1 = rep(0, TRM_2_dupli_clone_num),
  Tfh_1 = rep(0, TRM_2_dupli_clone_num),
  TRM_1 = rep(0, TRM_2_dupli_clone_num),
  naïve.T = rep(0, TRM_2_dupli_clone_num),
  TRM_2 = rep(1, TRM_2_dupli_clone_num),
  TCM = rep(0, TRM_2_dupli_clone_num),
  Tfh_2 = rep(0, TRM_2_dupli_clone_num),
  TEM_2 = rep(0, TRM_2_dupli_clone_num),
  naïve.Treg = rep(0, TRM_2_dupli_clone_num),
  eTreg_2 = rep(0, TRM_2_dupli_clone_num),
  TEMRA_1 = rep(0, TRM_2_dupli_clone_num),
  TEMRA_2 = rep(0, TRM_2_dupli_clone_num),
  row.names = noncommon_from_TRM_2_dupli_clone
)

reduced_data_filtered2 <- rbind(reduced_data_filtered, new_rows)

sets <- c("TRM_2", "TEM_1", "eTreg_1", "Tfh_1", "TRM_1",   
          "TCM", "Tfh_2", "TEM_2", "eTreg_2", "TEMRA_2")
sets <- rev(sets)

pdf("CITE_CD_TCR_upset_TRM_2_new_2_2.pdf", width = 9, height = 15)
upset(reduced_data_filtered2, sets = sets, order.by = "freq", nintersects = NA, keep.order = TRUE)
dev.off()

#Fig3F
p <- DimPlot(obj, reduction = "umap.rpca", label = F, group.by = "plotclone") + NoLegend() 
plot_data <- p$data


plot_data_no_others <- plot_data %>% 
  filter(plotclone != "others")

create_combinations <- function(data) {
  data %>%
    group_by(plotclone) %>%
    do({
      points <- .
      combs <- expand.grid(id1 = points$id, id2 = points$id) %>%
        filter(id1 < id2) %>%
        left_join(points, by = c("id1" = "id")) %>%
        rename(x1 = umaprpca_1, y1 = umaprpca_2) %>%
        left_join(points, by = c("id2" = "id")) %>%
        rename(x2 = umaprpca_1, y2 = umaprpca_2, plotclone = plotclone.y) %>%
        select(x1, y1, x2, y2, plotclone)
      combs
    }) %>%
    ungroup()
}

plot_data_no_others <- plot_data_no_others %>%
  mutate(id = row_number())

line_data <- create_combinations(plot_data_no_others)

p + 
  geom_point(aes(x = umaprpca_1, y = umaprpca_2, 
                 size = ifelse(plotclone == "others", 0.1, 3), 
                 color = plotclone)) +
  geom_segment(data = line_data, 
               aes(x = x1, y = y1, xend = x2, yend = y2, color = plotclone),
               size = 0.5, alpha = 0.5) +  
  scale_size_identity() +
  scale_color_manual(values = current_colors) + NoAxes() + ggtitle(NULL)

pdf("CITE_CD_TCR_clone_line_UMAP.pdf", width = 6, height = 6)
p +
  geom_point(data = plot_data %>% filter(plotclone == "others"),
             aes(x = umaprpca_1, y = umaprpca_2, 
                 size = 0.1, color = plotclone),
             alpha = 0.5) + 
  geom_segment(data = line_data, 
               aes(x = x1, y = y1, xend = x2, yend = y2, color = plotclone),
               size = 0.5, alpha = 0.5) + 
  geom_point(data = plot_data_no_others,
             aes(x = umaprpca_1, y = umaprpca_2, 
                 size = 2, color = plotclone)) + 
  scale_size_identity() +
  scale_color_manual(values = current_colors) + NoAxes() + ggtitle(NULL)
dev.off()






