setwd("/orange/zhou/projects/aging/SDIII/Flycellatlas/Flycellatlas_wholebody_relaxed/")

load("WB_Relaxed_ENV.RData")

library(Seurat)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table)
library(tidyverse)
library(clustree)
library(stringr)
#remotes::install_github("mojaveazure/seurat-disk")

loom.fn <- './r_fca_biohub_body_10x.loom'
s_cnct <- Connect(filename = loom.fn, mode = "r")



# gene list
n_genes <- s_cnct[['row_attrs']][['Gene']][['dims']]
gns <- s_cnct[['row_attrs']][['Gene']][1:n_genes]

# cell ID list
n_cells <- s_cnct[['col_attrs']][['CellID']][['dims']]
cellids <- s_cnct[['col_attrs']][['CellID']][1:n_cells]

# get raw counts matrix
raw.cnts <- as.data.frame(t(s_cnct[['matrix']][,] ))
names(raw.cnts) <- cellids
rownames(raw.cnts) <- gns



metadata <- data.frame(
  cellID = cellids,
  ClusterID = s_cnct[['col_attrs']][['ClusterID']][1:n_cells],
  annotation = s_cnct[['col_attrs']][['S_annotation']][1:n_cells]
)
rownames(metadata) <- cellids





s_obj <- CreateSeuratObject(counts = raw.cnts,
                            project = "fromLoom",
                            assay = "RNA",
                            meta.data = metadata)



save(s_obj,file = "./WB_relaxed_SeuratOBJ.RData")




rm(list=ls())


assign('WB_relaxed', get(load('./WB_relaxed_SeuratOBJ.RData')))



View(WB_relaxed@'meta.data')




WB_relaxed[["percent.mt"]] <- PercentageFeatureSet(WB_relaxed, pattern = "^mt:")
WB_relaxed[["percent.rb"]] <- PercentageFeatureSet(WB_relaxed, pattern = "^Rb")


# 
# 
# 
# 
# n.Feature.min <- 200  #least gene # >200
# n.Feature.max <- 4000  #max gene # <4000
# n.Count.min <-500 #least count >500
# n.Mt <- 10  #max mt <20%
# n.Rb <- 10   #max ribosome <10% *** changeable
# 
# cat("Before filter :",nrow(WB_relaxed@meta.data),"cells\n")
# 
# 
# WB_relaxed <- subset(WB_relaxed, 
#                    subset = 
#                      nFeature_RNA > n.Feature.min &
#                      nFeature_RNA<n.Feature.max &
#                      nCount_RNA >n.Count.min&
#                      percent.mt < n.Mt&
#                      percent.rb < n.Rb)

cat("After filter :",nrow(WB_relaxed@meta.data),"cells\n")






# run standard anlaysis workflow
WB_relaxed2 <- NormalizeData(WB_relaxed)
WB_relaxed2 <- FindVariableFeatures(WB_relaxed2)
WB_relaxed2 <- ScaleData(WB_relaxed2)

min(nrow(WB_relaxed2), ncol(WB_relaxed2))
WB_relaxed2 <- RunPCA(WB_relaxed2,npcs = 78,verbose = F)

#see turWB_relaxedng point
pdf("./WB_relaxed_elbowplot.pdf",width = 60,height = 40)
ElbowPlot(WB_relaxed2, ndims=100, reduction="pca")
dev.off()


WB_relaxed2 <- FindNeighbors(WB_relaxed2, dims = 1:78, reduction = "pca", verbose=TRUE)

WB_relaxed2 <- FindClusters(WB_relaxed2, resolution = 0.8, cluster.name = "WB_relaxed_clusters")

View(WB_relaxed2@meta.data)
library(clustree)
#see resolution
# pdf(paste0("./cluster_tree.pdf"),width = 60,height = 40)
# clustree(FindClusters(WB_relaxed2, resolution = seq(0.3,0.4,by=0.05)))
# dev.off()



WB_relaxed2 <- FindClusters(WB_relaxed2, resolution = 0.6, cluster.name = "WB_relaxed_clusters")


WB_relaxed3 <- RunTSNE(WB_relaxed3, dims = 1:78, reduction = "pca", reduction.name = "tSNE_WB_relaxed",verbose=T)
#WB_relaxed3 <- RunUMAP(WB_relaxed3, dims = 1:78, reduction = "pca", reduction.name = "UMAP_WB_relaxed",verbose=T)

#View(WB_relaxed@meta.data)




pdf(paste0("./DimPlot_WB_relaxed_tSNE.pdf"),width = 30,height = 20)
DimPlot(WB_relaxed3, reduction = "tSNE_WB_relaxed",group.by = "annotation",
        label = TRUE,
        label.size = 3,
        repel = TRUE,
        raster=FALSE,)
dev.off()




pdf(paste0("./DimPlot_WB_relaxed_UMAP.pdf"),width = 30,height = 20)
DimPlot(WB_relaxed3, reduction = "UMAP_WB_relaxed",group.by = "annotation",
        label = TRUE,
        label.size = 3,
        repel = TRUE,
        raster=FALSE,)
dev.off()



#UMAP feature plot
pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_SPH93.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("SPH93"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_CecA1.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("CecA1"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_AttA.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("AttA"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_AttB.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("AttB"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_CecA2.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("CecA2"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_AttC.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("AttC"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_DptA.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("DptA"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_DptB.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("DptB"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_Dro.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("Dro"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_Mtk.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("Mtk"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_IM18.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("IM18"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/UMAP_WB_relaxed_featureplot_edin.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("edin"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "UMAP_WB_relaxed",keep.scale = "feature")
dev.off()




#tSNE feature plot
pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_SPH93.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("SPH93"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_CecA1.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("CecA1"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_AttA.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("AttA"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_AttB.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("AttB"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_CecA2.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("CecA2"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_AttC.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("AttC"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_DptA.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("DptA"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_DptB.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("DptB"),  cols = c("lightblue2","darkred"),raster=FALSE, reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_Dro.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("Dro"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_Mtk.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("Mtk"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_IM18.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("IM18"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()

pdf(paste0("./FeaturePlot/tSNE_WB_relaxed_featureplot_edin.pdf"),width = 15,height = 10)
FeaturePlot(WB_relaxed3, features = c("edin"),  cols = c("lightblue2","darkred"), raster=FALSE,reduction = "tSNE_WB_relaxed",keep.scale = "feature")
dev.off()


WB_relaxed3@meta.data$"Sex" <- str_split(WB_relaxed3@meta.data$cellID,"_",simplify =T)[,4]

View(WB_relaxed3@meta.data)

unique(WB_relaxed3@meta.data$Sex)




Plot_compareWithBulkRNAseq <- function(SeuratOBJ, SEX, GENELIST){
  EXPR <- as.data.frame(
    AverageExpression(object = SeuratOBJ,
                      features = GENELIST,
                      group.by = c("Sex","annotation")
                      )$RNA
  );
  tmp1 <- EXPR[,str_split(colnames(EXPR),pattern = "_",simplify =T)[,1]== SEX];
  gplots::heatmap.2(as.matrix(tmp1),col = blues9,margins = c(20,20))
  
}

# 
# 
# test <- as.data.frame(
#   AverageExpression(object = WB_relaxed3,
#                     features = c("CecA1","AttA"),
#                     group.by = c("Sex","annotation")
#                     )$RNA
#   )
# 
# str_split(colnames(test),pattern = "_",simplify =T)[,1]=="Male"
# 
# 



pdf("./Heatmap_IIIMarkers_Male_celltype.pdf",width = 15,height = 10)
Plot_compareWithBulkRNAseq(SeuratOBJ=WB_relaxed3, SEX="Male", GENELIST=c("SPH93",
                                                                          "CecA1",
                                                                          "CecA2",
                                                                          "AttA",
                                                                          "AttB",
                                                                          "AttC",
                                                                          "DptA",
                                                                          "DptB",
                                                                          "Dro",
                                                                          "Mtk",
                                                                          "IM18",
                                                                          "edin"))
dev.off()




pdf("./Heatmap_IIIMarkers_Female_celltype.pdf",width = 15,height = 10)
Plot_compareWithBulkRNAseq(SeuratOBJ=WB_relaxed3, SEX="Female", GENELIST=c("SPH93",
                                                                         "CecA1",
                                                                         "CecA2",
                                                                         "AttA",
                                                                         "AttB",
                                                                         "AttC",
                                                                         "DptA",
                                                                         "DptB",
                                                                         "Dro",
                                                                         "Mtk",
                                                                         "IM18",
                                                                         "edin"))
dev.off()




#









save.image("WB_Relaxed_ENV.RData")




