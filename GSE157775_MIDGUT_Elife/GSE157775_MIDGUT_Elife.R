setwd("/orange/zhou/projects/II_Cancer/GSE157775_MIDGUT_Elife/")


library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(data.table)
library(monocle)
#library(DoubletFinder)
library(tidyverse)
library(patchwork)
library(clustree)
library(remotes)
library(data.table)
library(ggrepel)

load("./MIDGUT_Elife_wkimage.rdata")

# 
# 
# MIDGUT_Elife <- Read10X("./GSE157775_data/")
# 
# 
# MIDGUT_Elife <- CreateSeuratObject(counts = MIDGUT_Elife,
#                            project = "MIDGUT_Elife",
#                            min.cells = 0,
#                            min.features = 0)
# 
Author_metadata <- fread("./GSE157775_data/GSE157775_cell_annotation.csv")
# 
# View(MIDGUT_Elife@meta.data)
# 
# nrow(MIDGUT_Elife@meta.data)

########meta data does not align with seurat obj


##try loading author's work image of Seurat object

MIDGUT_Elife <- readRDS("./GSE157775_data/GSE157775_sce.RDS")

MIDGUT_Elife <- as.Seurat(MIDGUT_Elife)

MIDGUT_Elife@meta.data$annotation <- Author_metadata$assigned_cell_type

MIDGUT_Elife@meta.data$Age <- Author_metadata$age

MIDGUT_Elife@meta.data$UMAP_1 <- Author_metadata$UMAP_1
MIDGUT_Elife@meta.data$UMAP_2 <- Author_metadata$UMAP_2

#no Sex information in this study


View(MIDGUT_Elife@meta.data)






Idents(MIDGUT_Elife) <- MIDGUT_Elife@meta.data$annotation


##run the standard pipeline of Seurat to build reduction component in Seurat object



# run standard anlaysis workflow
MIDGUT_Elife2 <- NormalizeData(MIDGUT_Elife)
MIDGUT_Elife2 <- FindVariableFeatures(MIDGUT_Elife2)
MIDGUT_Elife2 <- ScaleData(MIDGUT_Elife2)

min(nrow(MIDGUT_Elife2), ncol(MIDGUT_Elife2))
MIDGUT_Elife2 <- RunPCA(MIDGUT_Elife2,npcs = 50,verbose = F)
# 
# #see turMIDGUT_Elifeng point
# pdf("./MIDGUT_Elife_elbowplot.pdf",width = 30,height = 20)
# ElbowPlot(MIDGUT_Elife2, ndims=50, reduction="pca")
# dev.off()


MIDGUT_Elife2 <- FindNeighbors(MIDGUT_Elife2, dims = 1:30, reduction = "pca", verbose=TRUE)

#MIDGUT_Elife2 <- FindClusters(MIDGUT_Elife2, resolution = 0.2, cluster.name = "MIDGUT_Elife_clusters")
# 
# #View(MIDGUT_Elife2@meta.data)
# library(clustree)
# #see resolution
# pdf(paste0("./cluster_tree.pdf"),width = 60,height = 40)
# clustree(FindClusters(MIDGUT_Elife2, resolution = seq(0.2,0.6,by=0.05)))
# dev.off()



MIDGUT_Elife2 <- FindClusters(MIDGUT_Elife2, resolution = 0.6, cluster.name = "MIDGUT_Elife_clusters")


#MIDGUT_Elife3 <- RunTSNE(MIDGUT_Elife2, dims = 1:30, reduction = "pca", reduction.name = "tSNE_MIDGUT_Elife",verbose=T)
MIDGUT_Elife3 <- RunUMAP(MIDGUT_Elife2, dims = 1:30, reduction = "pca", reduction.name = "UMAP_MIDGUT_Elife",verbose=T)



#change UMAP reduction coordinates to Author's

head(MIDGUT_Elife3@reductions[["UMAP_MIDGUT_Elife"]]@cell.embeddings)

MIDGUT_Elife3@reductions[["UMAP_MIDGUT_Elife"]]@cell.embeddings[,1] <- MIDGUT_Elife3@meta.data$UMAP_1

MIDGUT_Elife3@reductions[["UMAP_MIDGUT_Elife"]]@cell.embeddings[,2] <- MIDGUT_Elife3@meta.data$UMAP_2

head(MIDGUT_Elife3@reductions[["UMAP_MIDGUT_Elife"]]@cell.embeddings)


Idents(MIDGUT_Elife3) <- MIDGUT_Elife3@meta.data$annotation


#Now plot the UMAP Dimplot

pdf("./MIDGUT_Elife_UMAP_Dimplot.pdf",width = 9,height = 6)
DimPlot(MIDGUT_Elife3,reduction = "UMAP_MIDGUT_Elife",raster = F,label = T)
dev.off()


###################################################################

#convert FBgn to gene symbol in counts slot
rownames(MIDGUT_Elife3@assays$originalexp$counts)


dim(MIDGUT_Elife3@assays$originalexp$counts)
#[1] 12872 12149

MIDGUT_Elife3@assays$originalexp$counts



originalexp_counts <- as.data.frame(MIDGUT_Elife3@assays$originalexp$counts)
#########mapping FBGN to gene symbol

annotation_table <- fread("./FBgn2Symbol.csv",header = F)


ID_to_symbol <- function(expr_data){
  df=expr_data
  if(sum(colnames(df)=="gene_name", na.rm=TRUE)){
    df$gene_name <- ""
    for (i in 1:nrow(df)){
      if(df$ID[i]%in%annotation_table$V1==TRUE){
        df$gene_name[i] <-
          annotation_table$V2[which(df$ID[i]==annotation_table$V1)]
      }else{
        df$gene_name[i] <- "No_map"
      }
    }
    print(
      paste(sum(df$gene_name=="No_map", na.rm=TRUE)," _unmapped",sum(df$gene_name!="No_map", na.rm=TRUE)," mapped")
    );
  }
  else{
    print("there is no column called 'gene_name', use df$gene_name <- ''");
  }
  return(df)
}

originalexp_counts$ID <- rownames(originalexp_counts)
originalexp_counts$gene_name <- ''

originalexp_counts2 <- ID_to_symbol(originalexp_counts)
#"32  _unmapped 12840  mapped"


which(originalexp_counts2$gene_name=="No_map")

##get the FBGN of these unmapped ID

originalexp_counts2$ID[which(originalexp_counts2$gene_name=="No_map")]
####might need to manually map these with flybase ID validator

write.table(file = "./Unampped_ID.txt",
            originalexp_counts2$ID[which(originalexp_counts2$gene_name=="No_map")],
            sep = ",",
            quote = F,
            col.names = F,
            row.names = F
            )
####
#import unmapped2symbol done based on flybase ID validator

unmapped2symbol <- fread("./FlyBase_Fields_download_Unmapped2Symbol.txt",header = T)

originalexp_counts2$gene_name[which(originalexp_counts2$gene_name=="No_map")] <- unmapped2symbol$SYMBOL

originalexp_counts2$gene_name[which(originalexp_counts2$gene_name=="No_map")]
#successfully mapped


#convert FBgn in counts slot in midgut seurat obj into gene symbol, so feature plot can be made
MIDGUT_Elife3@assays[["originalexp"]]@counts@Dimnames[[1]] <- originalexp_counts2$gene_name

MIDGUT_Elife3@assays[["originalexp"]]@counts@Dimnames[[1]]



#let's plot a featureplot of hml to see if there is potentially a subset cluster of plasmatocytes.


FeaturePlot(MIDGUT_Elife3, features = c("SPH93"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")



FeaturePlot(MIDGUT_Elife3, features = c("CecA1"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")


FeaturePlot(MIDGUT_Elife3, features = c("AttA"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")


FeaturePlot(MIDGUT_Elife3, features = c("edin"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")


FeaturePlot(MIDGUT_Elife3, features = c("Dro"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")


FeaturePlot(MIDGUT_Elife3, features = c("IM18"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")


FeaturePlot(MIDGUT_Elife3, features = c("Hml"),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")


FeaturePlot(MIDGUT_Elife3, features = c(""),  cols = c("azure3","darkred"),raster=FALSE, reduction = "UMAP_MIDGUT_Elife",keep.scale = "feature",slot = "counts")












#########################################################################################

save.image("MIDGUT_Elife_wkimage.rdata")
