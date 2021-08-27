#load library ####
#install.packages('Seurat')
library(Seurat)
library(RNAransform)
library(dplyr)
library(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(GEOquery)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DropletUtils")
library(DropletUtils)
#instRound0.packages("devtools")
library(devtools)
#source("https://raw.githubusercontent.com/farrellja/URD/master/URD-InstRound0.R")
library(URD)
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(reticulate)
library(scales)
library(forcats)
library(cowplot)
#install.packages("magrittr")
library(magrittr)
#install.packages('varhandle')
library(varhandle)
#install.packages("googlesheets4")
library(googlesheets4)

#install.packages('BiocManager')
#BiocManager::install('limma')
library(limma)

#BiocManager::install('multtest')
#install.packages('metap')

library(multtest)
library(metap)
library(reshape2)

######################################################################################################
#Start ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
getwd()
list.files()
rm(list=ls())

load("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated.Rda")
load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated.Rda")

rm(list= ls()[!(ls() %in% c('Round0'))])
#####################################################################################################
## (SKIP) Perform an integrated analysis ####
DefaultAssay(Round0) <- "integrated"
colnames(Round0@meta.data)

# Run the standard workflow for visualization and clustering
Round0 <- ScaleData(Round0, verbose = FALSE)
Round0 <- RunPCA(Round0, npcs = 60, verbose = FALSE)

# t-SNE and Clustering
Round0 <- FindNeighbors(Round0, reduction = "pca", dims = 1:60)

Round0 <- FindClusters(Round0, resolution = 0.8)
#Round0 <- FindClusters(Round0, resolution = 0.7)

Round0 <- RunUMAP(object = Round0, reductiosdfn = "pca", dims = 1:60)
######################################################################################################
## (SKIP) Dimensionality reduction ####
## ElbowPlot
pdf('./combined_data/ElbowPlot_Round0.pdf',width=10,height=10)
ElbowPlot(object = Round0, ndims = 60)
dev.off()

## JackStraw plot
#Round0 <- JackStraw(object = Round0, num.replicate = 100)
#Round0 <- ScoreJackStraw(object = Round0, dims = 1:20)

#pdf('JackStrawPlot_Round0.pdf',width=10,height=10)
#JackStrawPlot(object = Round0, dims = 1:20)
#dev.off()

## (SKIP) Perform linear dimensional reduction 
pdf('./combined_data/VizDimLoadings_Round0.pdf', width=8,height=7)
VizDimLoadings(Round0, dims = 1:2, reduction = "pca")
dev.off()

pdf('./combined_data/DimPlots_Round0.pdf', width=8,height=7)
DimPlot(Round0, reduction = "pca")
dev.off()

#####################################################################################################
## chemistry cluster and "RNA" scaledata ####
colnames(Round0@meta.data)
Idents(object = Round0) <- ("number")
levels(Idents(Round0)) 

category <- as.data.frame(levels(Idents(Round0)) )
colnames(category) <- "Level"

phenotype <- read_sheet("https://https://docs.google.com/spreadsheets/d/1_cVcmP7AniLrhIE9XWw-T-_BqQXIkG8wi7wl1Rb3kl4/edit#gid=0&fvid=109363522")
phenotype <- phenotype[grep("Yes", phenotype$Third_project), ]
#Remove data with something wrong "Psoriasis09_postTx_week12_LS"
#phenotype <- phenotype[!phenotype$Rcoding_Number == "Psoriasis09_postTx_week12_LS" ,]

colnames(phenotype)
version <- phenotype[,c("Rcoding_Category","Rcoding_Number","Chemistry_2")]
#version$version <- paste(version$Rcoding_Category,version$Chemistry, sep="_")

version_sorted <- merge(category, version, by.x = "Level", by.y ="Rcoding_Number" , all.x = TRUE,sort=FALSE)

version_sorted$Level

version_sorted <- as.character(version_sorted$Chemistry_2)

new.cluster.ids <- version_sorted
                      
names(new.cluster.ids) <- levels(Round0)
Round0 <- RenameIdents(Round0, new.cluster.ids)
levels(Idents(Round0))
Idents(Round0) <- fct_relevel(Idents(Round0), sort)
Round0[["chemistry"]] <- Idents(object = Round0)

colnames(Round0@meta.data)
DefaultAssay(Round0) <- "RNA"
#Round0 <- ScaleData(Round0, verbose = FALSE)
#####################################################################################################
## Cluster level order change #1 ####
colnames(Round0@meta.data)
Idents(object = Round0) <- ("integrated_snn_res.0.8")
levels(Round0)

Idents(object = Round0) <- factor(Idents(object = Round0), levels = as.character(c(0:(nlevels(Round0)-1))))    

Round0[["integrated_snn_res.0.8"]] <- Idents(object = Round0)
#####################################################################################################
## Cluster level order change #2 ####
colnames(Round0@meta.data)
Idents(object = Round0) <- ("stim")
levels(Idents(Round0)) 
Idents(object = Round0) <- factor(Idents(object = Round0), levels = c("Control","Psoriasis_NL","Psoriasis_Pre","Psoriasis_Tx"))      
Round0[["stim"]] <- Idents(object = Round0)

#####################################################################################################
## Cluster level order change #3  ####
colnames(Round0@meta.data)
Idents(object = Round0) <- ("number")
levels(Round0)

Idents(object = Round0) <- factor(Idents(object = Round0), levels = c( "Control01"   ,
                                                                       "Control02"    ,
                                                                       "Control03"         ,
                                                                       "Control04"     ,
                                                                       "Control05"   ,
                                                                       "Control06"    ,
                                                                       "Control07"      , 
                                                                        "Control08"     ,
                                                                       "Control09"   ,
                                                                       "Control10"   ,
                                                                       "Psoriasis01_preTx_LS"   ,
                                                                       "Psoriasis02_preTx_LS"   ,
                                                                       "Psoriasis03_preTx_LS"  ,
                                                                       "Psoriasis04_preTx_LS"    ,    
                                                                       "Psoriasis04_postTx_week12_LS",
                                                                       "Psoriasis05_preTx_LS"   ,
                                                                       "Psoriasis05_postTx_week12_LS",
                                                                       "Psoriasis06_preTx_LS"    ,
                                                                       "Psoriasis07_preTx_LS"    ,
                                                                       "Psoriasis07_preTx_LS_02"   , 
                                                                       "Psoriasis07_preTx_NL_02"  ,
                                                                       "Psoriasis07_postTx_week12_LS" ,
                                                                       "Psoriasis07_postTx_week48_LS",
                                                                       "Psoriasis08_preTx_LS"  ,      
                                                                       "Psoriasis08_postTx_week12_LS",
                                                                       "Psoriasis08_postTx_week24_LS" ,
                                                                       "Psoriasis09_preTx_LS"    ,
                                                                       "Psoriasis10_preTx_LS"    
))

Round0[["number"]] <- Idents(object = Round0)

#####################################################################################################
## Visualization by clusters without labels####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

dir.create("./UMAP_clusters")
setwd("./UMAP_clusters")
getwd()
##cluster count 
Idents(object = Round0) <- ("integrated_snn_res.0.8")
colnames(Round0@meta.data)
cell.num <- table(Round0@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round0_cluster_count_integrated_snn_res.0.8.pdf', width=12,height=8)
DimPlot(object = Round0,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

pdf('UMAP_Round0_cluster_count_integrated_snn_res.0.8_noLabel.pdf', width=12,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##DimPlot
p1 <- DimPlot(Round0, reduction = "umap", group.by = "stim")+
  scale_color_manual(values=c('green','coral1','deeppink','deep sky blue'))

p2 <- DimPlot(Round0, reduction = "umap", label = FALSE)

pdf('UMAP_DimPlot_integrated_snn_res.0.8.pdf', width=14,height=7)
plot_grid(p1, p2)
dev.off()

##cluster count 
colnames(Round0@meta.data)
cell.num <- table(Round0@meta.data$integrated_snn_res.0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cluster_count_integrated_snn_res.0.8.pdf', width=10,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###
pdf('UMAP_DimPlot_split_integrated_snn_res.0.8.pdf', width=20,height=7)
DimPlot(Round0, reduction = "umap", split.by = "stim")
dev.off()


##Cell number
colnames(Round0@meta.data)
Idents(object = Round0) <- ("number")

cell.num <- table(Round0@meta.data$number)
ClusterLabels = paste(names(cell.num), paste0("(n=", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cellcount_integrated_snn_res.0.8.pdf', width=14,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

##
colnames(Round0@meta.data)
Idents(object = Round0) 

Idents(object = Round0) <- ("stim")
cell.numb <- table(Round0@meta.data$stim)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round0_category_integrated_snn_res.0.8.pdf',width=10,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()







#####################################################################################################
## Visualization by chemistry ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

dir.create("./UMAP_chemistry")
setwd("./UMAP_chemistry")

colnames(Round0@meta.data)
DefaultAssay(Round0) <- "RNA"
Idents(object = Round0) <- ("stim")
#Idents(object = Round0) <- factor(Idents(object = Round0), levels = c("Control","Psoriasis"))         
Idents(object = Round0) <- ("chemistry")
levels(Round0)

##DimPlot
p1 <- DimPlot(Round0, reduction = "umap", group.by = "stim")+
  scale_color_manual(values=c('green','coral1','deeppink','deep sky blue'))

p2 <- DimPlot(Round0, reduction = "umap", label = FALSE)

pdf('UMAP_Round0_DimPlot.pdf', width=15,height=7)
plot_grid(p1, p2)
dev.off()

##cluster count 
colnames(Round0@meta.data)
#cell.num <- table(Round0@meta.data$integrated_snn_res.0.8)
cell.num <- table(Round0@meta.data$chemistry)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round0_cluster_count.pdf', width=10,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###
pdf('UMAP_Round0_DimPlot_split.pdf', width=14,height=7)
DimPlot(Round0, reduction = "umap", split.by = "chemistry")
dev.off()


##
colnames(Round0@meta.data)
Idents(object = Round0) 

Idents(object = Round0) <- ("stim")
cell.numb <- table(Round0@meta.data$stim)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round0_stim.pdf',width=10,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


#####################################################################################################
## Dot plot -1 to check clusters ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

Round0@assays
colnames(Round0@meta.data)
#Round0@meta.data$integrated_snn_res.0.8
Round0@meta.data$integrated_snn_res.0.8
DefaultAssay(Round0)
levels(Idents(object = Round0))

#Idents(object = Round0) <- ("integrated_snn_res.0.8")
Idents(object = Round0) <- ("integrated_snn_res.0.8")
#Idents(object = Round0) <- ("RNA_snn_res.0.8")

DefaultAssay(Round0) <- "RNA"

features.plot <- c(  "SCGB2A2","DCD","COL1A1","DCN","CD34","LYVE1","CCL21",
                     "KRT16","KRT17","KRT6B","KRT5","KRT14", "KRT1","KRT10",
                     'FABP5',
                     "CDSN",
                     "LCE3D", 
                     "SPRR2G",
                     "LCE2C","LCE1A",
                     "IL1RL1","IL18",
                     "APOBEC3A","S100A9",
                     "SPINK5",
                     "CCL20","PRSS22",
                     "SPRR2D","SPRR2A",
                     "KRT23","KRT80",
                     "MLANA", "TYRP1", "DCT",
                     "THBD","CD1C","CLEC10A",
                     "CD163",
                     "CD14","LYZ",
                     "HLA-DRB5","HLA-DRA", "HLA-DRB1", 
                     "HLA-DQB1","HLA-DQA1",
                     "CCL22",
                     "CD40","CIITA","LY75","LAMP3",
                     "CTLA4","FOXP3","IL2RA","TIGIT", 
                     "IGLC2",
                     "PTPRC",
                     "CD8B","CD8A","GZMK","GZMH",
                     "TRBC1","TRAC",
                     "CD3D",
                     "GNLY", "KLRB1")

features.plot <- rev(features.plot)

#Idents(object = Round0) <- factor(Idents(object = Round0), 
#                                 levels = rev(c(0,4,10,14,15,19,17,6,9,16,1,2,3,5,7,8,11,12,13,18)))


Idents(object = Round0) <- fct_rev(Idents(object = Round0))

pdf('Dot plot_res0.8_initial.pdf', width=16,height=10)
DotPlot(object = Round0, features = features.plot) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=18))
dev.off()

pdf('Dot plot_res0.8_split.pdf', width=18,height=28)
DotPlot(Round0, features = features.plot, cols = c("blue", "red", "blue", "black"), dot.scale = 8, 
        split.by = "stim") + RotatedAxis()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=18))
dev.off()


#cluster3vs5 <- FindMarkers(Round0, ident.1 = 3, ident.2 = 5, min.pct = 0.25)
#write.csv(cluster3vs5, file = "cluster3vs5.csv")
Idents(object = Round0) <- fct_rev(Idents(object = Round0))

setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")



#####################################################################################################
## DEG without labeling ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

Round0@assays
colnames(Round0@meta.data)
DefaultAssay(Round0)
Idents(object = Round0) <- ("integrated_snn_res.0.8")
levels(Idents(object = Round0))

DefaultAssay(Round0) <- "RNA"
Round0 <- ScaleData(Round0, verbose = FALSE)

Combined.markers <- FindAllMarkers(object = Round0, only.pos = TRUE, min.pct = 0.25, 
                                   logfc.threshold = 0.25)

write.csv(Combined.markers, file = "Combined.markers_all.csv")

top10 <- Combined.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
pdf('Diff_Heatmp_Round0.pdf', width=24, height=30)
DoHeatmap(object = Round0, features = top10$gene) + NoLegend()
dev.off()



#####################################################################################################
## Feature plot by samples ####
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
dir.create("./Plots_samples")
setwd("./Plots_samples")
getwd()

DefaultAssay(Round0) <- "RNA"
colnames(Round0@meta.data)
Idents(object = Round0) <- ("number")
levels(Idents(object = Round0))


marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])


for (i in 1:length(marker)){
  pdf(paste(marker[i],"_VlnPlot_Round0.pdf",sep=''), width=14,height=4)
  
  p <- VlnPlot(Round0, features=as.character(marker[i]),  split.by = "stim", cols = c("cyan", "red", "green"))
  print(p)
  dev.off()
  
  pdf(paste(marker[i],"_FeaturePlot_Round0.pdf",sep=''), width=18,height=6)
  #p <- FeaturePlot(object = Round0, split.by = "stim", features = as.character(marker[i]),min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
  p <- FeaturePlot(object = Round0, split.by = "stim", features = as.character(marker[i]),label=FALSE, order=TRUE) 
  print(p)
  dev.off()
}


setwd("../")

########################################################################
# (From here) Dot plot with labels ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

Round0@assays
colnames(Round0@meta.data)
Round0@meta.data$integrated_snn_res.0.8
DefaultAssay(Round0)
levels(Idents(object = Round0))

Idents(object = Round0) <- ("integrated_snn_res.0.8")

DefaultAssay(Round0) <- "RNA"

new.cluster.ids <- c(
  "KC-S.Corneum.01","Mature_DC.01","KC-S.Corneum.02","CD4_T_cell.01","KC-S.Corneum.03","CD8_T_cell","KC-S.Basale","KC-S.Corneum.04","KC-S.Corneum.05","KC-S.Corneum.06","KC-S.Corneum.07","Treg","Melanocyte","KC-S.Spinosum","Endothelial_cell","Semimature_DC.01","Mature_DC.02","NK_cell","KC-S.Granulosum","CD4_T_cell.02","Doublet.01","CD161_T_cell","CD4_T_cell.03","KC-Hair_follicle","Fibroblast","Doublet.02","Eccrine_gland","Semimature_DC.02"
)	

names(new.cluster.ids) <- levels(Round0)
Round0 <- RenameIdents(Round0, new.cluster.ids)

levels(Idents(Round0))

Idents(object = Round0) <- factor(Idents(object = Round0), 
                                  levels = (c(
                                    "NK_cell","CD161_T_cell","CD8_T_cell","CD4_T_cell.01","CD4_T_cell.02","CD4_T_cell.03","Treg","Mature_DC.01","Mature_DC.02","Semimature_DC.01","Semimature_DC.02","Melanocyte","KC-S.Corneum.01","KC-S.Corneum.02","KC-S.Corneum.03","KC-S.Corneum.04","KC-S.Corneum.05","KC-S.Corneum.06","KC-S.Corneum.07","KC-S.Granulosum","KC-S.Spinosum","KC-S.Basale","KC-Hair_follicle","Endothelial_cell","Fibroblast","Eccrine_gland","Doublet.01","Doublet.02"
                                  )))

Round0[["ClusterNames_0.8"]] <- Idents(object = Round0)

####
Idents(object = Round0) <- ("ClusterNames_0.8")
features.plot <- c(  "KLRB1"   , "GNLY" ,  
                     "CD3D" ,    "TRAC" ,    "TRBC1" ,
                     "CD8A"   ,  "CD8B"   ,
                     "IFNG","IL17A","IL17F",
                     "GZMK"    ,"GZMH"   ,
                     # "PTPRC"  ,  
                     "IGLC2" ,  
                     "TIGIT" ,   "IL2RA"   ,  "FOXP3"   , "CTLA4" ,  
                     "LAMP3"  ,  "CCL22"  , "LY75"  ,   "CIITA"  ,  "CD40"  ,  
                     "HLA-DQA1","HLA-DQB1", "HLA-DRB1", "HLA-DRA" , "HLA-DRB5" ,
                     "LYZ"     ,  "CD14"   ,  "CD163"   ,  "THBD"    ,"IL10",
                     "CLEC10A" ,                      "CD1C"   , 
                     "DCT"    ,  "TYRP1"  ,  "MLANA"   ,
                     "CDSN" ,
                     "KRT80"   , "KRT23"   ,"SPINK5" , "PRSS22"  ,
                     "SPRR2A"  , "SPRR2D"  , 
                        "APOBEC3A", "IL1RL1" ,"IL18"   ,  
                      "LCE2C" ,     "LCE1A"  , 
                     "FABP5"  ,
                     "KRT10"   , "KRT1"   ,"CCL20"  ,  
                     "KRT14"  ,  "KRT5"  ,
                     "KRT6B"   , "KRT17"   , "KRT16"   ,
                     "CCL21"  ,  "LYVE1"   , "CD34"  ,   "DCN"    ,  "COL1A1" ,  "DCD"   ,   "SCGB2A2" )


Idents(object = Round0) <- fct_rev(Idents(object = Round0))


pdf('Dot plot_labeled01.pdf', width=20,height=10)
DotPlot(object = Round0, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=18))
dev.off()


pdf('Dot plot_labeled02.pdf', width=20,height=10)
DotPlot(object = Round0, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=18))
dev.off()


pdf('Dot plot_split_labeled.pdf', width=20,height=28)
DotPlot(Round0, features = features.plot, cols = c("black", "pink","red","blue"),  
        split.by = "stim") + RotatedAxis()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49))
dev.off()

Idents(object = Round0) <- fct_rev(Idents(object = Round0))
#


save(Round0, file = "./combined_data/Round0_integrated_analyzed_01.31.2021.Rda")


########################################################################
# Visualization with label ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

dir.create("./UMAP_Label")
setwd("./UMAP_Label")

Round0@assays
colnames(Round0@meta.data)
DefaultAssay(Round0) <- "RNA"

##cluster count 
Idents(object = Round0) <- ("ClusterNames_0.8")
colnames(Round0@meta.data)
cell.num <- table(Round0@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round0_cluster_count_ClusterNames_0.8.pdf', width=14,height=8)
DimPlot(object = Round0,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

pdf('UMAP_Round0_cluster_count_ClusterNames_0.8_noLabel.pdf', width=14,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##DimPlot
p1 <- DimPlot(Round0, reduction = "umap", group.by = "stim")+
  scale_color_manual(values=c('green','coral1','deeppink','deep sky blue'))

p2 <- DimPlot(Round0, reduction = "umap", label = FALSE)

pdf('UMAP_DimPlot_integrated_snn_res.0.8.pdf', width=18,height=7)
plot_grid(p1, p2)
dev.off()

##cluster count 
colnames(Round0@meta.data)
cell.num <- table(Round0@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cluster_count_integrated_snn_res.0.8.pdf', width=12,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 4, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###
pdf('UMAP_DimPlot_split_integrated_snn_res.0.8.pdf', width=20,height=7)
DimPlot(Round0, reduction = "umap", split.by = "stim")
dev.off()


##Cell number
colnames(Round0@meta.data)
Idents(object = Round0) <- ("number")

cell.num <- table(Round0@meta.data$number)
ClusterLabels = paste(names(cell.num), paste0("(n=", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_cellcount_integrated_snn_res.0.8.pdf', width=14,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()

##
colnames(Round0@meta.data)
Idents(object = Round0) 

Idents(object = Round0) <- ("stim")
cell.numb <- table(Round0@meta.data$stim)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round0_category_integrated_snn_res.0.8.pdf',width=10,height=8)
DimPlot(object = Round0,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


setwd("../")


#########################################
#gene expressing cell list #1 ####
#By samples
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

colnames(Round0@meta.data)
Idents(object = Round0) <- ("number")
levels(Idents(object = Round0))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])

#memory.limit()
#memory.limit(size=500000) 

expr <- FetchData(object = Round0, vars = marker[1])
list <- Round0[, which(x = expr > 0)]
list <- as.data.frame(table(list@active.ident))
colnames(list) <- c("number",as.character(marker[1]))
list1 <- list




for (i in 2:length(marker)){
  expr <- FetchData(object = Round0, vars = marker[i])
  list <- Round0[, which(x = expr > 0)]
  list <- as.data.frame(table(list@active.ident))
  colnames(list) <- c("number",as.character(marker[i]))
  list1 <- merge(list1,list, by="number", all = TRUE)
}

ncol(list1)
length(marker)

#list1 <- data.frame(t(list1))
#colnames(list1) <- as.character(unlist(list1[1,]))
#list1 <- list1[-1,]
list1[is.na(list1)] <- 0

write.csv(list1, file="gene_expressing_cells_cluster_expr_over_0.csv")

#rm(list)
#rm(list1)

#################################################################################################
##Feature plot with label ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
dir.create("./Plots_labels")
setwd("./Plots_labels")
getwd()

DefaultAssay(Round0) <- "RNA"

Idents(object = Round0) <- ("ClusterNames_0.8")
levels(Idents(object = Round0))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])



#for (i in 1:2){
for (i in 1:length(marker)){
  pdf(paste(marker[i],"_VlnPlot_Round0.pdf",sep=''), width=14,height=4)
  p <- VlnPlot(Round0, features=as.character(marker[i]),  split.by = "stim", cols = c("green", "red", "cyan"))
  print(p)
  dev.off()
  
  pdf(paste(marker[i],"_FeaturePlot_Round0.pdf",sep=''), width=18,height=6)
  #p <- FeaturePlot(object = Round0, split.by = "stim", features = as.character(marker[i]),min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
  p <- FeaturePlot(object = Round0, split.by = "stim", features = as.character(marker[i]),label=FALSE, order=TRUE) 
  print(p)
  dev.off()
}


setwd("../")

#################################################################################################
##Feature plot with label02 ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
dir.create("./Plots_labels_02")
setwd("./Plots_labels_02")
getwd()

DefaultAssay(Round0) <- "RNA"

Idents(object = Round0) <- ("ClusterNames_0.8")
levels(Idents(object = Round0))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])


#for (i in 1:2){
for (i in 1:length(marker)){
  pdf(paste(marker[i],"_VlnPlot_Round0.pdf",sep=''), width=14,height=4)
  p <- VlnPlot(Round0, features=as.character(marker[i])  )
  print(p)
  dev.off()
  
  pdf(paste(marker[i],"_FeaturePlot_Round0.pdf",sep=''), width=18,height=6)
  p <- FeaturePlot(object = Round0, split.by = "stim", features = as.character(marker[i]),min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
  #p <- FeaturePlot(object = Round0, split.by = "stim", features = as.character(marker[i]),label=FALSE, order=TRUE) 
  print(p)
  dev.off()
}


setwd("../")

#####################################################################################################
## DEG table and heatmap with labels ####
Round0@assays
colnames(Round0@meta.data)
DefaultAssay(Round0)
levels(Idents(object = Round0))

Idents(object = Round0) <- ("ClusterNames_0.8")
DefaultAssay(Round0) <- "RNA"

Combined.markers <- FindAllMarkers(object = Round0, only.pos = TRUE)


top1000 <- Combined.markers %>% group_by(cluster) %>% top_n(1000, avg_log2FC)
write.csv(top1000, file = "top1000_labeled.csv")

all_gene <- Combined.markers %>% group_by(cluster)
write.csv(all_gene, file = "all_gene_labeled.csv")



top10 <- Combined.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

Round0 <- ScaleData((Round0))
pdf('Diff_Heatmp_labeled.pdf', width=24, height=28)
DoHeatmap(object = Round0, features = top10$gene,raster = TRUE, group.bar = TRUE, draw.lines=TRUE) 
dev.off()

pdf('Diff_Heatmp_labeled_02.pdf', width=24, height=28)
DoHeatmap(object = Round0, features = top10$gene, raster = FALSE, size = 5, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)       + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
dev.off()

pdf('Diff_Heatmp_labeled_03.pdf', width=24, height=28)
DoHeatmap(object = Round0, features = top10$gene,raster = FALSE, group.bar = TRUE, draw.lines=TRUE) 
dev.off()

########################################################################
## Identify differential expressed genes across conditions ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

dir.create("./PsoriasisPre.Post.Control")
setwd("./PsoriasisPre.Post.Control")
getwd()

colnames(Round0@meta.data)
DefaultAssay(Round0) <- "RNA"
levels(Idents(object = Round0))
Idents(object = Round0) <- ("ClusterNames_0.8")

marker <- levels(Idents(object=Round0))


Round0$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(Round0), Round0$stim, sep = "_")
Idents(Round0) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(Round0) <- fct_relevel(Idents(Round0), sort)
Round0[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = Round0)
Idents(object = Round0) <- ("ClusterNames_0.8_PsoriasisvsControl")
levels(Idents(object = Round0))

#Psoriasis_PrevsControl
PsoriasisvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round0, ident.1 = paste(marker,"_Psoriasis_Pre",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round0_Psoriasis_PrevsControl.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(PsoriasisvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)



#Psoriasis_TxvsPsoriasis_Pre
PsoriasisTxvsPsoriasisPre.cell.marker <- function(marker){
  table <- FindMarkers(Round0, ident.1 = paste(marker,"_Psoriasis_Tx",sep=''), ident.2 = paste(marker,"_Psoriasis_Pre",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round0_Psoriasis_TxvsPsoriasis_Pre.csv",sep=''))
}

testFunction2 <- function (marker) {
  return(tryCatch(PsoriasisTxvsPsoriasisPre.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction2,USE.NAMES = TRUE, simplify = FALSE)

#Psoriasis_TxvsControl
PsoriasisTxvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round0, ident.1 = paste(marker,"_Psoriasis_Tx",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round0_Psoriasis_TxvsControl.csv",sep=''))
}

testFunction3 <- function (marker) {
  return(tryCatch(PsoriasisTxvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction3,USE.NAMES = TRUE, simplify = FALSE)


#Psoriasis_PrevsPsoriasis_NL
Psoriasis_PrevsPsoriasis_NL.cell.marker <- function(marker){
  table <- FindMarkers(Round0, ident.1 = paste(marker,"_Psoriasis_Pre",sep=''), ident.2 = paste(marker,"_Psoriasis_NL",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round0_Psoriasis_PrevsPsoriasis_NL.csv",sep=''))
}

testFunction4 <- function (marker) {
  return(tryCatch(Psoriasis_PrevsPsoriasis_NL.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction4,USE.NAMES = TRUE, simplify = FALSE)




#Psoriasis_NLvsControl
Psoriasis_NLvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round0, ident.1 = paste(marker,"_Psoriasis_NL",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round0_Psoriasis_NLvsControl.cell.marker.csv",sep=''))
}

testFunction5 <- function (marker) {
  return(tryCatch(Psoriasis_NLvsControl.cell.marker.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction5,USE.NAMES = TRUE, simplify = FALSE)


#Psoriasis_TxvsPsoriasis_NL
Psoriasis_TxvsPsoriasis_NL.cell.marker <- function(marker){
  table <- FindMarkers(Round0, ident.1 = paste(marker,"_Psoriasis_Tx",sep=''), ident.2 = paste(marker,"_Psoriasis_NL",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round0_Psoriasis_TxvsPsoriasis_NL.cell.marker.csv",sep=''))
}

testFunction6 <- function (marker) {
  return(tryCatch(Psoriasis_TxvsPsoriasis_NL.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction6,USE.NAMES = TRUE, simplify = FALSE)



setwd('../')


####################################################################################
##Cluster gene expression cells by sample and clusters ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")


Idents(Round0)
colnames(Round0@meta.data)
Idents(object = Round0) <- ("ClusterNames_0.8")
levels(Idents(object = Round0))

Round0$ClusterNames_0.8_Samples <- paste(Round0$number, Idents(Round0), sep = "_")

Idents(Round0) <- "ClusterNames_0.8_Samples"
Idents(Round0) <- fct_relevel(Idents(Round0), sort)
levels(Idents(object = Round0))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])

expr <- FetchData(object = Round0, vars = marker[1])
list <- Round0[, which(x = expr > 1)]
list <- as.data.frame(table(list@active.ident))
colnames(list) <- c("number",as.character(marker[1]))
list1 <- list

#for (i in 2:5){
for (i in 2:length(marker)){
  tryCatch(expr <- FetchData(object = Round0, vars = marker[i]), error=function(e) NULL)
  tryCatch(list <- Round0[, which(x = expr > 1)], error=function(e) NULL)
  tryCatch(list <- as.data.frame(table(list@active.ident)), error=function(e) NULL)
  colnames(list) <- c("number",as.character(marker[i]))
  list1 <- merge(list1,list, by="number", all = TRUE)
}

list1[is.na(list1)] <- 0

#list1 <- list1[,!grepl(".y", colnames(list1))]
#colnames(list1) = gsub(".x", "", colnames(list1))

cellnumbers <- as.data.frame(table(Round0@meta.data$ClusterNames_0.8_Samples))
cellnumbers[is.na(cellnumbers)] <- 0

colnames(cellnumbers) <- c("number","cell_number")

list2 <- merge(cellnumbers,list1, by="number", all = TRUE)
list2[is.na(list2)] <- 0
colnames(list2)[1] <- "Cluster_in_sample"



##percentage table
list3_classifier <- as.data.frame(list2[,1])
colnames(list3_classifier) <- "Cluster_in_sample"
list3_data <- list2[,-1]
list3_data <- data.frame(apply(list3_data, 2, function(x) as.numeric(as.character(x))))
list3_data <- list3_data / c(list3_data[,1])*100
colnames(list3_data) <- paste(colnames(list3_data), "_proportion", sep="")
list3_data <- list3_data[,-1]

list3 <- cbind(list3_classifier, list3_data)

list2 <- merge(list2, list3, by="Cluster_in_sample", all = TRUE)

colnames(list2)

write.csv(list2, file="gene_expressing_cells_sample.csv")

rm(list3)
rm(list3_classifier)
rm(list3_data)

####################################################################################
## Bar graph for samples ####
list <- list2[,c(1,2)]
list$Cluster_in_sample
list$Cluster <- list$Cluster_in_sample
levels(list$Cluster)


#list$Cluster= gsub("Dermis_", "", list$Cluster)
#list$Cluster= gsub("Epi_", "", list$Cluster)

list$Cluster= gsub(".*_LS_", "", list$Cluster)
list$Cluster= gsub(".*_NL_", "", list$Cluster)
list$Cluster= gsub("Control.._", "", list$Cluster)
#list$Cluster= gsub("repeat_", "", list$Cluster)
list$Cluster= gsub("02_", "", list$Cluster)

list$Cluster <- factor(list$Cluster)
levels(list$Cluster)

list$Cluster <- factor(list$Cluster, 
                       levels = (c(
                         "NK_cell","CD161_T_cell","CD8_T_cell","CD4_T_cell.01","CD4_T_cell.02","CD4_T_cell.03","Treg","Mature_DC.01","Mature_DC.02","Semimature_DC.01","Semimature_DC.02","Melanocyte","KC-S.Corneum.01","KC-S.Corneum.02","KC-S.Corneum.03","KC-S.Corneum.04","KC-S.Corneum.05","KC-S.Corneum.06","KC-S.Corneum.07","KC-S.Granulosum","KC-S.Spinosum","KC-S.Basale","KC-Hair_follicle","Endothelial_cell","Fibroblast","Eccrine_gland","Doublet.01","Doublet.02"
                                                )))

list$PsoriasisvsControl <- list$Cluster_in_sample
list$PsoriasisvsControl <- as.factor(list$PsoriasisvsControl)
levels(list$PsoriasisvsControl)

list$PsoriasisvsControl= gsub("_Eccrine_gland", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Endothelial_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Treg", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Semimature_DC...", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Corneum...", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Spinosum", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Basale", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Mature_DC...", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Fibroblast", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Melanocyte", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_NK_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_CD8_T_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_CD4_T_cell...", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_CD161_T_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-Hair_follicle", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Granulosum", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Doublet...", "", list$PsoriasisvsControl)

list$PsoriasisvsControl  <- factor(list$PsoriasisvsControl )
levels(list$PsoriasisvsControl)

list_for_ggplot <- melt(list)

colourCount = length(unique(list$PsoriasisvsControl))
getPalette = colorRampPalette(brewer.pal(12, "Accent"))

getPalette(colourCount)

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = PsoriasisvsControl)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  xlab("Samples") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual(values = getPalette(colourCount))


ggsave("Samples_in_cluster_bargraph_Round0.pdf", width = 12, height = 8)





####################################################################################
## Bar graph for control vs. preTx vs. postTx ####
list$simplified_category <- list$PsoriasisvsControl
list$simplified_category <- as.factor(list$simplified_category)
levels(list$simplified_category)

list$simplified_category= gsub("01", "", list$simplified_category)
list$simplified_category= gsub("_02", "", list$simplified_category)
list$simplified_category= gsub("02", "", list$simplified_category)
list$simplified_category= gsub("03", "", list$simplified_category)
list$simplified_category= gsub("04", "", list$simplified_category)
list$simplified_category= gsub("05", "", list$simplified_category)
list$simplified_category= gsub("06", "", list$simplified_category)
list$simplified_category= gsub("07", "", list$simplified_category)
list$simplified_category= gsub("08", "", list$simplified_category)
list$simplified_category= gsub("09", "", list$simplified_category)
list$simplified_category= gsub("10", "", list$simplified_category)


list$simplified_category <- factor(list$simplified_category, 
                       levels = (rev(c(
                         "Control"  ,
                         "Psoriasis_preTx_LS"     ,    "Psoriasis_preTx_NL"  ,
                         "Psoriasis_postTx_week12_LS", "Psoriasis_postTx_week24_LS", "Psoriasis_postTx_week48_LS"
                                                ))))
levels(list$simplified_category)

list_for_ggplot <- melt(list)

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = simplified_category)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  xlab("Samples") + ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual(values=rev(c('royalblue1','firebrick1','coral1','goldenrod1','khaki2','rosybrown3')))

ggsave("Cells_simplified_bargraph_Round0.pdf", width = 12, height = 8)

####################################################################################
##Cluster gene expression cells by psoriasis vs. control and clusters ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")


Idents(Round0)
colnames(Round0@meta.data)
Idents(object = Round0) <- ("stim")
levels(Idents(object = Round0))


Idents(object = Round0) <- ("ClusterNames_0.8")
levels(Idents(object = Round0))

Round0$ClusterNames_0.8_stim <- paste(Idents(Round0),Round0$stim, sep = "_")

Idents(Round0) <- "ClusterNames_0.8_stim"
Idents(Round0) <- fct_relevel(Idents(Round0), sort)
levels(Idents(object = Round0))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])

expr <- FetchData(object = Round0, vars = marker[1])
list <- Round0[, which(x = expr > 1)]
list <- as.data.frame(table(list@active.ident))
colnames(list) <- c("number",as.character(marker[1]))
list1 <- list

#for (i in 2:100){
for (i in 2:length(marker)){
  tryCatch(expr <- FetchData(object = Round0, vars = marker[i]), error=function(e) NULL)
  tryCatch(list <- Round0[, which(x = expr > 1)], error=function(e) NULL)
  tryCatch(list <- as.data.frame(table(list@active.ident)), error=function(e) NULL)
  colnames(list) <- c("number",as.character(marker[i]))
  list1 <- merge(list1,list, by="number", all = TRUE)
}

list1[is.na(list1)] <- 0

#list1 <- list1[,!grepl(".y", colnames(list1))]
#colnames(list1) = gsub(".x", "", colnames(list1))

cellnumbers <- as.data.frame(table(Round0@meta.data$ClusterNames_0.8_stim))
cellnumbers[is.na(cellnumbers)] <- 0

colnames(cellnumbers) <- c("number","cell_number")

list2 <- merge(cellnumbers,list1, by="number", all = TRUE)
list2[is.na(list2)] <- 0
colnames(list2)[1] <- "Cluster_in_sample"



##percentage table
list3_classifier <- as.data.frame(list2[,1])
colnames(list3_classifier) <- "Cluster_in_sample"
list3_data <- list2[,-1]
list3_data <- data.frame(apply(list3_data, 2, function(x) as.numeric(as.character(x))))
list3_data <- list3_data / c(list3_data[,1])*100
colnames(list3_data) <- paste(colnames(list3_data), "_proportion", sep="")
list3_data <- list3_data[,-1]

list3 <- cbind(list3_classifier, list3_data)

list2 <- merge(list2, list3, by="Cluster_in_sample", all = TRUE)

colnames(list2)

write.csv(list2, file="gene_expressing_cells_stim.csv")

rm(list3)
rm(list3_classifier)
rm(list3_data)#

#####################################################################################################
## Save
save(Round0, file = "~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_5.31.2021.Rda")
