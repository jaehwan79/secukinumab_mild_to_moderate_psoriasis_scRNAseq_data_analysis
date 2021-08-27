#load library ####
#instRound0.packages('Seurat')
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
#sourceters
#Start ####("https://raw.githubusercontent.com/farrellja/URD/master/URD-InstRound0.R")
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
#install.packages("viridis")
library(viridis)
#install.packages('data.table')
library(data.table) 
library(googlesheets4)
######################################################################################################3
# Round1 combining clusters ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021")



dir.create("./Round1")
setwd("./Round1")
getwd()
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

load("/Users/jaehwankim/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round1_integrated_analyzed_6.1.2021.Rda")

list.files()
#rm(list= ls()[!(ls() %in% c('Round0'))])
rm(list= ls()[!(ls() %in% c('Round1'))])

gene_subset  <- read_sheet("https://docs.google.com/spreadsheets/d/19nhrzGudrn8ihl4RASuyGlBe7JMHdWnnWKoZ0i2rZ4s/edit#gid=0")
gene_subset <- as.matrix(unique(gene_subset[1]))

#load("/Users/jkim05/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_01.31.2021.Rda")
#load("/Users/jkim05/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_01.31.2021.Rda")
#load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_02.09.2020.Rda")
#load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1/combined_data/Round1_integrated_analyzed_02.15.2020.Rda")
######################################################################################################
## (SKIP) Combine clustes ####
colnames(Round0@meta.data)
DefaultAssay(Round0) <- "RNA"
colnames(Round0@meta.data)
Idents(object = Round0) <- ("orig.ident")
Idents(object = Round0) <- ("number")
Idents(object = Round0) <- ("stim")

levels(Round0)

phenotype <- read_sheet("https://https://docs.google.com/spreadsheets/d/1_cVcmP7AniLrhIE9XWw-T-_BqQXIkG8wi7wl1Rb3kl4/edit#gid=0&fvid=109363522")
phenotype <- phenotype[grep("Yes", phenotype$Third_project), ]

#Serial
Idents(object = Round0) <- ("number")
Idents(object = Round0) <- factor(Idents(object = Round0), levels = phenotype$Rcoding_Number_05.01.2021)

new.cluster.ids <- phenotype$Rcoding_PrevsPost_05.01.2021
names(new.cluster.ids) <- levels(Round0)
Round0 <- RenameIdents(Round0, new.cluster.ids)

levels(Idents(Round0))

Idents(object = Round0) <- factor(Idents(object = Round0), 
                                  levels = (c("Control"   ,
                                              "Psoriasis_PreTx_LS_week0" ,
                                              "Psoriasis_PreTx_NL_week0"  ,
                                              "Psoriasis_PostTx_week12"  , 
                                              "Psoriasis_PostTx_LS_week24",
                                              "Psoriasis_PostTx_LS_week48"
                                              )))
Round0[["Serial"]] <- Idents(object = Round0)


#Round1
Round1 <- Round0
colnames(Round1@meta.data)
DefaultAssay(Round1) <- "RNA"
Idents(object = Round1) <- ("ClusterNames_0.8")
levels(Round1)

Round1 <- subset(Round1, idents = c( "CD4_T_cell.02" , "Doublet.01"   ,    "Doublet.02"),invert=TRUE)
levels(Round1)

new.cluster.ids <- c(
  "NK_cell"     ,
  "CD161_T_cell"   ,
  "CD8_T_cell"    ,
  "CD4_T_cell"   ,
  "CD4_T_cell"  ,
  "Treg"     ,       
  "Mature_DC"   ,
  "Mature_DC"  ,
  "Semimature_DC" ,
  "Semimature_DC" ,
  "Melanocyte"    ,
  "KC-S.Corneum" ,
  "KC-S.Corneum" ,
  "KC-S.Corneum" ,
  "KC-S.Corneum" ,
  "KC-S.Corneum" ,
  "KC-S.Corneum",
  "KC-S.Corneum" ,
  "KC-S.Granulosum" ,
  "KC-S.Spinosum"   ,
   "KC-S.Basale"   ,
  "KC-Hair_follicle" ,
  "Endothelial_cell",
  "Fibroblast"    ,
  "Eccrine_gland" 
)	


names(new.cluster.ids) <- levels(Round1)
Round1 <- RenameIdents(Round1, new.cluster.ids)
levels(Idents(Round1))
Round1[["ClusterNames_0.8"]] <- Idents(object = Round1)

######################################################################################################
## Save and Dimensionality reduction ####
getwd()
dir.create("./combined_data")
setwd("./combined_data")

## Save
#save(Round1, file = "Round1_integrated_analyzed_5.01.2021.Rda")

## ElbowPlot 
pdf('ElbowPlot_Round1.pdf',width=10,height=10)
ElbowPlot(object = Round1, ndims = 20)
dev.off()

## JackStraw plot
#Round1 <- JackStraw(object = Round1, num.replicate = 100)
#Round1 <- ScoreJackStraw(object = Round1, dims = 1:20)

#pdf('JackStrawPlot_Round1.pdf',width=10,height=10)
#JackStrawPlot(object = Round1, dims = 1:20)
#dev.off()

## (SKIP) Perform linear dimensional reduction
pdf('VizDimLoadings_Round1.pdf', width=8,height=7)
VizDimLoadings(Round1, dims = 1:2, reduction = "pca")
dev.off()

pdf('DimPlots_Round1.pdf', width=8,height=7)
DimPlot(Round1, reduction = "pca")
dev.off()

setwd("../")

########################################################################
# (From here) Dot plot with labels ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Round1@meta.data$ClusterNames_0.8
DefaultAssay(Round1)
levels(Idents(object = Round1))

Idents(object = Round1) <- ("ClusterNames_0.8")
#Idents(object = Round1) <-"ClusterNames_0.8_PsoriasisvsControl"
#Idents(object = Round1) <-"ClusterNames_0.8"

DefaultAssay(Round1) <- "RNA"

features.plot <- c(  "KLRB1"   , "GNLY" ,  
                     "CD3D" ,    "TRAC" ,    "TRBC1" ,
                     "CD8A"   ,  "CD8B"   ,
                     "IFNG","IL17A","IL17F",
                     "GZMK"    ,"GZMH"   ,
                     # "PTPRC"  ,  
                     #"IGLC2" ,  
                     "TIGIT" ,   "IL2RA"   ,  "FOXP3"   , "CTLA4" ,  
                     "LAMP3"  ,  "CCL22"  , "LY75"  ,   "CIITA"  ,  "CD40"  ,  
                     "HLA-DQA1","HLA-DQB1", "HLA-DRB1", "HLA-DRA" , "HLA-DRB5" ,
                     "LYZ"     ,  "CD14"   ,  "CD163"   ,  "THBD"    ,"IL10",
                    # "CLEC10A" ,                      "CD1C"   , 
                     "DCT"    ,  "TYRP1"  ,  "MLANA"   ,
                     #"CDSN" ,
                     #"KRT80"   , "KRT23"   ,"SPINK5" , "PRSS22"  ,
                     "SPRR2A"  , "SPRR2D"  , 
                    # "APOBEC3A", "IL1RL1" ,"IL18"   ,  
                    # "LCE2C" ,     "LCE1A"  , 
                     "FABP5"  ,
                     "KRT10"   , "KRT1"   ,"CCL20"  ,  
                     "KRT14"  ,  "KRT5"  ,
                     "KRT6B"   , "KRT17"   , "KRT16"   ,
                     "CCL21"  ,  "LYVE1"   , 
                    #"CD34"  ,  
                    "DCN"    ,  "COL1A1" ,  "DCD"   ,   "SCGB2A2" )




Idents(object = Round1) <- fct_rev(Idents(object = Round1))


pdf('Dot plot_Round1_labeled01.pdf', width=14,height=5)
DotPlot(object = Round1, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()


pdf('Dot plot_Round1_labeled02.pdf', width=14,height=5)
DotPlot(object = Round1, features = features.plot, cols = c("white","red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()

pdf('Dot plot_Round1_split_labeled.pdf', width=18,height=20)
DotPlot(Round1, features = features.plot, cols =c("black", "pink","red","blue"),  
        split.by = "stim") + RotatedAxis()+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=18))
dev.off()

Idents(object = Round1) <- fct_rev(Idents(object = Round1))

#######################################################################
# Visualization with label ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

dir.create("./UMAP_Label")
setwd("./UMAP_Label")

colnames(Round1@meta.data)
DefaultAssay(Round1)
DefaultAssay(Round1) <- "RNA"

Idents(object = Round1) <- ("Serial")
levels(Round1)
#Idents(object = Round1) <- factor(Idents(object = Round0), levels = c("Control","Psoriasis_NL","Psoriasis_Pre","Psoriasis"))       
Idents(object = Round1) <- ("ClusterNames_0.8")
Idents(object = Round1) 

color_code <- c('orange','gold4','firebrick1','salmon', 'mediumorchid1', 'green4','cyan3','indianred4','deepskyblue','royalblue1','navy','deeppink3','deepskyblue4','hotpink4','snow4','snow3')

##DimPlot
Idents(object = Round1) <- ("ClusterNames_0.8")

p1 <- DimPlot(Round1, reduction = "umap", group.by = "Serial")+
scale_color_manual(values=c('royalblue1','firebrick1','coral1','goldenrod1','khaki2','rosybrown3'))

p2 <- DimPlot(Round1, reduction = "umap", label = FALSE, cols=color_code)

pdf('UMAP_DimPlot_Round1_ClusterNames_0.8.pdf', width=16,height=7)
plot_grid(p1, p2)
dev.off()


##cluster count 
Idents(object = Round1) <- ("ClusterNames_0.8")
Idents(object = Round1) 
colnames(Round1@meta.data)
cell.num <- table(Round1@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

levels(Round1)

pdf('UMAP_Round1_cluster_count_ClusterNames_0.8.pdf', width=10,height=8)
DimPlot(object = Round1,label = TRUE, label.size = 4, reduction = "umap") +
  scale_colour_manual(values = color_code, breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()


###split by Serial
pdf('UMAP_Round1_DimPlot_split_ClusterNames_0.8.pdf', width=24,height=7)
DimPlot(Round1, reduction = "umap", split.by = "Serial", cols=color_code)
dev.off()


##Cell number
colnames(Round1@meta.data)
Idents(object = Round1) 
Idents(object = Round1) <- ("number")

cell.num <- table(Round1@meta.data$number)
ClusterLabels = paste(names(cell.num), paste0("(n=", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round1_cellcount_ClusterNames_0.8.pdf', width=13,height=8)
DimPlot(object = Round1,label = FALSE, label.size = 5, reduction = "umap") +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##control vs psoriasis cell count
colnames(Round1@meta.data)
Idents(object = Round1) 

Idents(object = Round1) <- ("Serial")
levels(Round1)
cell.numb <- table(Round1@meta.data$Serial)
ClusterLabels = paste(names(cell.numb), paste0("(n = ", cell.numb, ")"))
ClusterBreaks = names(cell.numb)

pdf('UMAP_Round1_category_ClusterNames_0.8.pdf',width=10,height=8)
DimPlot(object = Round1,label = FALSE, label.size = 5, reduction = "umap")+
  scale_colour_manual(values =c('royalblue1','firebrick1','coral1','goldenrod1','khaki2','rosybrown3'), breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")
dev.off()



##cluster count _modified
Idents(object = Round1) <- ("ClusterNames_0.8")
Idents(object = Round1) 
colnames(Round1@meta.data)
cell.num <- table(Round1@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round1_cluster_count_ClusterNames_0.8.pdf', width=10, height=8)
DimPlot(object = Round1,label = TRUE,  reduction = "umap") +
  scale_colour_manual(values = color_code, breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")


dev.off()

##cluster count _modified
Idents(object = Round1) <- ("ClusterNames_0.8")
Idents(object = Round1) 
colnames(Round1@meta.data)
cell.num <- table(Round1@meta.data$ClusterNames_0.8)
ClusterLabels = paste(names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)

pdf('UMAP_Round1_cluster_count_ClusterNames_0.8_nolabel.pdf', width=10,height=8)
DimPlot(object = Round1,label = FALSE,  reduction = "umap") +
  scale_colour_manual(values = color_code, breaks = ClusterBreaks, 
                      labels = ClusterLabels) +
  labs(x = "UMAP 1",
       y = "UMAP 2")


dev.off()


setwd("../")



#save(Round1, file = "./combined_data/Round1_integrated_analyzed_05.31.2021.Rda")

#################################################################################################
##Feature plot ####
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
dir.create("./Plots")
setwd("./Plots")
getwd()

DefaultAssay(Round1) <- "RNA"
Idents(object = Round1) <- ("ClusterNames_0.8")
levels(Idents(object = Round1))

#Idents(object = Round1) <- ("Serial")
#Idents(object = Round1) <- ("Serial")

marker <- gene_subset

for (i in 1:length(marker)){
  pdf(paste(marker[i],"_VlnPlot_Round1.pdf",sep=''), width=14,height=4)
  p <- VlnPlot(Round1, features=as.character(marker[i]),  split.by = "Serial", cols = c('royalblue1','firebrick1','coral1','goldenrod1','khaki2','rosybrown3'))
  print(p)
  dev.off()
  
  pdf(paste(marker[i],"_FeaturePlot_Round1.pdf",sep=''), width=18,height=6)
  #p <- FeaturePlot(object = Round1, split.by = "stim", features = as.character(marker[i]),min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
  p <- FeaturePlot(object = Round1, split.by = "Serial", features = as.character(marker[i]),label=FALSE, order=TRUE) 
  print(p)
  dev.off()
}


setwd("../")


#################################################################################################
##Feature plot with label02 ####
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
dir.create("./Plots_labels_02")
setwd("./Plots_labels_02")
getwd()

DefaultAssay(Round1) <- "RNA"

Idents(object = Round1) <- ("ClusterNames_0.8")
levels(Idents(object = Round1))

marker <- gene_subset

#for (i in 1:2){
for (i in 1:length(marker)){
  pdf(paste(marker[i],"_VlnPlot_Round1.pdf",sep=''), width=10,height=4)
  p <- VlnPlot(Round1, features=as.character(marker[i])  )
  print(p)
  dev.off()
  
  pdf(paste(marker[i],"_FeaturePlot_Round1.pdf",sep=''), width=18,height=6)
  p <- FeaturePlot(object = Round1, split.by = "Serial", features = as.character(marker[i]),min.cutoff = "q10", max.cutoff = "q90",label=FALSE, order=TRUE) 
  #p <- FeaturePlot(object = Round1, split.by = "stim", features = as.character(marker[i]),label=FALSE, order=TRUE) 
  print(p)
  dev.off()
}


setwd("../")

####################################################################################
##Cluster gene expression cells by sample ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Idents(Round1)
colnames(Round1@meta.data)
Idents(object = Round1) <- ("ClusterNames_0.8")
levels(Idents(object = Round1))

Round1$ClusterNames_0.8_Samples <- paste(Round1$number, Idents(Round1), sep = "_")

Idents(Round1) <- "ClusterNames_0.8_Samples"
Idents(Round1) <- fct_relevel(Idents(Round1), sort)
levels(Idents(object = Round1))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])

expr <- FetchData(object = Round1, vars = marker[1])
list <- Round1[, which(x = expr > 1)]
list <- as.data.frame(table(list@active.ident))
colnames(list) <- c("number",as.character(marker[1]))
list1 <- list

#for (i in 2:5){
for (i in 2:length(marker)){
  tryCatch(expr <- FetchData(object = Round1, vars = marker[i]), error=function(e) NULL)
  tryCatch(list <- Round1[, which(x = expr > 1)], error=function(e) NULL)
  tryCatch(list <- as.data.frame(table(list@active.ident)), error=function(e) NULL)
  colnames(list) <- c("number",as.character(marker[i]))
  list1 <- merge(list1,list, by="number", all = TRUE)
}

list1[is.na(list1)] <- 0

#list1 <- list1[,!grepl(".y", colnames(list1))]
#colnames(list1) = gsub(".x", "", colnames(list1))

cellnumbers <- as.data.frame(table(Round1@meta.data$ClusterNames_0.8_Samples))
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


## Bar graph for samples
list <- list2[,c(1,2)]
list$Cluster_in_sample
list$Cluster <- list$Cluster_in_sample
levels(list$Cluster)

list$Cluster= gsub(".*_LS_", "", list$Cluster)
list$Cluster= gsub(".*_NL_", "", list$Cluster)
list$Cluster= gsub("Control.._", "", list$Cluster)
#list$Cluster= gsub("repeat_", "", list$Cluster)
list$Cluster= gsub("02_", "", list$Cluster)

list$Cluster <- factor(list$Cluster)
levels(list$Cluster)

list$Cluster <- factor(list$Cluster, 
                       levels = (c(
                         "NK_cell","CD161_T_cell","CD8_T_cell","CD4_T_cell","Treg","Mature_DC","Semimature_DC",
                         "Melanocyte","KC-S.Corneum","KC-S.Granulosum","KC-S.Spinosum","KC-S.Basale","KC-Hair_follicle",
                         "Endothelial_cell","Fibroblast","Eccrine_gland"
                       )))

list$PsoriasisvsControl <- list$Cluster_in_sample
list$PsoriasisvsControl <- as.factor(list$PsoriasisvsControl)
levels(list$PsoriasisvsControl)

list$PsoriasisvsControl= gsub("_Eccrine_gland", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Endothelial_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Treg", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Semimature_DC", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Corneum", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Spinosum", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Basale", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Mature_DC", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Fibroblast", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_Melanocyte", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_NK_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_CD8_T_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_CD4_T_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_CD161_T_cell", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-Hair_follicle", "", list$PsoriasisvsControl)
list$PsoriasisvsControl= gsub("_KC-S.Granulosum", "", list$PsoriasisvsControl)

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

ggsave("Samples_in_cluster_bargraph_Round1.pdf", width = 12, height = 8)



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

ggsave("Cells_simplified_bargraph_Round1.pdf", width = 12, height = 8)

####################################################################################
##Cluster gene expression cells by psoriasis vs. control and clusters ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")


Idents(Round1)
colnames(Round1@meta.data)
Idents(object = Round1) <- ("Serial")
levels(Idents(object = Round1))


Idents(object = Round1) <- ("ClusterNames_0.8")
levels(Idents(object = Round1))

Round1$ClusterNames_0.8_Serial <- paste(Idents(Round1),Round1$Serial, sep = "_")

Idents(Round1) <- "ClusterNames_0.8_Serial"
Idents(Round1) <- fct_relevel(Idents(Round1), sort)
levels(Idents(object = Round1))

marker <- read.csv("~/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- read.csv("D:/Dropbox/SIngle_cell/Singlecell_Ref_Data/gene_subset.csv", header=TRUE)
marker <- unique(marker[,2])

expr <- FetchData(object = Round1, vars = marker[1])
list <- Round1[, which(x = expr > 1)]
list <- as.data.frame(table(list@active.ident))
colnames(list) <- c("number",as.character(marker[1]))
list1 <- list

#for (i in 2:5){
for (i in 2:length(marker)){
  tryCatch(expr <- FetchData(object = Round1, vars = marker[i]), error=function(e) NULL)
  tryCatch(list <- Round1[, which(x = expr > 1)], error=function(e) NULL)
  tryCatch(list <- as.data.frame(table(list@active.ident)), error=function(e) NULL)
  colnames(list) <- c("number",as.character(marker[i]))
  list1 <- merge(list1,list, by="number", all = TRUE)
}

list1[is.na(list1)] <- 0

#list1 <- list1[,!grepl(".y", colnames(list1))]
#colnames(list1) = gsub(".x", "", colnames(list1))

cellnumbers <- as.data.frame(table(Round1@meta.data$ClusterNames_0.8_Serial))
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

write.csv(list2, file="gene_expressing_cells_Serial.csv")

rm(list3)
rm(list3_classifier)
rm(list3_data)#
###################################################################################################
## DEG table and heatmap with labels ####
Round1@assays
colnames(Round1@meta.data)
DefaultAssay(Round1)
levels(Idents(object = Round1))

Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

Combined.markers <- FindAllMarkers(object = Round1, only.pos = TRUE, min.pct = 0, 
                                   logfc.threshold = 0.1)

top1000 <- Combined.markers %>% group_by(cluster) 

#top1000 <- Combined.markers %>% group_by(cluster) %>% top_n(1000, avg_logFC)
write.csv(top1000, file = "top1000_labeled.csv")


top10 <- Combined.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)



pdf('Diff_Heatmp_labeled.pdf', width=24, height=28)
DoHeatmap(object = Round1, features = top10$gene,raster = TRUE, group.bar = TRUE, draw.lines=TRUE) 
dev.off()

pdf('Diff_Heatmp_labeled_02.pdf', width=24, height=28)
DoHeatmap(object = Round1, features = top10$gene, raster = FALSE, size = 5, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)       + 
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")
dev.off()

pdf('Diff_Heatmp_labeled_03.pdf', width=24, height=28)
DoHeatmap(object = Round1, features = top10$gene,raster = FALSE, group.bar = TRUE, draw.lines=TRUE) 
dev.off()


########################################################################
## Identify differential expressed genes across conditions ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
dir.create("./PsoriasisPre.Post.Control")
setwd("./PsoriasisPre.Post.Control")
getwd()

colnames(Round1@meta.data)
DefaultAssay(Round1) <- "RNA"
levels(Idents(object = Round1))
Idents(object = Round1) <- ("ClusterNames_0.8")

marker <- levels(Idents(object=Round1))


Round1$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(Round1), Round1$Serial, sep = "_")
Idents(Round1) <- "ClusterNames_0.8_PsoriasisvsControl"
Idents(Round1) <- fct_relevel(Idents(Round1), sort)
Round1[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = Round1)
Idents(object = Round1) <- ("ClusterNames_0.8_PsoriasisvsControl")
levels(Idents(object = Round1))

# _Control
# _Psoriasis_PreTx_LS_week0
# _Psoriasis_PreTx_NL_week0
# _Psoriasis_PostTx_week12
# _Psoriasis_PostTx_LS_week24

#Psoriasis_PreTx_LS_week0 vs Control
PsoriasisvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PreTx_LS_week0",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PreTx_LS_week0vsControl.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(PsoriasisvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)


# _Psoriasis_PreTx_NL_week0 vs Control
PsoriasisvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PreTx_NL_week0",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PreTx_NL_week0vsControl.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(PsoriasisvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)


# _Psoriasis_PostTx_week12 vs Control
PsoriasisvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PostTx_week12",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PostTx_week12vsControl.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(PsoriasisvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)


# _Psoriasis_PostTx_LS_week24 vs Control
PsoriasisvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PostTx_LS_week24",sep=''), ident.2 = paste(marker,"_Control",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PostTx_LS_week24vsControl.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(PsoriasisvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)


# _Psoriasis_PreTx_LS_week0 vs _Psoriasis_PreTx_NL_week0
PsoriasisvsControl.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PreTx_LS_week0",sep=''), ident.2 = paste(marker,"_Psoriasis_PreTx_NL_week0",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PreTx_LS_week0_vs_Psoriasis_PreTx_NL_week0.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(PsoriasisvsControl.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)


##NOT WORKING

# _Psoriasis_PostTx_week12 vs _Psoriasis_PreTx_LS_week0
Psoriasisvs_Psoriasis_PreTx_LS_week0.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PostTx_week12",sep=''), ident.2 = paste(marker,"__Psoriasis_PreTx_LS_week0",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PostTx_week12vs_Psoriasis_PreTx_LS_week0.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(Psoriasisvs_Psoriasis_PreTx_LS_week0.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)



# _Psoriasis_PostTx_LS_week24 vs _Psoriasis_PreTx_LS_week0
Psoriasisvs_Psoriasis_PreTx_LS_week0.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PostTx_LS_week24",sep=''), ident.2 = paste(marker,"__Psoriasis_PreTx_LS_week0",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PostTx_LS_week24vs_Psoriasis_PreTx_LS_week0.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(Psoriasisvs_Psoriasis_PreTx_LS_week0.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)



# _Psoriasis_PostTx_LS_week24 vs _Psoriasis_PostTx_week12
Psoriasisvs_Psoriasis_PostTx_week12.cell.marker <- function(marker){
  table <- FindMarkers(Round1, ident.1 = paste(marker,"_Psoriasis_PostTx_LS_week24",sep=''), ident.2 = paste(marker,"__Psoriasis_PostTx_week12",sep=''), verbose = FALSE)
  write.csv(table, file=paste(marker,"_Round1_Psoriasis_PostTx_LS_week24vs_Psoriasis_PostTx_week12.csv",sep=''))
}

testFunction <- function (marker) {
  return(tryCatch(Psoriasisvs_Psoriasis_PostTx_week12.cell.marker(marker), error=function(e) NULL))
}

sapply(marker,testFunction,USE.NAMES = TRUE, simplify = FALSE)



setwd('../')


###################################################################################################
# Average expression - Dendritic cells ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Idents(object = Round1) <- ("ClusterNames_0.8")
DCs <- subset(Round1, idents = c("Mature_DC" ,   "Semimature_DC" ))

#Remove Non-responder
colnames(DCs@meta.data)
Idents(object = DCs) <- ("number")
levels(DCs)
DCs <- subset(DCs, idents = c("Psoriasis07_postTx_week12_LS"  ), inv=TRUE)

#Control vs Psoriasis_PreTx_LS_week0 vs Psoriasis_PostTx_week12
Idents(object = DCs) <- ("Serial")

levels(DCs)
DCs <- subset(DCs, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))

marker <- c(  "HLA-DRA" ,  "HLA-DPB1" ,"LAMP3"   ,  "CIITA" ,      "LY75"   ,   "CD86"   ,    "CD40"   , 
              "CD274"   ,  "PDCD1LG2" ,
              "CLEC4C",    "CD1C" ,  
              "ITGAX",  "AIF1"    ,"CD14"   ,      "THBD"  ,  
              "SIRPA",    "LILRB2" ,    
              "CLEC4A",
              "LILRB4"   ,  "LILRB1"  ,   
              "IL10","TGFB1",
              "FCGR1A", "IL23A"
             )


marker_shortened <- c(  "HLA-DPB1" ,"LAMP3"   ,
                        "ITGAX",
                         "THBD"  , "CLEC4A",
                        "IL10","TGFB1",
                        "IL23A"
)


Idents(object = DCs) <- ("ClusterNames_0.8")
DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$Serial, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(DCs)
Idents(DCs) <- factor(Idents(DCs), levels = c("Mature_DC_Control"  ,
                                              "Mature_DC_Psoriasis_PreTx_LS_week0"    ,
                                              "Mature_DC_Psoriasis_PostTx_week12"    ,
                                              "Semimature_DC_Control"    ,
                                             "Semimature_DC_Psoriasis_PreTx_LS_week0" ,
                                             "Semimature_DC_Psoriasis_PostTx_week12"
))

DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")


DCs.cluster.averages <- AverageExpression(DCs, return.seurat = TRUE)



pdf('DCs.averages_Round1_Heatmap.pdf', width=6, height=6)
DoHeatmap(object = DCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('DCs.averages_Round1_Heatmap_shortened.pdf', width=6, height=2)
DoHeatmap(object = DCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 2.5, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# Average expression - Dendritic cells : non-responder serial ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Idents(object = Round1) <- ("ClusterNames_0.8")
DCs <- subset(Round1, idents = c("Mature_DC" ,   "Semimature_DC" ))

# Non-responder
colnames(DCs@meta.data)
Idents(object = DCs) <- ("number")
levels(DCs)
DCs <- subset(DCs, idents = c(  "Psoriasis07_preTx_LS_02",
                              "Psoriasis07_postTx_week12_LS" ,
                              "Psoriasis07_postTx_week48_LS" 
))

#Control vs Psoriasis07 serial
Idents(object = DCs) <- ("Serial")
levels(DCs)
DCs <- subset(DCs, idents = c(
  "Psoriasis_PreTx_LS_week0" ,
  "Psoriasis_PostTx_week12"   ,
  "Psoriasis_PostTx_LS_week48"  ))


Idents(object = DCs) <- ("ClusterNames_0.8")
DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$Serial, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(DCs)
Idents(DCs) <- factor(Idents(DCs), levels = c(
  "Mature_DC_Psoriasis_PreTx_LS_week0"  ,
  "Mature_DC_Psoriasis_PostTx_week12"    ,
  "Mature_DC_Psoriasis_PostTx_LS_week48" ,
  
  "Semimature_DC_Psoriasis_PreTx_LS_week0"  ,
  "Semimature_DC_Psoriasis_PostTx_week12"  , 
  "Semimature_DC_Psoriasis_PostTx_LS_week48" 
))

DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")


DCs.cluster.averages <- AverageExpression(DCs, return.seurat = TRUE)

pdf('DCs.averages_Psoriasis07_nonresponder_Heatmap.pdf', width=6, height=6)
DoHeatmap(object = DCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('DCs.averages_Psoriasis07_nonresponder_Heatmap_shortened.pdf', width=6, height=2)
DoHeatmap(object = DCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('DCs.averages_Psoriasis07_nonresponder_Heatmap_shortened_expressionbar.pdf', width=12, height=4)
DoHeatmap(object = DCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# Average expression - Dendritic cells : responder serial ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Idents(object = Round1) <- ("ClusterNames_0.8")
DCs <- subset(Round1, idents = c("Mature_DC" ,   "Semimature_DC" ))

#Remove Non-responder
colnames(DCs@meta.data)
Idents(object = DCs) <- ("number")
levels(DCs)
DCs <- subset(DCs, idents = c(
  "Psoriasis08_preTx_LS"    ,
  "Psoriasis08_postTx_week12_LS",
  "Psoriasis08_postTx_week24_LS"
))

#Control vs Psoriasis07 serial
Idents(object = DCs) <- ("Serial")
levels(DCs)
DCs <- subset(DCs, idents = c(
  "Psoriasis_PreTx_LS_week0" ,
  "Psoriasis_PostTx_week12" ,
  "Psoriasis_PostTx_LS_week24" ))


Idents(object = DCs) <- ("ClusterNames_0.8")
DCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(DCs), DCs$Serial, sep = "_")
Idents(DCs) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(DCs)
Idents(DCs) <- factor(Idents(DCs), levels = c(
  "Mature_DC_Psoriasis_PreTx_LS_week0"   ,  
  "Mature_DC_Psoriasis_PostTx_week12"   ,
  "Mature_DC_Psoriasis_PostTx_LS_week24" ,
  
  "Semimature_DC_Psoriasis_PreTx_LS_week0" ,
  "Semimature_DC_Psoriasis_PostTx_week12"   ,
  "Semimature_DC_Psoriasis_PostTx_LS_week24"
))

DCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = DCs)
Idents(object = DCs) <- ("ClusterNames_0.8_PsoriasisvsControl")


DCs.cluster.averages <- AverageExpression(DCs, return.seurat = TRUE)

pdf('DCs.averages_Psoriasis08_responder_Heatmp.pdf', width=6, height=6)
DoHeatmap(object = DCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()



pdf('DCs.averages_Psoriasis08_responder_Heatmap_shortened.pdf', width=6, height=2)
DoHeatmap(object = DCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

pdf('DCs.averages_Psoriasis08_responder_Heatmap_shortened_expressionbar.pdf', width=12, height=4)
DoHeatmap(object = DCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


###################################################################################################
# Average expression - Keratinocytes ####

Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

memory.limit()
memory.limit(size=60000)

levels(Round1)
KCs <- subset(Round1, idents = c( c("KC-S.Corneum",  "KC-S.Granulosum" ,  "KC-S.Spinosum"  ,  "KC-S.Basale"  )))

#Control vs Psoriasis_PreTx_LS_week0 vs Psoriasis_PostTx_week12
Idents(object = KCs) <- ("Serial")
levels(KCs)
KCs <- subset(KCs, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))

#Remove Non-responder
colnames(KCs@meta.data)
Idents(object = KCs) <- ("number")
levels(KCs)
KCs <- subset(KCs, idents = c("Psoriasis07_postTx_week12_LS"  ), inv=TRUE)

Idents(object = KCs) <- ("Serial")
levels(KCs)

marker <- c("LCE3D", "CDSN" ,    
             "SPRR2G", "FABP5"  ,
            "S100A8","DEFB4B",
            "DEFB4A",
            "IL36G",
            #"IL17C",
              
            "KRT10" , 
            "KRT1",  
            "KRT5" , "KRT14" , "KRT15" 
                )


marker_shortened <- c(  "LCE3D", 
                        "IL36G",
                        "FABP5"  ,
                          "KRT10" , 
                        #"KRT1",  
                        "KRT14" ,
                        "IL36G"
)


Idents(object = KCs) <- ("ClusterNames_0.8")

KCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(KCs), KCs$Serial, sep = "_")
Idents(KCs) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(KCs)


Idents(object = KCs) <- factor(Idents(object = KCs), levels = (c("KC-S.Corneum_Control" ,
                                                                 "KC-S.Corneum_Psoriasis_PreTx_LS_week0"  ,  
                                                                 "KC-S.Corneum_Psoriasis_PostTx_week12" ,
                                                                 "KC-S.Granulosum_Control"   ,
                                                                 "KC-S.Granulosum_Psoriasis_PreTx_LS_week0" ,
                                                                 "KC-S.Granulosum_Psoriasis_PostTx_week12",
                                                                 "KC-S.Spinosum_Control"  ,
                                                                 "KC-S.Spinosum_Psoriasis_PreTx_LS_week0"  , 
                                                                 "KC-S.Spinosum_Psoriasis_PostTx_week12" ,
                                                                 "KC-S.Basale_Control"   ,                   
                                                                 "KC-S.Basale_Psoriasis_PreTx_LS_week0"  ,
                                                                 "KC-S.Basale_Psoriasis_PostTx_week12"
  )))
KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)

KCs.cluster.averages <- AverageExpression(KCs, return.seurat = TRUE)



pdf('KCs.cluster.averages_Round1_Heatmp01.pdf', width=8, height=3.2)
DoHeatmap(object = KCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

pdf('KCs.cluster.averages_Round1_Heatmp02.pdf', width=8, height=5)
DoHeatmap(object = KCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

pdf('KCs.cluster.averages_Round1_Heatmp01_shortened.pdf', width=8, height=2)
DoHeatmap(object = KCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

pdf('KCs.cluster.averages_Round1_Heatmp01_shortened02.pdf', width=8, height=10)
DoHeatmap(object = KCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


###################################################################################################
# Average expression - Keratinocytes : non-responder serial ####
Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

memory.limit()
memory.limit(size=60000)

levels(Round1)
KCs <- subset(Round1, idents = c( c("KC-S.Corneum",  "KC-S.Granulosum" ,  "KC-S.Spinosum"  ,  "KC-S.Basale"  )))

# Non-responder
colnames(KCs@meta.data)
Idents(object = KCs) <- ("number")
levels(KCs)
KCs <- subset(KCs, idents = c("Psoriasis07_preTx_LS_02",
                              "Psoriasis07_postTx_week12_LS" ,
                              "Psoriasis07_postTx_week48_LS" 
))

#Control vs Psoriasis07 serial
Idents(object = KCs) <- ("Serial")
levels(KCs)
KCs <- subset(KCs, idents = c(
  "Psoriasis_PreTx_LS_week0" ,
  "Psoriasis_PostTx_week12"   ,
  "Psoriasis_PostTx_LS_week48"  ))

Idents(object = KCs) <- ("Serial")
levels(KCs)
Idents(object = KCs) <- ("ClusterNames_0.8")

KCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(KCs), KCs$Serial, sep = "_")
Idents(KCs) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(KCs)

Idents(object = KCs) <- factor(Idents(object = KCs), levels = (c(
  "KC-S.Corneum_Psoriasis_PreTx_LS_week0" ,    
  "KC-S.Corneum_Psoriasis_PostTx_week12" ,
  "KC-S.Corneum_Psoriasis_PostTx_LS_week48" ,
  "KC-S.Granulosum_Psoriasis_PreTx_LS_week0" ,
  "KC-S.Granulosum_Psoriasis_PostTx_week12" ,
  "KC-S.Granulosum_Psoriasis_PostTx_LS_week48",
  "KC-S.Spinosum_Psoriasis_PreTx_LS_week0" ,   
  "KC-S.Spinosum_Psoriasis_PostTx_week12"   ,  
  "KC-S.Spinosum_Psoriasis_PostTx_LS_week48" ,
  "KC-S.Basale_Psoriasis_PreTx_LS_week0"  ,
  "KC-S.Basale_Psoriasis_PostTx_week12"  ,
  "KC-S.Basale_Psoriasis_PostTx_LS_week48"   
)))


KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)
KCs.cluster.averages <- AverageExpression(KCs, return.seurat = TRUE)

pdf('KCs.averages_Psoriasis07_nonresponder_Heatmp01.pdf', width=8, height=3.2)
DoHeatmap(object = KCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('KCs.averages_Psoriasis07_nonresponder_Heatmp02.pdf', width=8, height=5)
DoHeatmap(object = KCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

pdf('KCs.averages_Psoriasis07_nonresponder_Heatmp_shortened01.pdf', width=8, height=2)
DoHeatmap(object = KCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('KCs.averages_Psoriasis07_nonresponder_Heatmp_shortened02.pdf', width=8, height=5)
DoHeatmap(object = KCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()
###################################################################################################
# Average expression - Keratinocites : responder serial ####
Idents(object = Round1) <- ("ClusterNames_0.8")
KCs <- subset(Round1, idents = c( c("KC-S.Corneum",  "KC-S.Granulosum" ,  "KC-S.Spinosum"  ,  "KC-S.Basale"  )))

#Remove Non-responder
colnames(KCs@meta.data)
Idents(object = KCs) <- ("number")
levels(KCs)
KCs <- subset(KCs, idents = c(
  "Psoriasis08_preTx_LS"    ,
  "Psoriasis08_postTx_week12_LS",
  "Psoriasis08_postTx_week24_LS"
))

#Control vs Psoriasis07 serial
Idents(object = KCs) <- ("Serial")
levels(KCs)
KCs <- subset(KCs, idents = c(
  "Psoriasis_PreTx_LS_week0" ,
  "Psoriasis_PostTx_week12" ,
  "Psoriasis_PostTx_LS_week24" ))


Idents(object = KCs) <- ("ClusterNames_0.8")
KCs$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(KCs), KCs$Serial, sep = "_")
Idents(KCs) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(KCs)
Idents(KCs) <- factor(Idents(KCs), levels = c(
  "KC-S.Corneum_Psoriasis_PreTx_LS_week0"  ,
  "KC-S.Corneum_Psoriasis_PostTx_week12",      
  "KC-S.Corneum_Psoriasis_PostTx_LS_week24" ,  
  "KC-S.Granulosum_Psoriasis_PreTx_LS_week0",
  "KC-S.Granulosum_Psoriasis_PostTx_week12" ,
  "KC-S.Granulosum_Psoriasis_PostTx_LS_week24",
  "KC-S.Spinosum_Psoriasis_PreTx_LS_week0"  ,
  "KC-S.Spinosum_Psoriasis_PostTx_week12"  ,
  "KC-S.Spinosum_Psoriasis_PostTx_LS_week24" ,
  "KC-S.Basale_Psoriasis_PreTx_LS_week0"   ,   
  "KC-S.Basale_Psoriasis_PostTx_week12"  ,
  "KC-S.Basale_Psoriasis_PostTx_LS_week24" 
))

KCs[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = KCs)
Idents(object = KCs) <- ("ClusterNames_0.8_PsoriasisvsControl")


KCs.cluster.averages <- AverageExpression(KCs, return.seurat = TRUE)

pdf('KCs.averages_Psoriasis08_responder_Heatmp01.pdf', width=8, height=3.2)
DoHeatmap(object = KCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('KCs.averages_Psoriasis08_responder_Heatmp02.pdf', width=8, height=5)
DoHeatmap(object = KCs.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('KCs.averages_Psoriasis08_responder_shortened_Heatmp01.pdf', width=8, height=2)
DoHeatmap(object = KCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('KCs.averages_Psoriasis08_responder_shortened_Heatmp02.pdf', width=8, height=5)
DoHeatmap(object = KCs.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

#####################################################################
# Average expression - T-cells ####
Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"
levels(Round1)

T_cells <- subset(Round1, idents = c( "NK_cell"      ,    "CD161_T_cell"  ,   "CD8_T_cell"       ,"CD4_T_cell"     ,  "Treg"    ))

#Control vs Psoriasis_PreTx_LS_week0 vs Psoriasis_PostTx_week12
Idents(object = T_cells) <- ("Serial")
levels(T_cells)
T_cells <- subset(T_cells, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))

#Remove Non-responder
colnames(T_cells@meta.data)
Idents(object = T_cells) <- ("number")
levels(T_cells)
T_cells <- subset(T_cells, idents = c("Psoriasis07_postTx_week12_LS"  ), inv=TRUE)


Idents(object = T_cells) <- ("Serial")
levels(T_cells)

marker <- c("CD3D" , "KLRB1"  ,"IFNL1",  "TNFSF11",  "GNLY"  , "PRF1"  ,
            "CD8A"   , "CD8B"   ,      "CD8A"   , "CD8B"  ,
            "GZMK"   , "GZMH"   , 
            "IFNG"  ,  "IL26"  , "IL17A"  , "IL17F" ,
        #    "IL23R",
            #"HIF1A",
        #    "IL13","IL4", 
             "TIGIT"  ,  "FOXP3" ,  "IL2RA" 
        #"CTLA4"  
        #    "PDCD1"
                     #   "ITGAE",
        #    "ITGA1",
        #    "CXCR6",
        #     "CXCR3",
        #    "CD69",
        #    "CD40LG",
        #    "MME",
        #    "IL15"
                  )


marker_shortened <- c(
            "CD3D","KLRB1", "GZMK", "IL26"  , "IL17A"  , "IL17F" ,
   #         "IL23R",
            #    "IL13","IL4", 
            "FOXP3" ,  "IL2RA" 
            #"CTLA4"  
            #    "PDCD1"
            #   "ITGAE",
            #    "ITGA1",
            #    "CXCR6",
            #     "CXCR3",
            #    "CD69",
            #    "CD40LG",
            #    "MME",
            #    "IL15"
)
Idents(object = T_cells) <- ("ClusterNames_0.8")

T_cells$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(T_cells), T_cells$Serial, sep = "_")
Idents(T_cells) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(T_cells)
T_cells[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = T_cells)
levels(T_cells)

Idents(object = T_cells) <- factor(Idents(object = T_cells), levels = (c(
  "NK_cell_Control" ,"NK_cell_Psoriasis_PreTx_LS_week0" ,
                                                                          "NK_cell_Psoriasis_PostTx_week12"  ,     
  "CD161_T_cell_Control",
  "CD161_T_cell_Psoriasis_PreTx_LS_week0",
  "CD161_T_cell_Psoriasis_PostTx_week12"  ,
  "CD8_T_cell_Control" ,
  "CD8_T_cell_Psoriasis_PreTx_LS_week0",
  "CD8_T_cell_Psoriasis_PostTx_week12" ,
  "CD4_T_cell_Control",
  "CD4_T_cell_Psoriasis_PreTx_LS_week0" , 
  "CD4_T_cell_Psoriasis_PostTx_week12" ,
  "Treg_Control"   ,                      
                                                                         "Treg_Psoriasis_PreTx_LS_week0" ,
                                                                         "Treg_Psoriasis_PostTx_week12"     
  )))


T_cells_clusters.averages <- AverageExpression(T_cells, return.seurat = TRUE)


pdf('T_cells.averages_Round1_Heatmp01.pdf', width=7, height=4)
DoHeatmap(object = T_cells_clusters.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

pdf('T_cells.averages_Round1_Heatmp02.pdf', width=7, height=8)
DoHeatmap(object = T_cells_clusters.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Round1_shoretened_Heatmp01.pdf', width=7, height=2.5)
DoHeatmap(object = T_cells_clusters.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Round1_shoretened_Heatmp02.pdf', width=7, height=10)
DoHeatmap(object = T_cells_clusters.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()

###################################################################################################
# Average expression - T_cells : non-responder serial ####
Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

levels(Round1)
T_cells <- subset(Round1, idents = c( "NK_cell"      ,    "CD161_T_cell"  ,   "CD8_T_cell"       ,"CD4_T_cell"     ,  "Treg"    ))

# Non-responder
colnames(T_cells@meta.data)
Idents(object = T_cells) <- ("number")
levels(T_cells)
T_cells <- subset(T_cells, idents = c("Psoriasis07_preTx_LS_02",
                                      "Psoriasis07_postTx_week12_LS" ,
                                      "Psoriasis07_postTx_week48_LS" 
))

#Control vs Psoriasis07 serial
Idents(object = T_cells) <- ("Serial")
levels(T_cells)
T_cells <- subset(T_cells, idents = c(
  "Psoriasis_PreTx_LS_week0" ,
  "Psoriasis_PostTx_week12"   ,
  "Psoriasis_PostTx_LS_week48"  ))

Idents(object = T_cells) <- ("Serial")
levels(T_cells)
Idents(object = T_cells) <- ("ClusterNames_0.8")

T_cells$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(T_cells), T_cells$Serial, sep = "_")
Idents(T_cells) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(T_cells)

Idents(object = T_cells) <- factor(Idents(object = T_cells), levels = (c(
  "NK_cell_Psoriasis_PreTx_LS_week0" ,
  "NK_cell_Psoriasis_PostTx_week12"   ,
  "NK_cell_Psoriasis_PostTx_LS_week48"  ,
  "CD161_T_cell_Psoriasis_PreTx_LS_week0"    ,
  "CD161_T_cell_Psoriasis_PostTx_week12",
  "CD161_T_cell_Psoriasis_PostTx_LS_week48",
  "CD8_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD8_T_cell_Psoriasis_PostTx_week12"  ,
  "CD8_T_cell_Psoriasis_PostTx_LS_week48" ,
  "CD4_T_cell_Psoriasis_PreTx_LS_week0",
  "CD4_T_cell_Psoriasis_PostTx_week12" ,    
  "CD4_T_cell_Psoriasis_PostTx_LS_week48"  ,
  "Treg_Psoriasis_PreTx_LS_week0"          ,
  "Treg_Psoriasis_PostTx_week12"    ,
  "Treg_Psoriasis_PostTx_LS_week48" 
)))


T_cells[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = T_cells)
T_cells.cluster.averages <- AverageExpression(T_cells, return.seurat = TRUE)

pdf('T_cells.averages_Psoriasis07_nonresponder_Heatmap01.pdf', width=7, height=4)
DoHeatmap(object = T_cells.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis07_nonresponder_Heatmap02.pdf', width=7, height=8)
DoHeatmap(object = T_cells.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis07_nonresponder_shortened_Heatmap01.pdf', width=7, height=2.5)
DoHeatmap(object = T_cells.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis07_nonresponder_shortened_Heatmap02.pdf', width=7, height=10)
DoHeatmap(object = T_cells.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()
###################################################################################################
# Average expression - T_cells : responder serial ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Idents(object = Round1) <- ("ClusterNames_0.8")
T_cells <- subset(Round1, idents = c( "NK_cell"      ,    "CD161_T_cell"  ,   "CD8_T_cell"       ,"CD4_T_cell"     ,  "Treg"    ))

#Remove Non-responder
colnames(T_cells@meta.data)
Idents(object = T_cells) <- ("number")
levels(T_cells)
T_cells <- subset(T_cells, idents = c(
  "Psoriasis08_preTx_LS"    ,
  "Psoriasis08_postTx_week12_LS",
  "Psoriasis08_postTx_week24_LS"
))

#Control vs Psoriasis07 serial
Idents(object = T_cells) <- ("Serial")
levels(T_cells)
T_cells <- subset(T_cells, idents = c(
  "Psoriasis_PreTx_LS_week0" ,
  "Psoriasis_PostTx_week12" ,
  "Psoriasis_PostTx_LS_week24" ))


Idents(object = T_cells) <- ("ClusterNames_0.8")
T_cells$ClusterNames_0.8_PsoriasisvsControl <- paste(Idents(T_cells), T_cells$Serial, sep = "_")
Idents(T_cells) <- "ClusterNames_0.8_PsoriasisvsControl"
levels(T_cells)
Idents(T_cells) <- factor(Idents(T_cells), levels = c(
  "NK_cell_Psoriasis_PreTx_LS_week0"  ,
  "NK_cell_Psoriasis_PostTx_week12"    ,    
  "NK_cell_Psoriasis_PostTx_LS_week24"   ,
  "CD161_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD161_T_cell_Psoriasis_PostTx_week12" ,
  "CD161_T_cell_Psoriasis_PostTx_LS_week24",
  "CD8_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD8_T_cell_Psoriasis_PostTx_week12",     
  "CD8_T_cell_Psoriasis_PostTx_LS_week24" ,
  "CD4_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD4_T_cell_Psoriasis_PostTx_week12" ,
  "CD4_T_cell_Psoriasis_PostTx_LS_week24",  
  "Treg_Psoriasis_PreTx_LS_week0"    ,      
  "Treg_Psoriasis_PostTx_week12"  ,
  "Treg_Psoriasis_PostTx_LS_week24"    
))

T_cells[["ClusterNames_0.8_PsoriasisvsControl"]] <- Idents(object = T_cells)
Idents(object = T_cells) <- ("ClusterNames_0.8_PsoriasisvsControl")


T_cells.cluster.averages <- AverageExpression(T_cells, return.seurat = TRUE)

pdf('T_cells.averages_Psoriasis08_responder_Heatmap01.pdf', width=7, height=4)
DoHeatmap(object = T_cells.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis08_responder_Heatmap02.pdf', width=7, height=8)
DoHeatmap(object = T_cells.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis08_responder_shortened_Heatmap01.pdf', width=7, height=2.5)
DoHeatmap(object = T_cells.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis08_responder_shortened_Heatmap02.pdf', width=7, height=10)
DoHeatmap(object = T_cells.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
  # scale_fill_viridis(option="magma")
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 10, name = "RdBu"))
  scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()
#####################################################################################################
## Save
save(Round1, file = "~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round1_integrated_analyzed_6.1.2021.Rda")

#################################################################################################
