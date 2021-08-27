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

dir.create("./coexpression_T_cells")
setwd("./coexpression_T_cells")
getwd()
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/coexpression_T_cells")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/coexpression_T_cells")

list.files()
rm(list= ls()[!(ls() %in% c('Round1'))])

gene_subset  <- read_sheet("https://docs.google.com/spreadsheets/d/19nhrzGudrn8ihl4RASuyGlBe7JMHdWnnWKoZ0i2rZ4s/edit#gid=0")
gene_subset <- as.matrix(unique(gene_subset[1]))

#load("/Users/jkim05/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_01.31.2021.Rda")
#load("/Users/jkim05/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_01.31.2021.Rda")
#load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_02.09.2020.Rda")
#load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1/combined_data/Round1_integrated_analyzed_02.15.2020.Rda")
######################################################################################################
## T_cell IL17A vs IL17F co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

T_cell <- subset(Round2, idents = c("NK_cell"   ,
                                    "CD161_T_cell"  ,
                                    "CD8_T_cell"   ,
                                    "CD4_T_cell"     ,
                                    "Treg"   ))

IL17A_T_cell <- subset(T_cell, subset = IL17A > 1)
IL17A_T_cell <- subset(IL17A_T_cell, subset = IL17F < 1)
IL17A_T_cell_list <- as.data.frame(table(IL17A_T_cell@meta.data$ClusterNames_0.8_Serial))
IL17A_T_cell_list[is.na(IL17A_T_cell_list)] <- 0
colnames(IL17A_T_cell_list) <- c("Cluster","IL17A")

IL17A_IL17F_T_cell <- subset(T_cell, subset = IL17F > 1)
IL17A_IL17F_T_cell <- subset(IL17A_IL17F_T_cell, subset = IL17A > 1)
IL17A_IL17F_T_cell_list <- as.data.frame(table(IL17A_IL17F_T_cell@meta.data$ClusterNames_0.8_Serial))
IL17A_IL17F_T_cell_list[is.na(IL17A_IL17F_T_cell_list)] <- 0
colnames(IL17A_IL17F_T_cell_list) <- c("Cluster","IL17A&IL17F")


IL17F_T_cell <- subset(T_cell, subset = IL17F > 1)
IL17F_T_cell <- subset(IL17F_T_cell, subset = IL17A < 1)
IL17F_T_cell_list <- as.data.frame(table(IL17F_T_cell@meta.data$ClusterNames_0.8_Serial))
IL17F_T_cell_list[is.na(IL17F_T_cell_list)] <- 0
colnames(IL17F_T_cell_list) <- c("Cluster","IL17F")

list <- merge(IL17A_T_cell_list, IL17A_IL17F_T_cell_list, by ="Cluster", all=TRUE)
list <- merge(list, IL17F_T_cell_list, by  ="Cluster",all=TRUE)
list[is.na(list)] <- 0

total_cell_count <- as.data.frame(table(T_cell@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster", all.y = TRUE)
list[is.na(list)] <- 0

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_T_cell_IL17A_IL17F.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL17F","IL17A&IL17F", "IL17A" )))
levels(list_for_ggplot$Cluster)
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c( 
  "NK_cell_Control"   ,
  "NK_cell_Psoriasis_PreTx_LS_week0",
  "NK_cell_Psoriasis_PostTx_week12"  ,    
  "CD161_T_cell_Control"          ,   
  "CD161_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD161_T_cell_Psoriasis_PostTx_week12" ,
  "CD8_T_cell_Control"       ,
  "CD8_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD8_T_cell_Psoriasis_PostTx_week12"   ,
  "CD4_T_cell_Control"        ,           
  "CD4_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD4_T_cell_Psoriasis_PostTx_week12" ,
  "Treg_Control"  ,                       
  "Treg_Psoriasis_PreTx_LS_week0"     ,   
  "Treg_Psoriasis_PostTx_week12"      
)))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of IL17A and IL17F expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_T_cell_IL17A_IL17F_count.pdf", width = 7, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL17F","IL17A&IL17F", "IL17A" )))
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c( 
  "NK_cell_Control"   ,
  "NK_cell_Psoriasis_PreTx_LS_week0",
  "NK_cell_Psoriasis_PostTx_week12"  ,    
  "CD161_T_cell_Control"          ,   
  "CD161_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD161_T_cell_Psoriasis_PostTx_week12" ,
  "CD8_T_cell_Control"       ,
  "CD8_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD8_T_cell_Psoriasis_PostTx_week12"   ,
  "CD4_T_cell_Control"        ,           
  "CD4_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD4_T_cell_Psoriasis_PostTx_week12" ,
  "Treg_Control"  ,                       
  "Treg_Psoriasis_PreTx_LS_week0"     ,   
  "Treg_Psoriasis_PostTx_week12"      
)))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of IL17A and IL17F expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_IL17A_IL17F_proportion.pdf", width = 7, height = 6)

######################################################################################################
## T_cell CD8B vs KLRB1 co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

T_cell <- subset(Round2, idents = c("NK_cell"   ,
                                    "CD161_T_cell"  ,
                                    "CD8_T_cell"   ,
                                    "CD4_T_cell"     ,
                                    "Treg"   ))

CD8B_T_cell <- subset(T_cell, subset = CD8B > 1)
CD8B_T_cell <- subset(CD8B_T_cell, subset = KLRB1 < 1)
CD8B_T_cell_list <- as.data.frame(table(CD8B_T_cell@meta.data$ClusterNames_0.8_Serial))
CD8B_T_cell_list[is.na(CD8B_T_cell_list)] <- 0
colnames(CD8B_T_cell_list) <- c("Cluster","CD8B")

CD8B_KLRB1_T_cell <- subset(T_cell, subset = KLRB1 > 1)
CD8B_KLRB1_T_cell <- subset(CD8B_KLRB1_T_cell, subset = CD8B > 1)
CD8B_KLRB1_T_cell_list <- as.data.frame(table(CD8B_KLRB1_T_cell@meta.data$ClusterNames_0.8_Serial))
CD8B_KLRB1_T_cell_list[is.na(CD8B_KLRB1_T_cell_list)] <- 0
colnames(CD8B_KLRB1_T_cell_list) <- c("Cluster","CD8B&KLRB1")


KLRB1_T_cell <- subset(T_cell, subset = KLRB1 > 1)
KLRB1_T_cell <- subset(KLRB1_T_cell, subset = CD8B < 1)
KLRB1_T_cell_list <- as.data.frame(table(KLRB1_T_cell@meta.data$ClusterNames_0.8_Serial))
KLRB1_T_cell_list[is.na(KLRB1_T_cell_list)] <- 0
colnames(KLRB1_T_cell_list) <- c("Cluster","KLRB1")

list <- merge(CD8B_T_cell_list, CD8B_KLRB1_T_cell_list, by ="Cluster", all=TRUE)
list <- merge(list, KLRB1_T_cell_list, by  ="Cluster",all=TRUE)
list[is.na(list)] <- 0

total_cell_count <- as.data.frame(table(T_cell@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster", all.y = TRUE)
list[is.na(list)] <- 0

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_T_cell_CD8B_KLRB1.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("KLRB1","CD8B&KLRB1", "CD8B" )))
levels(list_for_ggplot$Cluster)
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c( 
  "NK_cell_Control"   ,
  "NK_cell_Psoriasis_PreTx_LS_week0",
  "NK_cell_Psoriasis_PostTx_week12"  ,    
  "CD161_T_cell_Control"          ,   
  "CD161_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD161_T_cell_Psoriasis_PostTx_week12" ,
  "CD8_T_cell_Control"       ,
  "CD8_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD8_T_cell_Psoriasis_PostTx_week12"   ,
  "CD4_T_cell_Control"        ,           
  "CD4_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD4_T_cell_Psoriasis_PostTx_week12" ,
  "Treg_Control"  ,                       
  "Treg_Psoriasis_PreTx_LS_week0"     ,   
  "Treg_Psoriasis_PostTx_week12"      
)))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of CD8B and KLRB1 expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_T_cell_CD8B_KLRB1_count.pdf", width = 7, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("KLRB1","CD8B&KLRB1", "CD8B" )))
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c( 
  "NK_cell_Control"   ,
  "NK_cell_Psoriasis_PreTx_LS_week0",
  "NK_cell_Psoriasis_PostTx_week12"  ,    
  "CD161_T_cell_Control"          ,   
  "CD161_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD161_T_cell_Psoriasis_PostTx_week12" ,
  "CD8_T_cell_Control"       ,
  "CD8_T_cell_Psoriasis_PreTx_LS_week0" ,
  "CD8_T_cell_Psoriasis_PostTx_week12"   ,
  "CD4_T_cell_Control"        ,           
  "CD4_T_cell_Psoriasis_PreTx_LS_week0"  ,
  "CD4_T_cell_Psoriasis_PostTx_week12" ,
  "Treg_Control"  ,                       
  "Treg_Psoriasis_PreTx_LS_week0"     ,   
  "Treg_Psoriasis_PostTx_week12"      
)))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of CD8B and KLRB1 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_CD8B_KLRB1_proportion.pdf", width = 7, height = 6)


######################################################################################################
## IL17A / IL17F  T_cell CD161 / CD8A co-expression ####

#IL17A / IL17F  T_cell
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
#Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Round2 <- subset(Round2, idents = c( "Psoriasis_PreTx_LS_week0"  ))
Idents(object = Round2) <- ("number")
#Control 10
#preTx_LS 11
#postTx_week12_LS 4

Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

#T_cell <- subset(Round2, idents = c("CD161_T_cell" ))

T_cell <- subset(Round2, idents = c(
  "CD161_T_cell"  ,
  "CD8_T_cell"   ,
  "CD4_T_cell"     ,
  "Treg"   ))


IL17A_T_cell <- subset(T_cell, subset = IL17A > 1)
IL17A_T_cell <- subset(IL17A_T_cell, subset = IL17F < 1)
levels(IL17A_T_cell)
new.cluster.ids <- c(  "IL17A_T_cell")
names(new.cluster.ids) <- levels(IL17A_T_cell)
IL17A_T_cell <- RenameIdents(IL17A_T_cell, new.cluster.ids)
levels(Idents(IL17A_T_cell))

IL17A_IL17F_T_cell <- subset(T_cell, subset = IL17F > 1)
IL17A_IL17F_T_cell <- subset(IL17A_IL17F_T_cell, subset = IL17A > 1)
levels(IL17A_IL17F_T_cell)
new.cluster.ids <- c(   "IL17A_IL17F_T_cell")
names(new.cluster.ids) <- levels(IL17A_IL17F_T_cell)
IL17A_IL17F_T_cell <- RenameIdents(IL17A_IL17F_T_cell, new.cluster.ids)
levels(Idents(IL17A_IL17F_T_cell))

IL17F_T_cell <- subset(T_cell, subset = IL17F > 1)
IL17F_T_cell <- subset(IL17F_T_cell, subset = IL17A < 1)
levels(IL17F_T_cell)
new.cluster.ids <- c(    "IL17F_T_cell")
names(new.cluster.ids) <- levels(IL17F_T_cell)
IL17F_T_cell <- RenameIdents(IL17F_T_cell, new.cluster.ids)
levels(Idents(IL17F_T_cell))

IL17A_IL17F_T_cell <- merge(x= IL17A_T_cell, y = c(IL17A_IL17F_T_cell, IL17F_T_cell ), project = "IL17A_IL17F_T_cell")

levels(Idents(IL17A_IL17F_T_cell))
new.cluster.ids <- c("IL17A_IL17F_T_cell",
                     "IL17A_IL17F_T_cell",
                     "IL17A_IL17F_T_cell")
names(new.cluster.ids) <- levels(IL17A_IL17F_T_cell)
IL17A_IL17F_T_cell <- RenameIdents(IL17A_IL17F_T_cell, new.cluster.ids)
levels(Idents(IL17A_IL17F_T_cell))
IL17A_IL17F_T_cell[["ClusterNames_0.8"]] <- Idents(object = IL17A_IL17F_T_cell)

# KLRB1 / CD8A co-expression

CD161_IL17A_IL17F_T_cell <- subset(IL17A_IL17F_T_cell, subset = KLRB1 > 1)
CD161_IL17A_IL17F_T_cell <- subset(CD161_IL17A_IL17F_T_cell, subset = CD8A < 1)
CD161_IL17A_IL17F_T_cell_list <- as.data.frame(table(CD161_IL17A_IL17F_T_cell@meta.data$ClusterNames_0.8))
CD161_IL17A_IL17F_T_cell_list[is.na(CD161_IL17A_IL17F_T_cell_list)] <- 0
colnames(CD161_IL17A_IL17F_T_cell_list) <- c("Cluster","CD161")

CD161_CD8A_IL17A_IL17F_T_cell <- subset(IL17A_IL17F_T_cell, subset = CD8A > 1)
CD161_CD8A_IL17A_IL17F_T_cell <- subset(CD161_CD8A_IL17A_IL17F_T_cell, subset = KLRB1 > 1)
CD161_CD8A_IL17A_IL17F_T_cell_list <- as.data.frame(table(CD161_CD8A_IL17A_IL17F_T_cell@meta.data$ClusterNames_0.8))
CD161_CD8A_IL17A_IL17F_T_cell_list[is.na(CD161_CD8A_IL17A_IL17F_T_cell_list)] <- 0
colnames(CD161_CD8A_IL17A_IL17F_T_cell_list) <- c("Cluster","CD161&CD8A")


CD8A_IL17A_IL17F_T_cell <- subset(IL17A_IL17F_T_cell, subset = CD8A > 1)
CD8A_IL17A_IL17F_T_cell <- subset(CD8A_IL17A_IL17F_T_cell, subset = KLRB1 < 1)
CD8A_IL17A_IL17F_T_cell_list <- as.data.frame(table(CD8A_IL17A_IL17F_T_cell@meta.data$ClusterNames_0.8))
CD8A_IL17A_IL17F_T_cell_list[is.na(CD8A_IL17A_IL17F_T_cell_list)] <- 0
colnames(CD8A_IL17A_IL17F_T_cell_list) <- c("Cluster","CD8A")

list <- merge(CD161_IL17A_IL17F_T_cell_list, CD161_CD8A_IL17A_IL17F_T_cell_list, by ="Cluster", all=TRUE)
list <- merge(list, CD8A_IL17A_IL17F_T_cell_list, by  ="Cluster",all=TRUE)
list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(IL17A_IL17F_T_cell@meta.data$ClusterNames_0.8))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster", all.y = TRUE)
list[is.na(list)] <- 0

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
#list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
#list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_IL17A_IL17F_T_cell_CD161_CD8A.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CD8A","CD161&CD8A", "CD161" )))
levels(list_for_ggplot$Cluster)
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c("Mature_DC_Control" ,
                                                                      "Mature_DC_Psoriasis_PreTx_LS_week0" ,
                                                                      "Mature_DC_Psoriasis_PostTx_week12",
                                                                      "T_cell_Control"  ,
                                                                      "T_cell_Psoriasis_PreTx_LS_week0",
                                                                      "T_cell_Psoriasis_PostTx_week12")))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of CD161 and CD8A expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_IL17A_IL17F_T_cell_CD161_CD8A_count.pdf", width = 3, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CD8A","CD161&CD8A", "CD161" )))
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c("Mature_DC_Control" ,
                                                                      "Mature_DC_Psoriasis_PreTx_LS_week0" ,
                                                                      "Mature_DC_Psoriasis_PostTx_week12",
                                                                      "T_cell_Control"  ,
                                                                      "T_cell_Psoriasis_PreTx_LS_week0",
                                                                      "T_cell_Psoriasis_PostTx_week12")))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of CD161 and CD8A expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_IL17A_IL17F_T_cell_CD161_CD8A_proportion.pdf", width = 3, height = 6)


######################################################################################################
## IL10  Semimature_DC THBD / CLEC4A co-expression ####

#IL10  Semimature_DC
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
#Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Round2 <- subset(Round2, idents = c( "Psoriasis_PreTx_LS_week0"  ))
Idents(object = Round2) <- ("number")
#Control 10
#preTx_LS 11
#postTx_week12_LS 4

Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Semimature_DC <- subset(Round2, idents = c("Semimature_DC" ))

IL10_Semimature_DC <- subset(Semimature_DC, subset = IL10 > 1)
levels(IL10_Semimature_DC)
new.cluster.ids <- c(  "IL10_Semimature_DC")
names(new.cluster.ids) <- levels(IL10_Semimature_DC)
IL10_Semimature_DC <- RenameIdents(IL10_Semimature_DC, new.cluster.ids)
levels(Idents(IL10_Semimature_DC))
IL10_Semimature_DC[["ClusterNames_0.8"]] <- Idents(object = IL10_Semimature_DC)


# THBD / CLEC4A co-expression
THBD_IL10_Semimature_DC <- subset(IL10_Semimature_DC, subset = THBD > 1)
THBD_IL10_Semimature_DC <- subset(THBD_IL10_Semimature_DC, subset = CLEC4A < 1)
THBD_IL10_Semimature_DC_list <- as.data.frame(table(THBD_IL10_Semimature_DC@meta.data$ClusterNames_0.8))
THBD_IL10_Semimature_DC_list[is.na(THBD_IL10_Semimature_DC_list)] <- 0
colnames(THBD_IL10_Semimature_DC_list) <- c("Cluster","THBD")

THBD_CLEC4A_IL10_Semimature_DC <- subset(IL10_Semimature_DC, subset = CLEC4A > 1)
THBD_CLEC4A_IL10_Semimature_DC <- subset(THBD_CLEC4A_IL10_Semimature_DC, subset = THBD > 1)
THBD_CLEC4A_IL10_Semimature_DC_list <- as.data.frame(table(THBD_CLEC4A_IL10_Semimature_DC@meta.data$ClusterNames_0.8))
THBD_CLEC4A_IL10_Semimature_DC_list[is.na(THBD_CLEC4A_IL10_Semimature_DC_list)] <- 0
colnames(THBD_CLEC4A_IL10_Semimature_DC_list) <- c("Cluster","THBD&CLEC4A")


CLEC4A_IL10_Semimature_DC <- subset(IL10_Semimature_DC, subset = CLEC4A > 1)
CLEC4A_IL10_Semimature_DC <- subset(CLEC4A_IL10_Semimature_DC, subset = THBD < 1)
CLEC4A_IL10_Semimature_DC_list <- as.data.frame(table(CLEC4A_IL10_Semimature_DC@meta.data$ClusterNames_0.8))
CLEC4A_IL10_Semimature_DC_list[is.na(CLEC4A_IL10_Semimature_DC_list)] <- 0
colnames(CLEC4A_IL10_Semimature_DC_list) <- c("Cluster","CLEC4A")

list <- merge(THBD_IL10_Semimature_DC_list, THBD_CLEC4A_IL10_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, CLEC4A_IL10_Semimature_DC_list, by  ="Cluster",all=TRUE)
list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(IL10_Semimature_DC@meta.data$ClusterNames_0.8))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster", all.y = TRUE)
list[is.na(list)] <- 0

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


#list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
#list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
#list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_IL10_Semimature_DC_THBD_CLEC4A.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CLEC4A","THBD&CLEC4A", "THBD" )))
levels(list_for_ggplot$Cluster)
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c("Mature_DC_Control" ,
                                                                      "Mature_DC_Psoriasis_PreTx_LS_week0" ,
                                                                      "Mature_DC_Psoriasis_PostTx_week12",
                                                                      "Semimature_DC_Control"  ,
                                                                      "Semimature_DC_Psoriasis_PreTx_LS_week0",
                                                                      "Semimature_DC_Psoriasis_PostTx_week12")))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of THBD and CLEC4A expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_IL10_Semimature_DC_THBD_CLEC4A_count.pdf", width = 3, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("CLEC4A","THBD&CLEC4A", "THBD" )))
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c("Mature_DC_Control" ,
                                                                      "Mature_DC_Psoriasis_PreTx_LS_week0" ,
                                                                      "Mature_DC_Psoriasis_PostTx_week12",
                                                                      "Semimature_DC_Control"  ,
                                                                      "Semimature_DC_Psoriasis_PreTx_LS_week0",
                                                                      "Semimature_DC_Psoriasis_PostTx_week12")))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of THBD and CLEC4A expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_IL10_Semimature_DC_THBD_CLEC4A_proportion.pdf", width = 3, height = 6)

######################################################################################################