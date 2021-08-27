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

dir.create("./coexpression_DCs")
setwd("./coexpression_DCs")
getwd()
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/coexpression_DCs")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/coexpression_DCs")

list.files()
rm(list= ls()[!(ls() %in% c('Round1'))])


gene_subset  <- read_sheet("https://docs.google.com/spreadsheets/d/19nhrzGudrn8ihl4RASuyGlBe7JMHdWnnWKoZ0i2rZ4s/edit#gid=0")
gene_subset <- as.matrix(unique(gene_subset[1]))

#load("/Users/jkim05/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_01.31.2021.Rda")
#load("/Users/jkim05/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_01.31.2021.Rda")
#load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/combined_data/Round0_integrated_analyzed_02.09.2020.Rda")
#load("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1/combined_data/Round1_integrated_analyzed_02.15.2020.Rda")

######################################################################################################3
# ClusterNames_0.8_Serial ####
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


######################################################################################################
## Mature_DC and Semimature_DC IL23A vs FCGR1A co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

DefaultAssay(Mature_DC_Semimature_DC) <- "RNA"

IL23A_Mature_DC_Semimature_DC <- Mature_DC_Semimature_DC

IL23A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL23A > 1)
IL23A_Mature_DC_Semimature_DC <- subset(IL23A_Mature_DC_Semimature_DC, subset = FCGR1A < 1)
IL23A_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL23A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL23A_Mature_DC_Semimature_DC_list[is.na(IL23A_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL23A_Mature_DC_Semimature_DC_list) <- c("Cluster","IL23A")

IL23A_FCGR1A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = FCGR1A > 1)
IL23A_FCGR1A_Mature_DC_Semimature_DC <- subset(IL23A_FCGR1A_Mature_DC_Semimature_DC, subset = IL23A > 1)
IL23A_FCGR1A_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL23A_FCGR1A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL23A_FCGR1A_Mature_DC_Semimature_DC_list[is.na(IL23A_FCGR1A_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL23A_FCGR1A_Mature_DC_Semimature_DC_list) <- c("Cluster","IL23A&FCGR1A")


FCGR1A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = FCGR1A > 1)
FCGR1A_Mature_DC_Semimature_DC <- subset(FCGR1A_Mature_DC_Semimature_DC, subset = IL23A < 1)
FCGR1A_Mature_DC_Semimature_DC_list <- as.data.frame(table(FCGR1A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
FCGR1A_Mature_DC_Semimature_DC_list[is.na(FCGR1A_Mature_DC_Semimature_DC_list)] <- 0
colnames(FCGR1A_Mature_DC_Semimature_DC_list) <- c("Cluster","FCGR1A")

list <- merge(IL23A_Mature_DC_Semimature_DC_list, IL23A_FCGR1A_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, FCGR1A_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_IL23A_FCGR1A.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FCGR1A","IL23A&FCGR1A", "IL23A" )))
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
  ggtitle("Number of IL23A and FCGR1A expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_IL23A_FCGR1A_count.pdf", width = 4, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FCGR1A","IL23A&FCGR1A", "IL23A" )))
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
  ggtitle("Proportion of IL23A and FCGR1A expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_IL23A_FCGR1A_proportion.pdf", width = 4, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC IL23A vs IL10 co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

DefaultAssay(Mature_DC_Semimature_DC) <- "RNA"

IL23A_Mature_DC_Semimature_DC <- Mature_DC_Semimature_DC

IL23A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL23A > 1)
IL23A_Mature_DC_Semimature_DC <- subset(IL23A_Mature_DC_Semimature_DC, subset = IL10 < 1)
IL23A_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL23A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL23A_Mature_DC_Semimature_DC_list[is.na(IL23A_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL23A_Mature_DC_Semimature_DC_list) <- c("Cluster","IL23A")

IL23A_IL10_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL10 > 1)
IL23A_IL10_Mature_DC_Semimature_DC <- subset(IL23A_IL10_Mature_DC_Semimature_DC, subset = IL23A > 1)
IL23A_IL10_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL23A_IL10_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL23A_IL10_Mature_DC_Semimature_DC_list[is.na(IL23A_IL10_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL23A_IL10_Mature_DC_Semimature_DC_list) <- c("Cluster","IL23A&IL10")


IL10_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL10 > 1)
IL10_Mature_DC_Semimature_DC <- subset(IL10_Mature_DC_Semimature_DC, subset = IL23A < 1)
IL10_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL10_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL10_Mature_DC_Semimature_DC_list[is.na(IL10_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL10_Mature_DC_Semimature_DC_list) <- c("Cluster","IL10")

list <- merge(IL23A_Mature_DC_Semimature_DC_list, IL23A_IL10_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, IL10_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_IL23A_IL10.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL10","IL23A&IL10", "IL23A" )))
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
  ggtitle("Number of IL23A and IL10 expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_IL23A_IL10_count.pdf", width = 4, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL10","IL23A&IL10", "IL23A" )))
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
  ggtitle("Proportion of IL23A and IL10 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_IL23A_IL10_proportion.pdf", width = 4, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC IL23A vs ITGAX co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

IL23A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL23A > 1)
IL23A_Mature_DC_Semimature_DC <- subset(IL23A_Mature_DC_Semimature_DC, subset = ITGAX < 1)
IL23A_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL23A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL23A_Mature_DC_Semimature_DC_list[is.na(IL23A_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL23A_Mature_DC_Semimature_DC_list) <- c("Cluster","IL23A")

IL23A_ITGAX_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = ITGAX > 1)
IL23A_ITGAX_Mature_DC_Semimature_DC <- subset(IL23A_ITGAX_Mature_DC_Semimature_DC, subset = IL23A > 1)
IL23A_ITGAX_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL23A_ITGAX_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL23A_ITGAX_Mature_DC_Semimature_DC_list[is.na(IL23A_ITGAX_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL23A_ITGAX_Mature_DC_Semimature_DC_list) <- c("Cluster","IL23A&ITGAX")


ITGAX_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = ITGAX > 1)
ITGAX_Mature_DC_Semimature_DC <- subset(ITGAX_Mature_DC_Semimature_DC, subset = IL23A < 1)
ITGAX_Mature_DC_Semimature_DC_list <- as.data.frame(table(ITGAX_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
ITGAX_Mature_DC_Semimature_DC_list[is.na(ITGAX_Mature_DC_Semimature_DC_list)] <- 0
colnames(ITGAX_Mature_DC_Semimature_DC_list) <- c("Cluster","ITGAX")

list <- merge(IL23A_Mature_DC_Semimature_DC_list, IL23A_ITGAX_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, ITGAX_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100


list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_IL23A_ITGAX.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("ITGAX","IL23A&ITGAX", "IL23A" )))
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
  ggtitle("Number of IL23A and ITGAX expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_IL23A_ITGAX_count.pdf", width = 5, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("ITGAX","IL23A&ITGAX", "IL23A" )))
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
  ggtitle("Proportion of IL23A and ITGAX expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_IL23A_ITGAX_proportion.pdf", width = 5, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC ITGAX vs THBD co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

ITGAX_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = ITGAX > 1)
ITGAX_Mature_DC_Semimature_DC <- subset(ITGAX_Mature_DC_Semimature_DC, subset = THBD < 1)
ITGAX_Mature_DC_Semimature_DC_list <- as.data.frame(table(ITGAX_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
ITGAX_Mature_DC_Semimature_DC_list[is.na(ITGAX_Mature_DC_Semimature_DC_list)] <- 0
colnames(ITGAX_Mature_DC_Semimature_DC_list) <- c("Cluster","ITGAX")

ITGAX_THBD_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = THBD > 1)
ITGAX_THBD_Mature_DC_Semimature_DC <- subset(ITGAX_THBD_Mature_DC_Semimature_DC, subset = ITGAX > 1)
ITGAX_THBD_Mature_DC_Semimature_DC_list <- as.data.frame(table(ITGAX_THBD_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
ITGAX_THBD_Mature_DC_Semimature_DC_list[is.na(ITGAX_THBD_Mature_DC_Semimature_DC_list)] <- 0
colnames(ITGAX_THBD_Mature_DC_Semimature_DC_list) <- c("Cluster","ITGAX&THBD")


THBD_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = THBD > 1)
THBD_Mature_DC_Semimature_DC <- subset(THBD_Mature_DC_Semimature_DC, subset = ITGAX < 1)
THBD_Mature_DC_Semimature_DC_list <- as.data.frame(table(THBD_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
THBD_Mature_DC_Semimature_DC_list[is.na(THBD_Mature_DC_Semimature_DC_list)] <- 0
colnames(THBD_Mature_DC_Semimature_DC_list) <- c("Cluster","THBD")

list <- merge(ITGAX_Mature_DC_Semimature_DC_list, ITGAX_THBD_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, THBD_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_ITGAX_THBD.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("THBD","ITGAX&THBD", "ITGAX" )))
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
  ggtitle("Number of ITGAX and THBD expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_ITGAX_THBD_count.pdf", width = 5, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("THBD","ITGAX&THBD", "ITGAX" )))
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
  ggtitle("Proportion of ITGAX and THBD expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_ITGAX_THBD_proportion.pdf", width = 5, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC THBD vs CLEC4A co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("number")
#Control 10
#preTx_LS 11
#postTx_week12_LS 4

Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

THBD_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = THBD > 1)
THBD_Mature_DC_Semimature_DC <- subset(THBD_Mature_DC_Semimature_DC, subset = CLEC4A < 1)
THBD_Mature_DC_Semimature_DC_list <- as.data.frame(table(THBD_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
THBD_Mature_DC_Semimature_DC_list[is.na(THBD_Mature_DC_Semimature_DC_list)] <- 0
colnames(THBD_Mature_DC_Semimature_DC_list) <- c("Cluster","THBD")

THBD_CLEC4A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = CLEC4A > 1)
THBD_CLEC4A_Mature_DC_Semimature_DC <- subset(THBD_CLEC4A_Mature_DC_Semimature_DC, subset = THBD > 1)
THBD_CLEC4A_Mature_DC_Semimature_DC_list <- as.data.frame(table(THBD_CLEC4A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
THBD_CLEC4A_Mature_DC_Semimature_DC_list[is.na(THBD_CLEC4A_Mature_DC_Semimature_DC_list)] <- 0
colnames(THBD_CLEC4A_Mature_DC_Semimature_DC_list) <- c("Cluster","THBD&CLEC4A")


CLEC4A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = CLEC4A > 1)
CLEC4A_Mature_DC_Semimature_DC <- subset(CLEC4A_Mature_DC_Semimature_DC, subset = THBD < 1)
CLEC4A_Mature_DC_Semimature_DC_list <- as.data.frame(table(CLEC4A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
CLEC4A_Mature_DC_Semimature_DC_list[is.na(CLEC4A_Mature_DC_Semimature_DC_list)] <- 0
colnames(CLEC4A_Mature_DC_Semimature_DC_list) <- c("Cluster","CLEC4A")

list <- merge(THBD_Mature_DC_Semimature_DC_list, THBD_CLEC4A_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, CLEC4A_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_THBD_CLEC4A.csv")

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

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_THBD_CLEC4A_count.pdf", width = 5, height = 6)

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

ggsave("cluster_bargraph_DC_THBD_CLEC4A_proportion.pdf", width = 5, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC THBD vs IL10 co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("number")
#Control 10
#preTx_LS 11
#postTx_week12_LS 4

Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

THBD_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = THBD > 1)
THBD_Mature_DC_Semimature_DC <- subset(THBD_Mature_DC_Semimature_DC, subset = IL10 < 1)
THBD_Mature_DC_Semimature_DC_list <- as.data.frame(table(THBD_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
THBD_Mature_DC_Semimature_DC_list[is.na(THBD_Mature_DC_Semimature_DC_list)] <- 0
colnames(THBD_Mature_DC_Semimature_DC_list) <- c("Cluster","THBD")

THBD_IL10_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL10 > 1)
THBD_IL10_Mature_DC_Semimature_DC <- subset(THBD_IL10_Mature_DC_Semimature_DC, subset = THBD > 1)
THBD_IL10_Mature_DC_Semimature_DC_list <- as.data.frame(table(THBD_IL10_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
THBD_IL10_Mature_DC_Semimature_DC_list[is.na(THBD_IL10_Mature_DC_Semimature_DC_list)] <- 0
colnames(THBD_IL10_Mature_DC_Semimature_DC_list) <- c("Cluster","THBD&IL10")


IL10_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL10 > 1)
IL10_Mature_DC_Semimature_DC <- subset(IL10_Mature_DC_Semimature_DC, subset = THBD < 1)
IL10_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL10_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL10_Mature_DC_Semimature_DC_list[is.na(IL10_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL10_Mature_DC_Semimature_DC_list) <- c("Cluster","IL10")

list <- merge(THBD_Mature_DC_Semimature_DC_list, THBD_IL10_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, IL10_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_THBD_IL10.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL10","THBD&IL10", "THBD" )))
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
  ggtitle("Number of THBD and IL10 expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_THBD_IL10_count.pdf", width = 5, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("IL10","THBD&IL10", "THBD" )))
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
  ggtitle("Proportion of THBD and IL10 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_THBD_IL10_proportion.pdf", width = 5, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC ITGAX vs FCGR1A co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("number")
#Control 10
#preTx_LS 11
#postTx_week12_LS 4

Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

ITGAX_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = ITGAX > 1)
ITGAX_Mature_DC_Semimature_DC <- subset(ITGAX_Mature_DC_Semimature_DC, subset = FCGR1A < 1)
ITGAX_Mature_DC_Semimature_DC_list <- as.data.frame(table(ITGAX_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
ITGAX_Mature_DC_Semimature_DC_list[is.na(ITGAX_Mature_DC_Semimature_DC_list)] <- 0
colnames(ITGAX_Mature_DC_Semimature_DC_list) <- c("Cluster","ITGAX")

ITGAX_FCGR1A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = FCGR1A > 1)
ITGAX_FCGR1A_Mature_DC_Semimature_DC <- subset(ITGAX_FCGR1A_Mature_DC_Semimature_DC, subset = ITGAX > 1)
ITGAX_FCGR1A_Mature_DC_Semimature_DC_list <- as.data.frame(table(ITGAX_FCGR1A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
ITGAX_FCGR1A_Mature_DC_Semimature_DC_list[is.na(ITGAX_FCGR1A_Mature_DC_Semimature_DC_list)] <- 0
colnames(ITGAX_FCGR1A_Mature_DC_Semimature_DC_list) <- c("Cluster","ITGAX&FCGR1A")


FCGR1A_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = FCGR1A > 1)
FCGR1A_Mature_DC_Semimature_DC <- subset(FCGR1A_Mature_DC_Semimature_DC, subset = ITGAX < 1)
FCGR1A_Mature_DC_Semimature_DC_list <- as.data.frame(table(FCGR1A_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
FCGR1A_Mature_DC_Semimature_DC_list[is.na(FCGR1A_Mature_DC_Semimature_DC_list)] <- 0
colnames(FCGR1A_Mature_DC_Semimature_DC_list) <- c("Cluster","FCGR1A")

list <- merge(ITGAX_Mature_DC_Semimature_DC_list, ITGAX_FCGR1A_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, FCGR1A_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_ITGAX_FCGR1A.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FCGR1A","ITGAX&FCGR1A", "ITGAX" )))
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
  ggtitle("Number of ITGAX and FCGR1A expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_ITGAX_FCGR1A_count.pdf", width = 5, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FCGR1A","ITGAX&FCGR1A", "ITGAX" )))
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
  ggtitle("Proportion of ITGAX and FCGR1A expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_ITGAX_FCGR1A_proportion.pdf", width = 5, height = 6)

######################################################################################################
## Mature_DC and Semimature_DC IL10 vs TGFB1 co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("number")
#Control 10
#preTx_LS 11
#postTx_week12_LS 4

Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

Mature_DC_Semimature_DC <- subset(Round2, idents = c("Mature_DC","Semimature_DC" ))

IL10_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = IL10 > 1)
IL10_Mature_DC_Semimature_DC <- subset(IL10_Mature_DC_Semimature_DC, subset = TGFB1 < 1)
IL10_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL10_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL10_Mature_DC_Semimature_DC_list[is.na(IL10_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL10_Mature_DC_Semimature_DC_list) <- c("Cluster","IL10")

IL10_TGFB1_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = TGFB1 > 1)
IL10_TGFB1_Mature_DC_Semimature_DC <- subset(IL10_TGFB1_Mature_DC_Semimature_DC, subset = IL10 > 1)
IL10_TGFB1_Mature_DC_Semimature_DC_list <- as.data.frame(table(IL10_TGFB1_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
IL10_TGFB1_Mature_DC_Semimature_DC_list[is.na(IL10_TGFB1_Mature_DC_Semimature_DC_list)] <- 0
colnames(IL10_TGFB1_Mature_DC_Semimature_DC_list) <- c("Cluster","IL10&TGFB1")


TGFB1_Mature_DC_Semimature_DC <- subset(Mature_DC_Semimature_DC, subset = TGFB1 > 1)
TGFB1_Mature_DC_Semimature_DC <- subset(TGFB1_Mature_DC_Semimature_DC, subset = IL10 < 1)
TGFB1_Mature_DC_Semimature_DC_list <- as.data.frame(table(TGFB1_Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
TGFB1_Mature_DC_Semimature_DC_list[is.na(TGFB1_Mature_DC_Semimature_DC_list)] <- 0
colnames(TGFB1_Mature_DC_Semimature_DC_list) <- c("Cluster","TGFB1")

list <- merge(IL10_Mature_DC_Semimature_DC_list, IL10_TGFB1_Mature_DC_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, TGFB1_Mature_DC_Semimature_DC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(Mature_DC_Semimature_DC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster")

list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11

write.csv(list, file="Cellnumbers_Mature_DC_Semimature_DC_IL10_TGFB1.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("TGFB1","IL10&TGFB1", "IL10" )))
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
  ggtitle("Number of IL10 and TGFB1 expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_Mature_DC_Semimature_DC_IL10_TGFB1_count.pdf", width = 5, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("TGFB1","IL10&TGFB1", "IL10" )))
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
  ggtitle("Proportion of IL10 and TGFB1 expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_IL10_TGFB1_proportion.pdf", width = 5, height = 6)

######################################################################################################
## IL10 / TGFB1  Semimature_DC THBD / CLEC4A co-expression ####

#IL10 / TGFB1  Semimature_DC
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
IL10_Semimature_DC <- subset(IL10_Semimature_DC, subset = TGFB1 < 1)
levels(IL10_Semimature_DC)
new.cluster.ids <- c(  "IL10_Semimature_DC")
names(new.cluster.ids) <- levels(IL10_Semimature_DC)
IL10_Semimature_DC <- RenameIdents(IL10_Semimature_DC, new.cluster.ids)
levels(Idents(IL10_Semimature_DC))

IL10_TGFB1_Semimature_DC <- subset(Semimature_DC, subset = TGFB1 > 1)
IL10_TGFB1_Semimature_DC <- subset(IL10_TGFB1_Semimature_DC, subset = IL10 > 1)
levels(IL10_TGFB1_Semimature_DC)
new.cluster.ids <- c(   "IL10_TGFB1_Semimature_DC")
names(new.cluster.ids) <- levels(IL10_TGFB1_Semimature_DC)
IL10_TGFB1_Semimature_DC <- RenameIdents(IL10_TGFB1_Semimature_DC, new.cluster.ids)
levels(Idents(IL10_TGFB1_Semimature_DC))

TGFB1_Semimature_DC <- subset(Semimature_DC, subset = TGFB1 > 1)
TGFB1_Semimature_DC <- subset(TGFB1_Semimature_DC, subset = IL10 < 1)
levels(TGFB1_Semimature_DC)
new.cluster.ids <- c(    "TGFB1_Semimature_DC")
names(new.cluster.ids) <- levels(TGFB1_Semimature_DC)
TGFB1_Semimature_DC <- RenameIdents(TGFB1_Semimature_DC, new.cluster.ids)
levels(Idents(TGFB1_Semimature_DC))

IL10_TGFB1_Semimature_DC <- merge(x= IL10_Semimature_DC, y = c(IL10_TGFB1_Semimature_DC, TGFB1_Semimature_DC ), project = "IL10_TGFB1_Semimature_DC")

levels(Idents(IL10_TGFB1_Semimature_DC))
new.cluster.ids <- c("IL10_TGFB1_Semimature_DC",
                     "IL10_TGFB1_Semimature_DC",
                     "IL10_TGFB1_Semimature_DC")
names(new.cluster.ids) <- levels(IL10_TGFB1_Semimature_DC)
IL10_TGFB1_Semimature_DC <- RenameIdents(IL10_TGFB1_Semimature_DC, new.cluster.ids)
levels(Idents(IL10_TGFB1_Semimature_DC))
IL10_TGFB1_Semimature_DC[["ClusterNames_0.8"]] <- Idents(object = IL10_TGFB1_Semimature_DC)

# THBD / CLEC4A co-expression

THBD_IL10_TGFB1_Semimature_DC <- subset(IL10_TGFB1_Semimature_DC, subset = THBD > 1)
THBD_IL10_TGFB1_Semimature_DC <- subset(THBD_IL10_TGFB1_Semimature_DC, subset = CLEC4A < 1)
THBD_IL10_TGFB1_Semimature_DC_list <- as.data.frame(table(THBD_IL10_TGFB1_Semimature_DC@meta.data$ClusterNames_0.8))
THBD_IL10_TGFB1_Semimature_DC_list[is.na(THBD_IL10_TGFB1_Semimature_DC_list)] <- 0
colnames(THBD_IL10_TGFB1_Semimature_DC_list) <- c("Cluster","THBD")

THBD_CLEC4A_IL10_TGFB1_Semimature_DC <- subset(IL10_TGFB1_Semimature_DC, subset = CLEC4A > 1)
THBD_CLEC4A_IL10_TGFB1_Semimature_DC <- subset(THBD_CLEC4A_IL10_TGFB1_Semimature_DC, subset = THBD > 1)
THBD_CLEC4A_IL10_TGFB1_Semimature_DC_list <- as.data.frame(table(THBD_CLEC4A_IL10_TGFB1_Semimature_DC@meta.data$ClusterNames_0.8))
THBD_CLEC4A_IL10_TGFB1_Semimature_DC_list[is.na(THBD_CLEC4A_IL10_TGFB1_Semimature_DC_list)] <- 0
colnames(THBD_CLEC4A_IL10_TGFB1_Semimature_DC_list) <- c("Cluster","THBD&CLEC4A")


CLEC4A_IL10_TGFB1_Semimature_DC <- subset(IL10_TGFB1_Semimature_DC, subset = CLEC4A > 1)
CLEC4A_IL10_TGFB1_Semimature_DC <- subset(CLEC4A_IL10_TGFB1_Semimature_DC, subset = THBD < 1)
CLEC4A_IL10_TGFB1_Semimature_DC_list <- as.data.frame(table(CLEC4A_IL10_TGFB1_Semimature_DC@meta.data$ClusterNames_0.8))
CLEC4A_IL10_TGFB1_Semimature_DC_list[is.na(CLEC4A_IL10_TGFB1_Semimature_DC_list)] <- 0
colnames(CLEC4A_IL10_TGFB1_Semimature_DC_list) <- c("Cluster","CLEC4A")

list <- merge(THBD_IL10_TGFB1_Semimature_DC_list, THBD_CLEC4A_IL10_TGFB1_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, CLEC4A_IL10_TGFB1_Semimature_DC_list, by  ="Cluster",all=TRUE)
list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(IL10_TGFB1_Semimature_DC@meta.data$ClusterNames_0.8))
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

write.csv(list, file="Cellnumbers_IL10_TGFB1_Semimature_DC_THBD_CLEC4A.csv")

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

ggsave("cluster_bargraph_IL10_TGFB1_Semimature_DC_THBD_CLEC4A_count.pdf", width = 3, height = 6)

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

ggsave("cluster_bargraph_IL10_TGFB1_Semimature_DC_THBD_CLEC4A_proportion.pdf", width = 3, height = 6)





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
## IL23A  Semimature_DC ITGAX / FCGR1A co-expression ####

#IL23A  Semimature_DC
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

IL23A_Semimature_DC <- subset(Semimature_DC, subset = IL23A > 1)
levels(IL23A_Semimature_DC)
new.cluster.ids <- c(  "IL23A_Semimature_DC")
names(new.cluster.ids) <- levels(IL23A_Semimature_DC)
IL23A_Semimature_DC <- RenameIdents(IL23A_Semimature_DC, new.cluster.ids)
levels(Idents(IL23A_Semimature_DC))
IL23A_Semimature_DC[["ClusterNames_0.8"]] <- Idents(object = IL23A_Semimature_DC)


# ITGAX / FCGR1A co-expression
ITGAX_IL23A_Semimature_DC <- subset(IL23A_Semimature_DC, subset = ITGAX > 1)
ITGAX_IL23A_Semimature_DC <- subset(ITGAX_IL23A_Semimature_DC, subset = FCGR1A < 1)
ITGAX_IL23A_Semimature_DC_list <- as.data.frame(table(ITGAX_IL23A_Semimature_DC@meta.data$ClusterNames_0.8))
ITGAX_IL23A_Semimature_DC_list[is.na(ITGAX_IL23A_Semimature_DC_list)] <- 0
colnames(ITGAX_IL23A_Semimature_DC_list) <- c("Cluster","ITGAX")

ITGAX_FCGR1A_IL23A_Semimature_DC <- subset(IL23A_Semimature_DC, subset = FCGR1A > 1)
ITGAX_FCGR1A_IL23A_Semimature_DC <- subset(ITGAX_FCGR1A_IL23A_Semimature_DC, subset = ITGAX > 1)
ITGAX_FCGR1A_IL23A_Semimature_DC_list <- as.data.frame(table(ITGAX_FCGR1A_IL23A_Semimature_DC@meta.data$ClusterNames_0.8))
ITGAX_FCGR1A_IL23A_Semimature_DC_list[is.na(ITGAX_FCGR1A_IL23A_Semimature_DC_list)] <- 0
colnames(ITGAX_FCGR1A_IL23A_Semimature_DC_list) <- c("Cluster","ITGAX&FCGR1A")


FCGR1A_IL23A_Semimature_DC <- subset(IL23A_Semimature_DC, subset = FCGR1A > 1)
FCGR1A_IL23A_Semimature_DC <- subset(FCGR1A_IL23A_Semimature_DC, subset = ITGAX < 1)
FCGR1A_IL23A_Semimature_DC_list <- as.data.frame(table(FCGR1A_IL23A_Semimature_DC@meta.data$ClusterNames_0.8))
FCGR1A_IL23A_Semimature_DC_list[is.na(FCGR1A_IL23A_Semimature_DC_list)] <- 0
FCGR1A_IL23A_Semimature_DC_list$Freq <- 0 #0 is correct
colnames(FCGR1A_IL23A_Semimature_DC_list) <- c("Cluster","FCGR1A")

list <- merge(ITGAX_IL23A_Semimature_DC_list, ITGAX_FCGR1A_IL23A_Semimature_DC_list, by ="Cluster", all=TRUE)
list <- merge(list, FCGR1A_IL23A_Semimature_DC_list, by  ="Cluster",all=TRUE)
list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(IL23A_Semimature_DC@meta.data$ClusterNames_0.8))
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

write.csv(list, file="Cellnumbers_IL23A_Semimature_DC_ITGAX_FCGR1A.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FCGR1A","ITGAX&FCGR1A", "ITGAX" )))
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
  ggtitle("Number of ITGAX and FCGR1A expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_IL23A_Semimature_DC_ITGAX_FCGR1A_count.pdf", width = 3, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("FCGR1A","ITGAX&FCGR1A", "ITGAX" )))
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
  ggtitle("Proportion of ITGAX and FCGR1A expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_IL23A_Semimature_DC_ITGAX_FCGR1A_proportion.pdf", width = 3, height = 6)

######################################################################################################
