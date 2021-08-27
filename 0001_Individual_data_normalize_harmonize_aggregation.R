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
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(reticulate)
library(scales)
library(forcats)
library(cowplot)
#install.packages("googlesheets4")
library(googlesheets4)
library(patchwork)
######################################################################################################3
#Start ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
rm(list = ls())
getwd()
list.files()

dir.create("./individual_data_analysis")
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021/individual_data_analysis")

#Read single cell data summary
phenotype <- read_sheet("https://https://docs.google.com/spreadsheets/d/1_cVcmP7AniLrhIE9XWw-T-_BqQXIkG8wi7wl1Rb3kl4/edit#gid=0&fvid=109363522")
phenotype <- phenotype[grep("Yes", phenotype$Third_project), ]

phenotype <- phenotype[order(phenotype$Rcoding_Number_05.01.2021),]

##Remove cluster *** barcodes and HighMT bc
remove.genes<-function(sample, BCtoremove){
  sample <- sample[, !colnames(sample) %in% BCtoremove]
  ubiq <- c(rownames(sample[grep("^MT-", rownames(sample)), ]),
            rownames(sample[grep("^RPS", rownames(sample)), ]),
            rownames(sample[grep("^RPL", rownames(sample)), ]),
            rownames(sample[grep("^RBP", rownames(sample)), ]),
            rownames(sample[grep("^MIR", rownames(sample)), ]),
            rownames(sample[grep("^SNOR", rownames(sample)), ]),
            rownames(sample[grep("^MTRNR", rownames(sample)), ]),
            rownames(sample[grep("^MALAT1", rownames(sample)), ]))
  sample <- sample[!rownames(sample) %in% ubiq, ]
  return(sample)
}## Remove genes by Mitochondrial percentage

#####################################################################################
### Load data  ####
#rm(list = ls()[which("phenotype" != ls())])
i=1

for (i in 1:nrow(phenotype)){
setwd("/Volumes/JKIM_20TB/Single_cell_analysis_GEO/Psoriasis_posttx_05.01.2021/individual_data")
dir.create("./10XRead")
new.folder <- "/Volumes/JKIM_20TB/Single_cell_analysis_GEO/Psoriasis_posttx_05.01.2021/individual_data/10XRead"

copy <- file.copy(paste(phenotype$Rcoding_Number_05.01.2021[i],"_barcodes.tsv.gz",sep=""), new.folder)
copy <- file.copy(paste(phenotype$Rcoding_Number_05.01.2021[i],"_features.tsv.gz",sep=""), new.folder)
copy <- file.copy(paste(phenotype$Rcoding_Number_05.01.2021[i],"_matrix.mtx.gz",sep=""), new.folder)

setwd("/Volumes/JKIM_20TB/Single_cell_analysis_GEO/Psoriasis_posttx_05.01.2021/individual_data/10XRead")
file.rename(paste(phenotype$Rcoding_Number_05.01.2021[i],"_barcodes.tsv.gz",sep=""),"barcodes.tsv.gz" )
file.rename(paste(phenotype$Rcoding_Number_05.01.2021[i],"_features.tsv.gz",sep=""),"features.tsv.gz" )
file.rename(paste(phenotype$Rcoding_Number_05.01.2021[i],"_matrix.mtx.gz",sep=""),"matrix.mtx.gz" )

#
assign("data",Read10X(data.dir = new.folder))

data <- CreateSeuratObject(counts = data, project = phenotype$Rcoding_Number_05.01.2021[i], min.cells = 3, min.features = 100)
data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")
data@meta.data$stim <- phenotype$Rcoding_Category[i]
data@meta.data$number <- phenotype$Rcoding_Number_05.01.2021[i]

pdf(paste("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021/individual_data_analysis/Vnplot_", phenotype$Rcoding_Number_05.01.2021[i], ".pdf",sep=""),width=10,height=5) ## Vnplot 
print(VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

plot1 <- FeatureScatter(object = data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(paste("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021/individual_data_analysis/FeatureScatter_", phenotype$Rcoding_Number_05.01.2021[i], ".pdf",sep=""),width=12,height=5)
print(plot1 + plot2)
dev.off()

data_morethan25percentMT <- NA
data_morethan25percentMTbc <- NA
tryCatch(data_morethan25percentMT <- subset(data, subset = percent.mt > 25), error=function(e) NULL)
tryCatch(data_morethan25percentMTbc <- colnames(data_morethan25percentMT@assays$RNA@data), error=function(e) NULL)
rm(data)
rm(data_morethan25percentMT)


#
assign("data",Read10X(data.dir = new.folder))
data<-remove.genes(data, data_morethan25percentMTbc)

data <- CreateSeuratObject(counts = data, project = phenotype$Rcoding_Number_05.01.2021[i], min.cells = 3, min.features = 100)
data[["percent.mt"]] <- PercentageFeatureSet(object = data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
data@meta.data$stim <- phenotype$Rcoding_Category[i]
data@meta.data$number <- phenotype$Rcoding_Number_05.01.2021[i]
rm(data_morethan25percentMTbc)

data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)

assign(phenotype$Rcoding_Number_05.01.2021[i],data)
rm(all.genes)
rm(data)
rm(plot1)
rm(plot2)

unlink("/Volumes/JKIM_20TB/Single_cell_analysis_GEO/Psoriasis_posttx_05.01.2021/individual_data/10XRead", recursive = T)
}



######################################################################################################3
##Merge data set ####
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")
dir.create("combined_data")
getwd()
list.files()

#V2.0 & NextSeq - 13
list_V2.0 <- phenotype[phenotype$Chemistry == "V2.0",]
list_V2.0 <- list_V2.0$Rcoding_Number_05.01.2021
as.factor(list_V2.0)
length(list_V2.0)

V2.0 <- merge(x = Control01, y = c(
  Control02          ,          Control03       ,             Control04         ,           Control05    ,               
  Psoriasis01_preTx_LS    ,     Psoriasis02_preTx_LS   ,      Psoriasis03_preTx_LS    ,     Psoriasis04_postTx_week12_LS, Psoriasis04_preTx_LS        ,
  Psoriasis05_postTx_week12_LS, Psoriasis05_preTx_LS   ,      Psoriasis06_preTx_LS   
    ), 
  add.cell.ids = list_V2.0, 
  project = "V2.0")


#V3.0 & NextSeq - 3
list_V3.0 <-phenotype[phenotype$Chemistry == "V3.0",]
list_V3.0 <- list_V3.0$Rcoding_Number_05.01.2021
as.factor(list_V3.0)
length(list_V3.0)

V3.0 <- merge(x = Psoriasis07_preTx_LS, y = c(Psoriasis08_preTx_LS ,Psoriasis09_preTx_LS
), 
add.cell.ids = list_V3.0
, project = "V3.0")


#V3.1 & NextSeq - 8
list_V3.1_NextSeq <- phenotype[phenotype$Chemistry == "V3.1" & phenotype$Sequencer == "NextSeq",]
list_V3.1_NextSeq <- list_V3.1_NextSeq$Rcoding_Number_05.01.2021
as.factor(list_V3.1_NextSeq)
length(list_V3.1_NextSeq)

V3.1_NextSeq <- merge(x = Control06,  y = c(Psoriasis07_postTx_week12_LS, Psoriasis07_postTx_week48_LS, Psoriasis07_preTx_LS_02,
                                            Psoriasis07_preTx_NL_02, Psoriasis08_postTx_week12_LS ,Psoriasis08_postTx_week24_LS ,Psoriasis10_preTx_LS
                                           
), 
add.cell.ids = list_V3.1_NextSeq
, project = "V3.1_NextSeq")

#V3.1 & NovaSeq - 4
list_V3.1_NovaSeq <- phenotype[phenotype$Chemistry == "V3.1" & phenotype$Sequencer == "NovaSeq",]
list_V3.1_NovaSeq <- list_V3.1_NovaSeq$Rcoding_Number_05.01.2021
as.factor(list_V3.1_NovaSeq)
length(list_V3.1_NovaSeq)

V3.1_NovaSeq <- merge(x = Control07,  y = c(Control08  ,  Control09 ,Control10
), 
add.cell.ids = list_V3.1_NovaSeq
, project = "V3.1_NovaSeq")

######################################################################################################
### Data normalization
V2.0 <- NormalizeData(V2.0)
V2.0 <- FindVariableFeatures(V2.0, selection.method = "vst", nfeatures = 2000)

V3.0 <- NormalizeData(V3.0)
V3.0 <- FindVariableFeatures(V3.0, selection.method = "vst", nfeatures = 2000)

V3.1_NextSeq <- NormalizeData(V3.1_NextSeq)
V3.1_NextSeq<- FindVariableFeatures(V3.1_NextSeq, selection.method = "vst", nfeatures = 2000)

V3.1_NovaSeq <- NormalizeData(V3.1_NovaSeq)
V3.1_NovaSeq<- FindVariableFeatures(V3.1_NovaSeq, selection.method = "vst", nfeatures = 2000)

######################################################################################################
### Data inegration
immune.anchors <- FindIntegrationAnchors(object.list = list(V2.0,V3.0,V3.1_NextSeq ,V3.1_NovaSeq), dims = 1:30)
Round0 <- IntegrateData(anchorset = immune.anchors, dims = 1:30)


######################################################################################################
### Perform an integrated analysis
setwd("~/Dropbox/10X_JKIM/aggregation/Psoriasis_posttx_05.01.2021")

Round0@assays
DefaultAssay(Round0) <- "integrated"

## Vnplot 
pdf('./combined_data/Vnplot_Round0.pdf',width=30,height=5)
VlnPlot(object = Round0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
dev.off()

#FeatureScatter
plot1 <- FeatureScatter(object = Round0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = Round0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf('./combined_data/FeatureScatter_Round0.pdf',width=20,height=5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

save(Round0, file = "./combined_data/Round0_integrated.Rda")



