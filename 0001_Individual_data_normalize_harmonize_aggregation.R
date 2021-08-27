#load library ####
# We used Seurat R package (version 4.0, https://satijalab.org/seurat/) 
# installed in R (version 4.0.2, https://www.r-project.org/) for the downstream single-cell RNA sequencing data analy
library(Seurat)
######################################################################################################3
#Function to filter out cells with >25% mitochondrial genes 
#and ubiquitously expressed ribosomal protein-coding (RPS and RPL) and MALAT noncoding RNA, miRNA, and snRNA genes 
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
}

#####################################################################################
# Individual single-cell data loading & quality control ####
# Genes expressed in <3 cells, and cells with <100 or > 5,000 genes, 
# and a mitochondrial gene percentage of >25% were filtered out to eliminate partial cells and doublets. 
# Ubiquitously expressed ribosomal protein-coding (RPS and RPL) and MALAT noncoding RNA, miRNA, and snoRNA genes 
# were also excluded from analysis. 
# Seurat objects were created, followed by normalizing data, scaling data, and finding variable 2,000 genes. 

# This strategy has been used for the coauthors' previous skin single-cell RNA sequencing studies
# 1. Der, E. et al. Tubular cell and keratinocyte single-cell transcriptomics applied to lupus nephritis reveal type I IFN and fibrosis relevant pathways. Nature immunology 20, 915-927 (2019)
# 2. He, H. et al. Single-cell transcriptome analysis of human skin identifies novel fibroblast subpopulation and enrichment of immune subsets in atopic dermatitis. The Journal of allergy and clinical immunology (2020).
# 3. Kim J, Lee J, Kim HJ, Kameyama N, Nazarian R, Der E, et al. Single-cell transcriptomics applied to emigrating cells from psoriasis elucidate pathogenic versus regulatory immune cell subsets. J Allergy Clin Immunol 2021

# All the single-cell data used for the analysis are publicly availalbe on NCBI Gene Expression Omnibus (GEO) with phenotype information
# The following is the example of data loading "Control01".

Control01 <- Read10X(data.dir = "./Control01")

Control01 <- CreateSeuratObject(counts = Control01, project = "Control01", min.cells = 3, min.features = 100)
Control01[["percent.mt"]] <- PercentageFeatureSet(object = Control01, pattern = "^MT-")
Control01@meta.data$stim <- as.character(phenotype[grep('Control01', phenotype$ID),2])
Control01@meta.data$number <- as.character(phenotype[grep('Control01', phenotype$ID),1])

Control01_morethan25percentMT <- subset(Control01, subset = percent.mt > 25)
Control01_morethan25percentMTbc <- colnames(Control01_morethan25percentMT@assays$RNA@data)

Control01 <- Read10X(data.dir = "./Control01")
Control01<-remove.genes(Control01, Control01_morethan25percentMTbc)

Control01 <- CreateSeuratObject(counts = Control01, project = "Control01", min.cells = 3, min.features = 100)
Control01[["percent.mt"]] <- PercentageFeatureSet(object = Control01, pattern = "^MT-")
Control01 <- subset(Control01, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
Control01@meta.data$stim <- as.character(phenotype[grep('Control01', phenotype$ID),2])
Control01@meta.data$number <- as.character(phenotype[grep('Control01', phenotype$ID),1])
rm(Control01_morethan25percentMTbc)

Control01 <- NormalizeData(Control01)
Control01 <- FindVariableFeatures(Control01, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Control01)
Control01 <- ScaleData(Control01, features = all.genes)





######################################################################################################3
#Data Integration
# We used Seurat R package (version 4.0, https://satijalab.org/seurat/) for the single-cell data integration. 
# We first merged single-cell data of samples with the identical reagent kit version 
# and the identical sequencer: 
#(1) psoriasis samples with reagent kit V2.0 - 10 samples, \
#(2) psoriasis samples with reagent kit V3.0 - 3 samples, and 
#(3) control samples with reagent kit V2.0 - 5 samples. 

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
# Merged data sets are individualized normalized
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
# To harmonize merged groups into a single dataset without batch effects, 
# correspondences between cells in three merged datasets were identified by the FIndIntegrationAnchors function, 
# and used for data integration with the IntegratedData function as detailed by Butler et al. 
immune.anchors <- FindIntegrationAnchors(object.list = list(V2.0,V3.0,V3.1_NextSeq ,V3.1_NovaSeq), dims = 1:30)
Round0 <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
