#####################################################################################################
## Perform an integrated analysis ####
DefaultAssay(Round0) <- "integrated"
colnames(Round0@meta.data)

# Run the standard workflow for visualization and clustering
Round0 <- ScaleData(Round0, verbose = FALSE)
Round0 <- RunPCA(Round0, npcs = 60, verbose = FALSE)

# t-SNE and Clustering
Round0 <- FindNeighbors(Round0, reduction = "pca", dims = 1:60)
Round0 <- FindClusters(Round0, resolution = 0.8)
Round0 <- RunUMAP(object = Round0, reductiosdfn = "pca", dims = 1:60)

######################################################################################################
## Dimensionality reduction ####
## ElbowPlot
pdf('./combined_data/ElbowPlot_Round0.pdf',width=10,height=10)
ElbowPlot(object = Round0, ndims = 60)
dev.off()

##  Perform linear dimensional reduction 
pdf('./combined_data/VizDimLoadings_Round0.pdf', width=8,height=7)
VizDimLoadings(Round0, dims = 1:2, reduction = "pca")
dev.off()

pdf('./combined_data/DimPlots_Round0.pdf', width=8,height=7)
DimPlot(Round0, reduction = "pca")
dev.off()

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
## Dot plot -1 to check clusters ####
Round0@assays
colnames(Round0@meta.data)
Round0@meta.data$integrated_snn_res.0.8
DefaultAssay(Round0)
levels(Idents(object = Round0))

Idents(object = Round0) <- ("integrated_snn_res.0.8")
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


Idents(object = Round0) <- fct_rev(Idents(object = Round0))

pdf('Dot plot_res0.8_initial.pdf', width=16,height=10)
DotPlot(object = Round0, features = features.plot) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=18))
dev.off()

Idents(object = Round0) <- fct_rev(Idents(object = Round0))

#####################################################################################################
## DEG without labeling ####
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


########################################################################
#  Dot plot with labels ####
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

########################################################################
# Visualization with label ####

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
