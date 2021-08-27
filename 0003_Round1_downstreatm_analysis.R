######################################################################################################
## Combine clustes ####
Round1 <- Round0
colnames(Round1@meta.data)
DefaultAssay(Round1) <- "RNA"
Idents(object = Round1) <- ("ClusterNames_0.8")
levels(Round1)

Round1 <- subset(Round1, idents = c( "CD4_T_cell.02" , "Doublet.01"   ,    "Doublet.02"),invert=TRUE) #Remove doublet clusters for downstream analysis

#Combine clusters wiht identical immune cell subset
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


########################################################################
# Dot plot with labels ####
Idents(object = Round1) <- ("ClusterNames_0.8")
DefaultAssay(Round1) <- "RNA"

features.plot <- c(  "KLRB1"   , "GNLY" ,  
                     "CD3D" ,    "TRAC" ,    "TRBC1" ,
                     "CD8A"   ,  "CD8B"   ,
                     "IFNG","IL17A","IL17F",
                     "GZMK"    ,"GZMH"   ,
                     "TIGIT" ,   "IL2RA"   ,  "FOXP3"   , "CTLA4" ,  
                     "LAMP3"  ,  "CCL22"  , "LY75"  ,   "CIITA"  ,  "CD40"  ,  
                     "HLA-DQA1","HLA-DQB1", "HLA-DRB1", "HLA-DRA" , "HLA-DRB5" ,
                     "LYZ"     ,  "CD14"   ,  "CD163"   ,  "THBD"    ,"IL10",
                     "DCT"    ,  "TYRP1"  ,  "MLANA"   ,
                     "SPRR2A"  , "SPRR2D"  , 
                     "FABP5"  ,
                     "KRT10"   , "KRT1"   ,"CCL20"  ,  
                     "KRT14"  ,  "KRT5"  ,
                     "KRT6B"   , "KRT17"   , "KRT16"   ,
                     "CCL21"  ,  "LYVE1"   , 
                    "DCN"    ,  "COL1A1" ,  "DCD"   ,   "SCGB2A2" )

Idents(object = Round1) <- fct_rev(Idents(object = Round1))

pdf('Dot plot_Round1_labeled01.pdf', width=14,height=5)
DotPlot(object = Round1, features = features.plot) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.49)) +
  theme(axis.text=element_text(size=12))
dev.off()

Idents(object = Round1) <- fct_rev(Idents(object = Round1))

#######################################################################
# Visualization with label ####

colnames(Round1@meta.data)
DefaultAssay(Round1)
DefaultAssay(Round1) <- "RNA"

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


###################################################################################################
# Average expression - Dendritic cells ####
Idents(object = Round1) <- ("ClusterNames_0.8")
DCs <- subset(Round1, idents = c("Mature_DC" ,   "Semimature_DC" ))

#Remove Non-responder (non-responder is presented in a separate analysis)
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


###################################################################################################
# Average expression - Dendritic cells : non-responder serial ####
setwd("~/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")
setwd("D:/Dropbox/NSK_0937/NSK_0937_Experiment_Data/NSK_0937_scRNAseq/NSK_0937_scRNAseq_Data_processed/Psoriasis_posttx_05.01.2021/Round1")

Idents(object = Round1) <- ("ClusterNames_0.8")
DCs <- subset(Round1, idents = c("Mature_DC" ,   "Semimature_DC" ))

#Remove Non-responder (non-responder is presented in a separate analysis)
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
               "TIGIT"  ,  "FOXP3" ,  "IL2RA" 
                   )


marker_shortened <- c(
            "CD3D","KLRB1", "GZMK", "IL26"  , "IL17A"  , "IL17F" ,
            "FOXP3" ,  "IL2RA" 
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
   scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis08_responder_Heatmap02.pdf', width=7, height=8)
DoHeatmap(object = T_cells.cluster.averages,features = marker, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
   scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis08_responder_shortened_Heatmap01.pdf', width=7, height=2.5)
DoHeatmap(object = T_cells.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
   scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()


pdf('T_cells.averages_Psoriasis08_responder_shortened_Heatmap02.pdf', width=7, height=10)
DoHeatmap(object = T_cells.cluster.averages,features = marker_shortened, raster = FALSE, size = 3, group.bar = TRUE, slot = "scale.data", draw.lines=FALSE)      +
   scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", "#fbfcbd", "#ff0000"))(256))
dev.off()
