######################################################################################################
## KC IL36G vs DEFB4A co-expression ####
Round2 <- Round1
colnames(Round2@meta.data)
DefaultAssay(Round2) <- "RNA"
colnames(Round2@meta.data)

Idents(object = Round2) <- ("Serial")
levels(Round2)
Round2 <- subset(Round2, idents = c("Control" , "Psoriasis_PreTx_LS_week0" ,"Psoriasis_PostTx_week12"  ))
Idents(object = Round2) <- ("ClusterNames_0.8")
levels(Round2)

KC <- subset(Round2, idents = c("KC-S.Corneum"   ,  "KC-S.Granulosum" , "KC-S.Spinosum"  ,  "KC-S.Basale"  ))

DefaultAssay(KC) <- "RNA"

IL36G_KC <- KC

IL36G_KC <- subset(KC, subset = IL36G > 1)
IL36G_KC <- subset(IL36G_KC, subset = DEFB4A < 1)
IL36G_KC_list <- as.data.frame(table(IL36G_KC@meta.data$ClusterNames_0.8_Serial))
IL36G_KC_list[is.na(IL36G_KC_list)] <- 0
colnames(IL36G_KC_list) <- c("Cluster","IL36G")

IL36G_DEFB4A_KC <- subset(KC, subset = DEFB4A > 1)
IL36G_DEFB4A_KC <- subset(IL36G_DEFB4A_KC, subset = IL36G > 1)
IL36G_DEFB4A_KC_list <- as.data.frame(table(IL36G_DEFB4A_KC@meta.data$ClusterNames_0.8_Serial))
IL36G_DEFB4A_KC_list[is.na(IL36G_DEFB4A_KC_list)] <- 0
colnames(IL36G_DEFB4A_KC_list) <- c("Cluster","IL36G&DEFB4A")


DEFB4A_KC <- subset(KC, subset = DEFB4A > 1)
DEFB4A_KC <- subset(DEFB4A_KC, subset = IL36G < 1)
DEFB4A_KC_list <- as.data.frame(table(DEFB4A_KC@meta.data$ClusterNames_0.8_Serial))
DEFB4A_KC_list[is.na(DEFB4A_KC_list)] <- 0
colnames(DEFB4A_KC_list) <- c("Cluster","DEFB4A")

list <- merge(IL36G_KC_list, IL36G_DEFB4A_KC_list, by ="Cluster", all=TRUE)
list <- merge(list, DEFB4A_KC_list, by  ="Cluster",all=TRUE)

list[is.na(list)] <- 0


total_cell_count <- as.data.frame(table(KC@meta.data$ClusterNames_0.8_Serial))
total_cell_count[is.na(total_cell_count)] <- 0
colnames(total_cell_count) <- c("Cluster","total_cell_number")
list <- merge(list, total_cell_count, by  ="Cluster",all=TRUE)

#view(list)
list[,c(6,7,8)] <- list[,c(2,3,4)]
colnames(list)[c(6,7,8)] <- paste(colnames(list[,c(2,3,4)]), "_Proportion", sep = "")
list[,c(6,7,8)] <- list[,c(6,7,8)]/list[,5]*100

#Control sample number 10
#PosstTx sample number 4
#Pretx sample number 11
list[grep("Control", list$Cluster),c(2:4) ] <- list[grep("Control", list$Cluster),c(2:4) ]/10
list[grep("PostTx", list$Cluster),c(2:4) ] <- list[grep("PostTx", list$Cluster),c(2:4) ]/4
list[grep("PreTx", list$Cluster),c(2:4) ] <- list[grep("PreTx", list$Cluster),c(2:4) ]/11
list[is.na(list)] <- 0

write.csv(list, file="Cellnumbers_KC_IL36G_DEFB4A.csv")

#Bar graph for clusters (count)
list_for_ggplot <- list[c(1:4)]
list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("DEFB4A","IL36G&DEFB4A", "IL36G" )))
levels(list_for_ggplot$Cluster)
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c(  "KC-S.Corneum_Control"        ,            
                                                                        "KC-S.Corneum_Psoriasis_PreTx_LS_week0"   ,
                                                                        "KC-S.Corneum_Psoriasis_PostTx_week12" ,
                                                                        "KC-S.Granulosum_Control" ,
                                                                        "KC-S.Granulosum_Psoriasis_PreTx_LS_week0",
                                                                        "KC-S.Granulosum_Psoriasis_PostTx_week12" ,
                                                                        "KC-S.Spinosum_Control" ,
                                                                        "KC-S.Spinosum_Psoriasis_PreTx_LS_week0" ,
                                                                        "KC-S.Spinosum_Psoriasis_PostTx_week12"   , "KC-S.Basale_Control"     ,
                                                                        "KC-S.Basale_Psoriasis_PreTx_LS_week0" ,
                                                                        "KC-S.Basale_Psoriasis_PostTx_week12"     )))


ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Number of IL36G and DEFB4A expressing cells per sample")+
  xlab("Clusters") + ylab("Number of cells per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red")) 

ggsave("cluster_bargraph_KC_IL36G_DEFB4A_count.pdf", width = 6, height = 6)

#Bar graph for clusters (proportion)
list_for_ggplot <- list[c(1,6:8)]
colnames(list_for_ggplot)[c(2,3,4)] <-  sub("_Proportion", "", colnames(list_for_ggplot)[c(2,3,4)])

list_for_ggplot <- reshape2::melt(list_for_ggplot)
list_for_ggplot$variable <- factor(list_for_ggplot$variable, levels = (c("DEFB4A","IL36G&DEFB4A", "IL36G" )))
list_for_ggplot$Cluster<- factor(list_for_ggplot$Cluster, levels = (c(  "KC-S.Corneum_Control"        ,            
                                                                        "KC-S.Corneum_Psoriasis_PreTx_LS_week0"   ,
                                                                        "KC-S.Corneum_Psoriasis_PostTx_week12" ,
                                                                        "KC-S.Granulosum_Control" ,
                                                                        "KC-S.Granulosum_Psoriasis_PreTx_LS_week0",
                                                                        "KC-S.Granulosum_Psoriasis_PostTx_week12" ,
                                                                        "KC-S.Spinosum_Control" ,
                                                                        "KC-S.Spinosum_Psoriasis_PreTx_LS_week0" ,
                                                                        "KC-S.Spinosum_Psoriasis_PostTx_week12"   , "KC-S.Basale_Control"     ,
                                                                        "KC-S.Basale_Psoriasis_PreTx_LS_week0" ,
                                                                        "KC-S.Basale_Psoriasis_PostTx_week12"    )))

ggplot(data = list_for_ggplot, aes(Cluster, value, fill = variable)) +
  geom_bar(stat = 'identity') +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Proportion of IL36G and DEFB4A expressing cells")+
  xlab("Clusters") + ylab("Proportion of cells (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black")) +
  scale_fill_manual("Expression", values = c("green","yellow", "red"))

ggsave("cluster_bargraph_DC_IL36G_DEFB4A_proportion.pdf", width = 6, height = 6)
