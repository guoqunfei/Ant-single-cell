setwd("00.Rdata")
rm(list = ls())

library(data.table)
library(Matrix)
library(Seurat)
library(magrittr) 
library(dplyr)    
library(cowplot)

if(FALSE){
load("MergeV2_Worker_4000_B1_B2_202003gene200_log.RData")
workercom <-ant.combined
rm(ant.combined)
rm(ant.anchors)

load("MergeV2_Queen_4000_B1_B2_202003gene200.RData")
queencom <-ant.combined
rm(ant.combined)
rm(ant.anchors)

load("MergeV2_Gyne_4000_B1_B2_202003gene200_log.RData")
gynecom <-ant.combined
rm(ant.combined)
rm(ant.anchors)

load("MergeV2_Male_4000_B1_B2_202003gene200.RData")
malecom <-ant.combined
rm(ant.combined)
rm(ant.anchors)


#queencom <- NormalizeData(queencom,verbose = FALSE)

queencom <- FindVariableFeatures(queencom, selection.method = "vst", nfeatures = 4000)
workercom <- FindVariableFeatures(workercom, selection.method = "vst", nfeatures = 4000)
gynecom <- FindVariableFeatures(gynecom, selection.method = "vst", nfeatures = 4000)
malecom <- FindVariableFeatures(malecom, selection.method = "vst", nfeatures = 4000)

workercom$stim <- "Worker"
queencom$stim <- "Queen"
gynecom$stim <- "Gyne"
malecom$stim <- "Male"

ant.anchors <- FindIntegrationAnchors(object.list = list(workercom,queencom,gynecom,malecom), dims = 1:20,anchor.features = 4000 )
ant.combined <- IntegrateData(anchorset = ant.anchors, dims = 1:20)

DefaultAssay(ant.combined) <- "integrated"

save(ant.anchors,ant.combined, file="Merge_log_inter_W_Q_G_M_4000_1.5_40_202004gene200.RData")

#}

load("Merge_log_inter_W_Q_G_M_4000_1.5_40_202004gene200.RData")
# Run the standard workflow for visualization and clustering
ant.combined <- ScaleData(ant.combined, verbose = FALSE)
ant.combined <- RunPCA(ant.combined, npcs = 25, verbose = FALSE)
# t-SNE and Clustering
ant.combined <- RunUMAP(ant.combined, reduction = "pca", dims = 1:20)
ant.combined <- RunTSNE(ant.combined, reduction = "pca", dims = 1:20)
ant.combined <- FindNeighbors(ant.combined, reduction = "pca", dims = 1:20)
ant.combined <- FindClusters(ant.combined, resolution = 1.5)
#ant.combined <- FindClusters(ant.combined, resolution = 1)
# Visualization

#pdf("MergeV3_Queen_Worker_Gyne_Male_B1B2_202003gene200_4000_20_0.8.pdf",height=10,width=20)

best_color<- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFE4E1","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#1B4400","#4FC601","#3B5DFF","#FF2F80","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#DDEFFF","#7B4F4B","#A1C299","#0AA6D8","#00A087FF","#4DBBD5FF","#E64B35FF","#3C5488FF","#F38400","#A1CAF1", "#C2B280","#848482","#E68FAC", "#0067A5","#F99379", "#604E97","#F6A600", "#B3446C","#DCD300","#882D17", "#8DB600","#654522", "#E25822", "#2B3D26","#191970","#000080","#6495ED","#1E90FF","#00BFFF","#00FFFF","#FF1493","#FF00FF","#A020F0","#63B8FF","#008B8B","#54FF9F","#00FF00","#76EE00","#FFF68F","Yellow1","Gold1","DarkGoldenrod4","#FF6A6A","#FF8247","#FFA54F","#FF7F24","#FF3030","#FFA500","#FF7F00","#FF7256","#FF6347","#FF4500","#FF1493","#FF6EB4","#EE30A7","#8B008B")

p1 <- DimPlot(ant.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(ant.combined, reduction = "umap", label = TRUE,cols = best_color)
p3 <- DimPlot(ant.combined, reduction = "tsne", group.by = "stim")
p4 <- DimPlot(ant.combined, reduction = "tsne", label = TRUE,cols = best_color)
cellnumber<- dim(as.data.frame(ant.combined@reductions$umap@cell.embeddings))
print(paste0("cellnumber:  ", cellnumber[1]))

pdf("Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_umap.pdf",height=10,width=22)
plot_grid(p1, p2)
dev.off()
pdf("Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_tsne.pdf",height=10,width=22)
plot_grid(p3, p4)
dev.off()
pdf("Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_tsne.pdf",height=10,width=12)
p4

dev.off()


DefaultAssay(ant.combined) <- "RNA"

ant.combined.markers <- FindAllMarkers(ant.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ant.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top20 <- ant.combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

top40 <- ant.combined.markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_logFC)


write.table(top20,"Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_top20.xls",sep = "\t",row.names = F,col.names =TRUE)
write.table(top40,"Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_top40.xls",sep = "\t",row.names = F,col.names =TRUE)
write.table(ant.combined.markers,"Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_all.xls",sep = "\t",row.names = F,col.names =TRUE)


save(ant.combined.markers,ant.combined, file="Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_RNA.RData")

dev.off()

}

load("Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_RNA.RData")


cluster.averages <- AverageExpression(ant.combined)
ant.esper_RNA <- cluster.averages[["RNA"]]
cor_test <- cor(ant.esper_RNA,ant.esper_RNA)
library(pheatmap)

pdf("Merge_log_inter_W_Q_G_M_4000_1.5_25_202004gene200_totalgeneAndvargene.pdf",height=15,width=16)
pheatmap(cor_test)
var_featuregene <- as.data.frame(ant.combined@assays$integrated@var.features)
colnames(var_featuregene) <- "vargene"

var_ant.expr <- data.frame(matrix(nrow = dim(var_featuregene)[1],ncol=dim(ant.esper_RNA)[2]))
rownames(var_ant.expr)<- var_featuregene$vargene
colnames(var_ant.expr)<-colnames(ant.esper_RNA)

for (i in rownames(var_ant.expr)) {
	  if(i%in%rownames(ant.esper_RNA))var_ant.expr[i,]<- ant.esper_RNA[i,]
}

var_ant.expr <-  na.omit(var_ant.expr)
cor_var_total_ant_ant.exp <- cor(var_ant.expr,var_ant.expr)
pheatmap(cor_var_total_ant_ant.exp)

dev.off()





