library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)


Mpha <- readRDS("Mpha.rds")

ant.sub.sc <- as.SingleCellExperiment(Mpha,assay="integrated")

ant.milo <- Milo(ant.sub.sc)

# d should be the paramaters used in seurat
ant.milo <- buildGraph(ant.milo, k = 50, d = 20, reduced.dim = "PCA")
# @int_colData@listData$reducedDims
ant.milo <- makeNhoods(ant.milo, prop = 0.1, k = 50, d=20, refined = TRUE, reduced_dims = "PCA")

#ant.milo$orig.ident<-paste(ant.milo$caste, ant.milo$experiment, sep=".")

ant.milo <- countCells(ant.milo, meta.data = as.data.frame(colData(ant.milo)), sample="rep")
ant.design <- data.frame(colData(ant.milo))[,c("rep", "caste")]

ant.design$caste <- as.factor(ant.design$caste)
ant.design <- distinct(ant.design)
rownames(ant.design) <- ant.design$rep

ant.milo <- calcNhoodDistance(ant.milo, d=20, reduced.dim = "PCA")

da_results <- testNhoods(ant.milo, design = ~ caste, design.df = ant.design, reduced.dim="PCA")

ant.milo <- buildNhoodGraph(ant.milo)
umap_pl <- plotReducedDim(ant.milo, colour_by = "caste",text_by="seurat_clusters",text_size=3, point_size=0.5,dimred = "UMAP") + guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(ant.milo, da_results, layout="UMAP", alpha=0.1)
umap_pl + nh_graph_pl + plot_layout(guides="collect")

da_results <- annotateNhoods(ant.milo, da_results, coldata_col = "seurat_clusters")
ggplot(da_results, aes(seurat_clusters_fraction)) + geom_histogram(bins=50)
da_results$seurat_clusters <- ifelse(da_results$seurat_clusters_fraction < 0.7, "Mixed", da_results$seurat_clusters)
plotDAbeeswarm(da_results, group.by = "seurat_clusters")

cluster_sort <- c("8","14","18","23","26","31","33","35","5","13","21","37","6","15","16","20","25","27","42","11","34","38","0","1","2","3","4","7","9","10","12","19","28","30","41","17","36","40","22","29","39","24","32","Mixed")

QueenGyne_results <- testNhoods(ant.milo, design = ~ 0 + caste, design.df = ant.design, reduced.dim="PCA", model.contrasts = "casteQueen - casteGyne")
QueenGyne_results <- annotateNhoods(ant.milo, QueenGyne_results, coldata_col = "seurat_clusters")
QueenGyne_results$seurat_clusters <- ifelse(QueenGyne_results$seurat_clusters_fraction < 0.7, "Mixed", QueenGyne_results$seurat_clusters)
QueenGyne_results$seurat_clusters <- factor(QueenGyne_results$seurat_clusters, levels=cluster_sort)
plotDAbeeswarm(QueenGyne_results, group.by = "seurat_clusters") + ggsave("QueenGyne.DAbeeswarm.pdf", height=15, width=5)
write.table(QueenGyne_results, "QueenGyne.rep.res.txt", sep="\t", quote=F)

QueenMale_results <- testNhoods(ant.milo, design = ~ 0 + caste, design.df = ant.design, reduced.dim="PCA", model.contrasts = "casteQueen - casteMale")
QueenMale_results <- annotateNhoods(ant.milo, QueenMale_results, coldata_col = "seurat_clusters")
QueenMale_results$seurat_clusters <- ifelse(QueenMale_results$seurat_clusters_fraction < 0.7, "Mixed", QueenMale_results$seurat_clusters)
QueenMale_results$seurat_clusters <- factor(QueenMale_results$seurat_clusters, levels=cluster_sort)
plotDAbeeswarm(QueenMale_results, group.by = "seurat_clusters") + ggsave("QueenMale.DAbeeswarm.pdf", height=15, width=5)
write.table(QueenMale_results, "QueenMale.rep.res.txt", sep="\t", quote=F)

WorkerQueen_results <- testNhoods(ant.milo, design = ~ 0 + caste, design.df = ant.design, reduced.dim="PCA", model.contrasts = "casteWorker - casteQueen")
WorkerQueen_results <- annotateNhoods(ant.milo, WorkerQueen_results, coldata_col = "seurat_clusters")
WorkerQueen_results$seurat_clusters <- ifelse(WorkerQueen_results$seurat_clusters_fraction < 0.7, "Mixed", WorkerQueen_results$seurat_clusters)
WorkerQueen_results$seurat_clusters <- factor(WorkerQueen_results$seurat_clusters, levels=cluster_sort)
plotDAbeeswarm(WorkerQueen_results, group.by = "seurat_clusters") + ggsave("WorkerQueen.DAbeeswarm.pdf", height=15, width=5)
write.table(WorkerQueen_results, "WorkerQueen.rep.res.txt", sep="\t", quote=F)

WorkerGyne_results <- testNhoods(ant.milo, design = ~ 0 + caste, design.df = ant.design, reduced.dim="PCA", model.contrasts = "casteWorker - casteGyne")
WorkerGyne_results <- annotateNhoods(ant.milo, WorkerGyne_results, coldata_col = "seurat_clusters")
WorkerGyne_results$seurat_clusters <- ifelse(WorkerGyne_results$seurat_clusters_fraction < 0.7, "Mixed", WorkerGyne_results$seurat_clusters)
WorkerGyne_results$seurat_clusters <- factor(WorkerGyne_results$seurat_clusters, levels=cluster_sort)
plotDAbeeswarm(WorkerGyne_results, group.by = "seurat_clusters") + ggsave("WorkerGyne.DAbeeswarm.pdf", height=15, width=5)
write.table(WorkerGyne_results, "WorkerGyne.rep.res.txt", sep="\t", quote=F)

WorkerMale_results <- testNhoods(ant.milo, design = ~ 0 + caste, design.df = ant.design, reduced.dim="PCA", model.contrasts = "casteWorker - casteMale")
WorkerMale_results <- annotateNhoods(ant.milo, WorkerMale_results, coldata_col = "seurat_clusters")
WorkerMale_results$seurat_clusters <- ifelse(WorkerMale_results$seurat_clusters_fraction < 0.7, "Mixed", WorkerMale_results$seurat_clusters)
WorkerMale_results$seurat_clusters <- factor(WorkerMale_results$seurat_clusters, levels=cluster_sort)
plotDAbeeswarm(WorkerMale_results, group.by = "seurat_clusters") + ggsave("WorkerMale.DAbeeswarm.pdf", height=15, width=5)
write.table(WorkerMale_results, "WorkerMale.rep.res.txt", sep="\t", quote=F)

GyneMale_results <- testNhoods(ant.milo, design = ~ 0 + caste, design.df = ant.design, reduced.dim="PCA", model.contrasts = "casteGyne - casteMale")
GyneMale_results <- annotateNhoods(ant.milo, GyneMale_results, coldata_col = "seurat_clusters")
GyneMale_results$seurat_clusters <- ifelse(GyneMale_results$seurat_clusters_fraction < 0.7, "Mixed", GyneMale_results$seurat_clusters)
GyneMale_results$seurat_clusters <- factor(GyneMale_results$seurat_clusters, levels=cluster_sort)
plotDAbeeswarm(GyneMale_results, group.by = "seurat_clusters") + ggsave("GyneMale.DAbeeswarm.pdf", height=15, width=5)
write.table(GyneMale_results, "GyneMale.rep.res.txt", sep="\t", quote=F)
