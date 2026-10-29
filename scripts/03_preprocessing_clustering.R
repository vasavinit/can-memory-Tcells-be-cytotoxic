# Loading libraries
library(Seurat)
library(tidyverse)

# read filtered RDS File
cd4_temra_filtered <- readRDS("data/processed/cd4_temra_filtered.rds")

# checking for batch effects if any due to 2 different subjects
cd4_temra_filtered <- NormalizeData(object = cd4_temra_filtered)

#scaling data
cd4_temra_filtered <- ScaleData(object = cd4_temra_filtered,
                                features = VariableFeatures(cd4_temra_filtered))

#finding most variable features
cd4_temra_filtered <- FindVariableFeatures(object = cd4_temra_filtered)
plot_variables <- VariableFeaturePlot(cd4_temra_filtered)
top10 <- head(VariableFeatures(cd4_temra_filtered), 10)
variablefeat_top10 <- LabelPoints(plot = plot_variables, points = top10, repel = TRUE)
ggsave("results/figures/variable_top10.png",
       plot = variablefeat_top10)

#pca
cd4_temra_filtered <- RunPCA(object = cd4_temra_filtered)
elbow_plot <- ElbowPlot(cd4_temra_filtered)
ggsave("results/figures/elbow_plot.png",
       plot = elbow_plot)

#neighbours
cd4_temra_filtered <- FindNeighbors(object = cd4_temra_filtered, dims = 1:15)

#clusters
cd4_temra_filtered <- FindClusters(object = cd4_temra_filtered)

# umap
cd4_temra_filtered <- RunUMAP(object = cd4_temra_filtered, dims = 1:15)
umap_plot <- DimPlot(cd4_temra_filtered, reduction = "umap", group.by = "subject")
ggsave(filename = "results/figures/unintegrated_umap_plot.png",
        plot = umap_plot)

# checking cluster composition by subject
cluster_subject <- table(Cluster = cd4_temra_filtered$seurat_clusters,
                         Subject = cd4_temra_filtered$subject)

# Percentage of each cluster belonging to each subject
cluster_subject_pct <- round(100 * prop.table(cluster_subject, margin = 1),2)

write.csv(cluster_subject_pct, "results/tables/unintegrated_cluster_subject_pct.csv")

saveRDS(cd4_temra_filtered,
        "data/processed/cd4_temra_clustered.rds")






