# library loading
library(Seurat)
library(harmony)

cd4_temra_integrated <- readRDS("data/processed/cd4_temra_clustered.rds")

# keeping original unintegrated clusters
cd4_temra_integrated$unintegrated_cluster <- 
  as.character(cd4_temra_integrated$seurat_clusters)

# integrating
set.seed(42)

cd4_temra_integrated <- harmony::RunHarmony(object = cd4_temra_integrated,
                                            group.by.vars = "subject",
                                            reduction.use = "pca",
                                            dims.use = 1:15,
                                            reduction.save = "harmony",
                                            project.dim = FALSE)


# cluster again
cd4_temra_integrated <- FindNeighbors(cd4_temra_integrated,
                                      reduction = "harmony",
                                      dims = 1:15)
cd4_temra_integrated <- FindClusters(cd4_temra_integrated)

cd4_temra_integrated$integrated_cluster <-
  as.character(Idents(cd4_temra_integrated))

cd4_temra_integrated <- RunUMAP(cd4_temra_integrated,
                                reduction = "harmony",
                                dims = 1:15,
                                reduction.name = "umap.harmony",
                                seed.use = 42)


integrated_umap <- DimPlot(cd4_temra_integrated, reduction = "umap.harmony",
       group.by = "subject")

ggsave(filename = "results/figures/integrated_umap.png",
       plot = integrated_umap)

saveRDS(cd4_temra_integrated,
        file = "data/processed/cd4_temra_integrated.rds")
