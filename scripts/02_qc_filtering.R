# Loading libraries
library(Seurat)
library(tidyverse)

# read RDS file
cd4_temra <- readRDS("data/processed/cd4_temra_loaded.rds")

# -----------------
# 1. Quality control and filtering

#calculating mitochondrial gene %
cd4_temra[["mt_percent"]] <- PercentageFeatureSet(cd4_temra, pattern = "^MT-")

# to check distribution of features/counts/mt_percent
qc_violin <- VlnPlot(cd4_temra, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"), 
        group.by = "subject",
        ncol = 3)

ggsave(filename = "results/figures/qc_violin_before.png",
       plot = qc_violin)

# scatter plot to see relationship between nCount_RNA and nFeature_RNA
qc_scatter_before <- FeatureScatter(cd4_temra, feature1 = "nCount_RNA", 
                                    feature2 = "nFeature_RNA") +
                                    geom_smooth(method = "lm")

ggsave(filename = "results/figures/before_qc_scatter.png",
       plot = qc_scatter_before)

# checking summary statistics
qc <- cd4_temra@meta.data

# Overall summary
qc_features <- c("nFeature_RNA", "nCount_RNA", "mt_percent")
qc_summary <- sapply(qc[, qc_features], summary)

write.csv(qc_summary, "results/tables/qc_summary.csv")

# Distribution quantiles
qc_quantiles <- sapply(qc[, qc_features], quantile,
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
)

write.csv(qc_quantiles, "results/tables/qc_quantiles.csv")

# Filter
cd4_temra_filtered <- subset(cd4_temra, subset = nFeature_RNA > 300 &
                               nFeature_RNA < 2500 &
                               mt_percent < 10)

# checking summary statistics after filtering
qc_filtered <- cd4_temra_filtered@meta.data

qc_summary_filtered <- sapply(qc_filtered[, qc_features], summary)

write.csv(qc_summary_filtered, "results/tables/qc_summary_filtered.csv")

# checking scatter plot after filtering
qc_scatter_after <- FeatureScatter(cd4_temra_filtered, feature1 = "nCount_RNA", 
                                   feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

ggsave(filename = "results/figures/after_qc_scatter.png",
       plot = qc_scatter_after)

# saving current seurat object
saveRDS(cd4_temra_filtered,
        "data/processed/cd4_temra_filtered.rds")

