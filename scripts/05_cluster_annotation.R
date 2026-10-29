# load library
library(Seurat)
library(dplyr)
library(ggplot2)

#read data
cd4_temra_integrated <- readRDS("data/processed/cd4_temra_integrated.rds")

#finding clusters
Idents(cd4_temra_integrated) <- "integrated_cluster"
all_markers <- FindAllMarkers(object = cd4_temra_integrated,
                          logfc.threshold = 0.25,
                          min.pct = 0.25,
                          only.pos = TRUE)

top10_markers <- all_markers |>
  filter(p_val_adj < 0.05, avg_log2FC > 0) |>
  group_by(cluster) |>
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) |>
  ungroup()

top10_markers |>
  select(
    cluster, gene, avg_log2FC,
    pct.1, pct.2, p_val_adj
  ) |>
  print(n = Inf)

cluster_annotations <- c(
  "0" = "Cytotoxic-associated",
  "1" = "Cytotoxic-associated",
  "2" = "Cytotoxic-associated",
  "3" = "Memory-associated",
  "4" = "Cytotoxic-associated",
  "5" = "Cytotoxic-associated",
  "6" = "Unresolved",
  "7" = "MAIT-like",
  "8" = "Cycling"
)

cd4_temra_integrated$broad_state <- unname(
  cluster_annotations[
    as.character(cd4_temra_integrated$integrated_cluster)
  ]
)

DimPlot(
  cd4_temra_integrated,
  reduction = "umap.harmony",
  group.by = "broad_state",
  label = TRUE,
  repel = TRUE
)

state_composition <- cd4_temra_integrated@meta.data |>
  count(subject, broad_state) |>
  group_by(subject) |>
  mutate(percentage = 100 * n / sum(n)) |>
  ungroup()

ggplot(
  state_composition,
  aes(x = subject, y = percentage, fill = broad_state)
) +
  geom_col(width = 0.7) +
  labs(
    x = "Subject",
    y = "Percentage of cells",
    fill = "Transcriptional state",
    title = "CD4+ TEMRA cell-state composition"
  ) +
  theme_classic()








