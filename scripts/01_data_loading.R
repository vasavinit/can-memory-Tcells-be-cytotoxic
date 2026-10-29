# Loading libraries
library(Seurat)
library(Matrix)

# loading the data
# as it is a .txt file
counts <- read.delim("data/raw/GSE106543_SC_raw_Counts_10x.txt",
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

# converting counts to sparse matrix to save memory
counts_sparse <- Matrix(as.matrix(counts),
                        sparse = TRUE)
rm(counts)

# converting it to a seurat object
cd4_temra <- CreateSeuratObject(counts = counts_sparse,
                                project = "CD4_TEMRA")


# creating subject column as subject info is part of the barcode
cd4_temra$subject <- sub("SCRNA_(Subject[12])_.*",
                         "\\1",
                         colnames(cd4_temra))

# save current seurat object
saveRDS(cd4_temra,
        "data/processed/cd4_temra_loaded.rds")




