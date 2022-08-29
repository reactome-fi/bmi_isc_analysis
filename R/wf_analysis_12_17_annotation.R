# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# Load 12.5 and 17.5 Seurat Object
isc <- readRDS("rds/seurat_12_17.rds")
Idents(isc) <- "seurat_clusters"

# Annotate cells based on figure 4 and seurat clusters
temp <- rep(NA, ncol(isc))
names(temp) <- colnames(isc)
temp[WhichCells(isc, idents = c(18, 24, 2))] <- "17.5 Stem"
temp[WhichCells(isc, idents = c(14, 19, 8))] <- "Progenitor"
temp[WhichCells(isc, idents = c(13, 6, 21, 10, 20, 9, 23))] <- "EE Prog"
temp[WhichCells(isc, idents = c(15, 4))] <- "Unk-A"
temp[WhichCells(isc, idents = c(16, 22))] <- "Unk-B"
temp[WhichCells(isc, idents = c(0, 3, 1, 7, 12, 25, 11, 17, 5))] <- "12.5 Stem"
isc[['annotated_clusters']] <- temp
rm(temp)
Idents(isc) <- "annotated_clusters"

# Make and write df with information for table
df <- data.frame(
  "cell" = colnames(isc), 
  "timepoint" = as.character(isc$orig.ident),
  "seurat_clusters" = as.character(as.numeric(isc$seurat_clusters)),
  "annotated_clusters" = as.character(isc$annotated_clusters)
)
write.csv(df, "12_17_cell_annotations.csv", row.names = FALSE)

# Find seurat and annotated cluster DEGs
Idents(isc) <- "seurat_clusters"
temp <- FindAllMarkers(isc, 
                       only.pos = TRUE)
write.csv(temp, "de_genes_12_17_seurat_clusters.csv")

Idents(isc) <- "annotated_clusters"
temp <- FindAllMarkers(isc, 
                       only.pos = TRUE)
write.csv(temp, "de_genes_12_17_annotated_clusters.csv")