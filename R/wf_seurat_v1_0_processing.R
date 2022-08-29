# 

# These workflows are adapted from two vignettes provided by Seurat
# https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# https://satijalab.org/seurat/v3.1/sctransform_vignette.html

# ---- Global Vars and Methods ----

# -- Create directories for output files --
setwd('/path/to/bmi_isc_analysis/R')
dir.create(file.path(getwd(), "rds"), showWarnings = FALSE)

# ---- Adult Hom Workflow ----
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("Hom" = "../raw/Adult_Hom_gfp/filtered_feature_bc_matrix")

isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 3, min.features = 200)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
# VlnPlot(isc, features = "percent.mt")
isc <- subset(isc, subset = percent.mt < 20)

# -- Normalization, Scaling, and Find Variable Features for PCA --
isc <- SCTransform(isc, verbose = FALSE)

# -- Linear Dimensional reduction and PC selection --
# Perform Linear dimensional reduction
isc <- RunPCA(isc, features = VariableFeatures(object = isc))

# -- Graph-based clustering using PCs --
isc <- FindNeighbors(isc, dims = 1:11)
isc <- FindClusters(isc, resolution = 1.2)

# -- Non linear dimensional reduction -- 
isc <- RunUMAP(isc, dims = 1:11)

# -- Add Determined groupings to dataset -- 
temp <- rep(NA, ncol(isc))
names(temp) <- colnames(isc)
temp[WhichCells(isc, idents = c(17, 9, 10, 8, 4, 17))] <- "Tuft"
temp[WhichCells(isc, idents = c(13))] <- "EE-D"
temp[WhichCells(isc, idents = c(21, 1, 7))] <- "EE-EC"
temp[WhichCells(isc, idents = c(0, 2, 12))] <- "EE-K"
temp[WhichCells(isc, idents = c(6, 19))] <- "EE-I"
temp[WhichCells(isc, idents = c(11))] <- "Pan"
temp[WhichCells(isc, idents = c(14))] <- "Ent"
temp[WhichCells(isc, idents = c(16))] <- "Stem"
temp[WhichCells(isc, idents = c(18, 3, 20, 5, 15))] <- "Sec Pro"
isc[['cell_type']] <- temp
rm(temp)

# -- Save Seurat Obj for results generation -- 
saveRDS(isc, file = "rds/seurat_aHom.rds")

# ---- Adult Hom and Irr Workflow ----
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("Hom" = "../raw/Adult_Hom_gfp/filtered_feature_bc_matrix",
             "Irr" = "../raw/Adult_Irr_gfp/filtered_feature_bc_matrix")
isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 3, min.features = 200)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
isc <- subset(isc, subset = percent.mt < 20)

# -- Normalization, Scaling, and Find Variable Features for PCA --
isc <- SCTransform(isc, verbose = FALSE)

# -- Linear Dimensional reduction and PC selection --
# Perform Linear dimensional reduction
isc <- RunPCA(isc, features = VariableFeatures(object = isc))

# -- Graph-based clustering using PCs --
isc <- FindNeighbors(isc, dims = 1:11)
isc <- FindClusters(isc, resolution = 1.8)

# -- Non linear dimensional reduction -- 
isc <- RunUMAP(isc, dims = 1:11)

# -- Add determined groupings to dataset --
isc <- SetIdent(object = isc, value = isc$orig.ident)
hom <- WhichCells(isc, idents = "Hom")
irr <- WhichCells(isc, idents = "Irr")

# 1.8
isc <- SetIdent(object = isc, value = isc$seurat_clusters)
ee_ec <- WhichCells(isc, idents = c(8, 26, 0, 5, 35, 14, 22))
ee_i <- WhichCells(isc, idents = c(6, 17, 25, 29))
ee_d <- WhichCells(isc, idents = c(12, 27))
ee_k <- WhichCells(isc, idents = c(1, 11, 21, 15, 9, 36))
ee <- c(ee_ec, ee_i, ee_d, ee_k)
tuft <- WhichCells(isc, idents = c(23, 7, 4, 2, 24))
sec_int <- WhichCells(isc, idents = c(3, 10, 13, 33, 16))
sec_prec <- WhichCells(isc, idents = c(30, 19))
pan <- WhichCells(isc, idents = c(32, 18, 31))
ent <- WhichCells(isc, idents = c(34))
stem <- WhichCells(isc, idents = c(20, 28))

temp <- rep(NA, ncol(isc))
names(temp) <- colnames(isc)
temp[hom] <- "Other H"
temp[irr] <- "Other IR"
temp[intersect(hom, stem)] <- "Stem H"
temp[intersect(irr, stem)] <- "Stem IR"
temp[sec_prec] <- "Sec Pre IR"
isc[['cell_type_combined']] <- temp

temp[intersect(hom, ee)] <- "EE H"
temp[intersect(irr, ee)] <- "EE IR"
temp[intersect(hom, tuft)] <- "Tuft H"
temp[intersect(irr, tuft)] <- "Tuft IR"
temp[intersect(hom, sec_int)] <- "Sec Pro H"
temp[intersect(irr, sec_int)] <- "Sec Pro IR"
temp[intersect(hom, pan)] <- "Pan H"
temp[intersect(irr, pan)] <- "Pan IR"
temp[intersect(hom, ent)] <- "Ent H"
temp[intersect(irr, ent)] <- "Ent IR"
temp[intersect(hom, stem)] <- "Stem H"
temp[intersect(irr, stem)] <- "Stem IR"
temp[sec_prec] <- "Sec Pre IR"
isc[['cell_type_ee_combined']] <- temp

temp[intersect(hom, ee_ec)] <- "EE-EC H"
temp[intersect(irr, ee_ec)] <- "EE-EC IR"
temp[intersect(hom, ee_i)] <- "EE-I H"
temp[intersect(irr, ee_i)] <- "EE-I IR"
temp[intersect(hom, ee_d)] <- "EE-D H"
temp[intersect(irr, ee_d)] <- "EE-D IR"
temp[intersect(hom, ee_k)] <- "EE-K H"
temp[intersect(irr, ee_k)] <- "EE-K IR"
isc[['cell_type_ee_separate']] <- temp

# -- Save Seurat Obj for results generation -- 
saveRDS(isc, file = "rds/seurat_aHom_aIrr.rds")

# ---- E12.5 and Adult Hom Workflow ----
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("12_5" = "../raw/12_5_gfp/filtered_feature_bc_matrix",
             "Hom" = "../raw/Adult_Hom_gfp/filtered_feature_bc_matrix")
isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 3, min.features = 200)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
isc <- subset(isc, subset = percent.mt < 20)

# -- Normalization, Scaling, and Find Variable Features for PCA --
isc <- SCTransform(isc, verbose = FALSE)

# -- Save Seurat Obj for results generation -- 
saveRDS(isc, file = "rds/seurat_12_aHom.rds")

# ---- E17.5 and Adult Hom Workflow ----
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("17_5" = "../raw/17_5_gfp/filtered_feature_bc_matrix",
             "Hom" = "../raw/Adult_Hom_gfp/filtered_feature_bc_matrix")
isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 3, min.features = 200)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
isc <- subset(isc, subset = percent.mt < 20)

# -- Normalization, Scaling, and Find Variable Features for PCA --
isc <- SCTransform(isc, verbose = FALSE)

# -- Linear Dimensional reduction and PC selection --
# Perform Linear dimensional reduction
isc <- RunPCA(isc, features = VariableFeatures(object = isc))

# -- Save Seurat Obj for results generation -- 
saveRDS(isc, file = "rds/seurat_17_aHom.rds")

# ---- E12.5 and E17.5 Workflow ----
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("12_5" = "../raw/12_5_gfp/filtered_feature_bc_matrix",
             "17_5" = "../raw/17_5_gfp/filtered_feature_bc_matrix")
isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 3, min.features = 200)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
isc <- subset(isc, subset = percent.mt < 20)

# -- Normalization, Scaling, and Find Variable Features for PCA --
isc <- SCTransform(isc, verbose = FALSE)

# -- Linear Dimensional reduction and PC selection --
# Perform Linear dimensional reduction
isc <- RunPCA(isc, features = VariableFeatures(object = isc))

# -- Graph-based clustering using PCs --
isc <- FindNeighbors(isc, dims = 1:12)
isc <- FindClusters(isc, resolution = 1.2)

# -- Non linear dimensional reduction -- 
isc <- RunUMAP(isc, dims = 1:12)

# -- Save Seurat Obj for results generation -- 
saveRDS(isc, file = "rds/seurat_12_17.rds")

# ---- Yan All Cells Dataset ---- 
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("Bmi1_1" = "../raw/Yan/Bmi1_1",
             "Bmi1_2" = "../raw/Yan/Bmi1_2",
             "Lgr5Neg_1" = "../raw/Yan/Lgr5_Neg_1",
             "Lgr5Neg_2" = "../raw/Yan/Lgr5_Neg_2",
             "Lgr5Pos_1" = "../raw/Yan/Lgr5_Pos_1",
             "Lgr5Pos_2" = "../raw/Yan/Lgr5_Pos_2",
             "Prox1_1" = "../raw/Yan/Prox1_1",
             "Prox1_2" = "../raw/Yan/Prox1_2")

isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 2, min.features = 100)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
# VlnPlot(isc, features = "percent.mt")
isc <- subset(isc, subset = nFeature_RNA >= 100)
isc <- subset(isc, subset = nFeature_RNA <= 5000)
isc <- subset(isc, subset = percent.mt <= 10)

isc <- NormalizeData(isc, normalization.method = "LogNormalize", scale.factor = 10000)
isc <- FindVariableFeatures(isc, selection.method = "disp", mean.cutoff = c(-0.5, 0.5), dispersion.cutoff = c(-0.5, 0.5), nfeatures = 2419)
all.genes <- rownames(isc)
isc <- ScaleData(isc, features = all.genes)
isc <- RunPCA(isc, features = VariableFeatures(object =isc))
isc <- RunTSNE(isc, dims = 1:20)
# -- Graph-based clustering using PCs --
isc <- FindNeighbors(isc, dims = 1:20, k.param = 100)
isc <- FindClusters(isc, resolution = 0.4)

# -- Save File -- 
saveRDS(isc, "rds/yan_all.rds")

# ---- Yan Bmi1-GFP Cells ----
# -- Load Libraries -- 
library(Seurat)
library(sctransform)
library(ggplot2)

# -- Import data and set up object --
dirList <- c("Bmi1_1" = "../raw/Yan/Bmi1_1",
             "Bmi1_2" = "../raw/Yan/Bmi1_2")
isc.data <- Read10X(data.dir = dirList)
isc <- CreateSeuratObject(counts = isc.data, project = "test", min.cells = 2, min.features = 100)

# -- Processing and QC --
isc[["percent.mt"]] <- PercentageFeatureSet(isc, pattern = "^mt-")
isc <- subset(isc, subset = nFeature_RNA >= 100)
isc <- subset(isc, subset = nFeature_RNA <= 5000)
isc <- subset(isc, subset = percent.mt <= 10)

# -- Normalization, Scaling, and Find Variable Features for PCA --
isc <- SCTransform(isc, verbose = FALSE)

# -- Linear Dimensional reduction and PC selection --
# Perform Linear dimensional reduction
isc <- RunPCA(isc, features = VariableFeatures(object = isc))

# -- Graph-based clustering using PCs --
isc <- FindNeighbors(isc, dims = 1:9)
isc <- FindClusters(isc, resolution = 0.4)

# -- Non linear dimensional reduction -- 
isc <- RunUMAP(isc, dims = 1:9)

# -- Save File -- 
saveRDS(isc, "rds/yan_bmi1.rds")

# ---- Yan Bmi1-GFP Integration ----
# -- Load Libraries --
library(Seurat)
library(ggplot2)
library(sctransform)
options(future.globals.maxSize = 4000 * 1024^2)

# -- Integrate Data --
yan <- readRDS("rds/yan_bmi1.rds")
isc <- readRDS("rds/seurat_aHom.rds")

isc.list <- list(yan, isc)
for (i in 1:length(isc.list)) {
  isc.list[[i]] <- SCTransform(isc.list[[i]], verbose = FALSE)
}
isc.features <- SelectIntegrationFeatures(object.list = isc.list, nfeatures = 3000)
isc.list <- PrepSCTIntegration(object.list = isc.list, anchor.features = isc.features, 
                               verbose = FALSE)
isc.anchors <- FindIntegrationAnchors(object.list = isc.list, normalization.method = "SCT", 
                                      anchor.features = isc.features, verbose = FALSE)
isc.integrated <- IntegrateData(anchorset = isc.anchors, normalization.method = "SCT", 
                                verbose = FALSE)
isc.integrated <- RunPCA(isc.integrated, verbose = FALSE)
isc.integrated <- RunUMAP(isc.integrated, dims = 1:9)

# -- Add Cluster Information -- 
temp1 <- as.numeric(as.character(yan$seurat_clusters))
names(temp1) <- names(yan$seurat_clusters)
temp1[which(temp1 == 5)] <- "Bmi1-2 Stem"
temp1[which(temp1 != "Bmi1-2 Stem")] <- "Bmi1-2 Other"
temp2 <- as.character(isc$cell_type)
names(temp2) <- names(isc$cell_type)
temp2[which(temp2 == "Stem")] <- "Bmi1-1 Stem"
temp2[which(temp2 != "Bmi1-1 Stem")] <- "Bmi1-1 Other"

isc.integrated[['idents']] <- c(temp1, temp2)
isc.integrated <- SetIdent(isc.integrated, value =isc.integrated$idents)

# -- Save Dataset -- 
saveRDS(isc.integrated, "rds/seurat_aHom_yan.rds")

# ---- Li CRC Dataset ----
# -- Load Libraries -- 
library(dplyr)
library(Seurat)
library(umap)
library(ggplot2)
library(cowplot)
library(plyr)

# -- Process CV into matrix -- 
data_nm <- read.csv("../raw/Li_CRC/GSE81861_CRC_NM_epithelial_cells_COUNT.csv", header = FALSE,stringsAsFactors = FALSE)
crc.data <- data_nm[2:nrow(data_nm),2:ncol(data_nm)]
crc.data <- as.matrix(crc.data)
storage.mode(crc.data) <- "numeric"
full_rownames <- data_nm[2:nrow(data_nm),1]
full_colnames <- as.list(data_nm[1,2:ncol(data_nm)])
gene_locs <- character(0)
gene_names <- character(0)
gene_ids <- character(0)
cell_ids <- character(0)
cell_type <- character(0)
cell_sample <- character(0)
cell_group <- character(0)

for(obj in full_rownames) {
  temp <- strsplit(obj, "_")
  gene_locs <- c(gene_locs, temp[[1]][1])
  gene_names <- c(gene_names, temp[[1]][2])
  gene_ids <- c(gene_ids, temp[[1]][3])
}

for(obj in full_colnames) {
  temp <- strsplit(obj, "__")
  cell_ids <- c(cell_ids, temp[[1]][1])
  cell_type <- c(cell_type, temp[[1]][2])
  cell_sample <- c(cell_sample, "Normal")
  cell_group <- c(cell_group, paste0("Normal ", temp[[1]][2]))
}

data_tm <- read.csv("../raw/Li_CRC/GSE81861_CRC_tumor_epithelial_cells_COUNT.csv", header = FALSE,stringsAsFactors = FALSE)
crc.data_tm  <- data_tm[2:nrow(data_tm),2:ncol(data_tm)]
crc.data_tm  <- as.matrix(crc.data_tm)
storage.mode(crc.data_tm) <- "numeric"
crc.data <- cbind2(crc.data, crc.data_tm)
full_colnames <- as.list(data_tm[1,2:ncol(data_tm)])

for(obj in full_colnames) {
  temp <- strsplit(obj, "__")
  cell_ids <- c(cell_ids, temp[[1]][1])
  cell_type <- c(cell_type, temp[[1]][2])
  cell_sample <- c(cell_sample, "Tumor")
  cell_group <- c(cell_group, paste0("Tumor ", temp[[1]][2]))
}
names(cell_type) <- cell_ids
names(cell_sample) <- cell_ids
names(cell_group) <- cell_ids
rownames(crc.data) <- gene_names
colnames(crc.data) <- cell_ids
rm(data_nm, data_tm, full_rownames, full_colnames, crc.data_tm)

# -- Seurat WF --
# Load data and initialize Seurat object with non-normalized data
crc <- CreateSeuratObject(counts = crc.data, project = "test", min.cells = 3, min.features = 200)
crc[['type']] <- cell_type
crc[['sample']] <- cell_sample
crc[['group']] <- cell_group

# -- Normalization, Scaling, and Find Variable Features for PCA --
crc <- SCTransform(crc, verbose = FALSE)

# -- Linear Dimensional reduction and PC selection --
# Perform Linear dimensional reduction
crc <- RunPCA(crc, features = VariableFeatures(object = crc))

# -- Graph-based clustering using PCs --
crc <- FindNeighbors(crc, dims = 1:13)
crc <- FindClusters(crc, resolution = 1.2)

# -- Non linear dimensional reduction -- 
crc <- RunUMAP(crc, dims = 1:13)

# -- Add Groupings -- 
bmi1_status <- rep("BMI1-", ncol(crc))
names(bmi1_status) <- colnames(crc)
bmi1_status[names(which(crc@assays$RNA@counts["BMI1",] > 0))] <-"BMI1+"
crc[['bmi1_status']] <- bmi1_status
lgr5_status <- rep("LGR5-", ncol(crc))
names(lgr5_status) <- colnames(crc)
lgr5_status[names(which(crc@assays$RNA@counts["LGR5",] > 0))] <-"LGR5+"
crc[['lgr5_status']] <- lgr5_status
ascl2_status <- rep("ASCL2-", ncol(crc))
names(ascl2_status) <- colnames(crc)
ascl2_status[names(which(crc@assays$RNA@counts["ASCL2",] > 0))] <-"ASCL2+"
crc[['ascl2_status']] <- ascl2_status
crc$group[which(crc$group == "Normal stemTA")] <- "Normal Stem/TA"
crc$group[which(crc$group == "Normal GobletA")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal GobletB")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal GobletC")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal EnterocyteType1A")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal EnterocyteType1B")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal EnterocyteType2")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal NonStem1")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal NonStem2")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Normal Unknown")] <- "Normal Non-Stem/TA"
crc$group[which(crc$group == "Tumor stemTA")] <- "Tumor Stem/TA"
crc$group[which(crc$group == "Tumor EnterocyteLike")] <- "Tumor Non-Stem/TA"
crc$group[which(crc$group == "Tumor GobletLike")] <- "Tumor Non-Stem/TA"

# -- Save Seurat Obj for results generation -- 
saveRDS(crc, file = "rds/seurat_li_crc.rds")

