# ---- Global Vars and Methods ----

# -- Load Libraries --
library(Seurat)
library(ggplot2)
library(cowplot)

# -- Adult Hom Colors -- 
isc <- readRDS("rds/seurat_aHom.rds")
col_aHom <- rep("#000000", length(unique(isc$cell_type)))
names(col_aHom) <- c("EE-EC", 
                     "EE-I",
                     "EE-D",
                     "EE-K",
                     "Sec Pro",
                     "Tuft",
                     "Pan",
                     "Ent",
                     "Stem")

col_aHom["EE-EC"] <-      "#0055aa"
col_aHom["EE-I"] <-       "#0000ff"
col_aHom["EE-D"] <-       "#5555ff"
col_aHom["EE-K"] <-       "#aaaaff"
col_aHom["Tuft"] <-         "#ff551c"
col_aHom["Sec Pro"] <-      "#646464"
col_aHom["Pan"] <-       "#ffa216"
col_aHom["Ent"] <-   "#006400"
col_aHom["Stem"] <-         "#e30000"

# -- Adult Hom and Irr Colors --
isc <- readRDS("rds/seurat_aHom_aIrr.rds")
col_aHom_aIrr <- rep("#000000", length(unique(isc$cell_type_ee_combined)))
names(col_aHom_aIrr) <- c("EE H", 
                          "EE IR",
                          "Sec Pro H",
                          "Sec Pro IR",
                          "Tuft H",
                          "Tuft IR",
                          "Pan H",
                          "Pan IR",
                          "Ent Hom", 
                          "Ent IR",
                          "Sec Pre IR",
                          "Stem H",
                          "Stem IR")

col_aHom_aIrr["EE H"] <-          "#0055aa"
col_aHom_aIrr["EE IR"] <-          "#71b8ff"
col_aHom_aIrr["Tuft H"] <-        "#ff551c"
col_aHom_aIrr["Tuft IR"] <-        "#ff7f00"
col_aHom_aIrr["Sec Pro H"] <-     "#616161"
col_aHom_aIrr["Sec Pro IR"] <-     "#aaaaaa"
col_aHom_aIrr["Pan H"] <-      "#ffa216"
col_aHom_aIrr["Pan IR"] <-      "#ffd900"
col_aHom_aIrr["Ent H"] <-  "#006400"
col_aHom_aIrr["Ent IR"] <-  "#228b22"
col_aHom_aIrr["Sec Pre IR"] <-    "#800080"
col_aHom_aIrr["Stem H"] <-        "#e30000"
col_aHom_aIrr["Stem IR"] <-        "#ff7184"

# ---- Figure 5 - Annotated Adult Hom UMAP --------
# -- Load Libraries --
library(Seurat)
library(ggplot2)

# -- Load/Prep Seurat Obj -- 
isc <- readRDS("rds/seurat_aHom.rds")
isc <- SetIdent(object = isc, value = isc$cell_type)
group_order <- c("EE-EC", 
                 "EE-I", 
                 "EE-D", 
                 "EE-K",
                 "Tuft",
                 "Sec Pro",
                 "Pan",
                 "Ent",
                 "Stem")
isc@active.ident <- factor(x = isc@active.ident, levels = group_order)

# -- Plotting -- 
tiff("figures/aHom_umap_annotated.tiff", units="in", width=8, height=6, res=1000)
DimPlot(isc, 
        reduction = "umap", 
        cols = unname(col_aHom[group_order]),
        pt.size = 1) + coord_equal()
dev.off()

# ---- Figure 5 - Adult Hom Heatmap --------
# -- Load Libraries --
library(Seurat)
library(ggplot2)

# -- Load/Prep Seurat Obj -- 
isc <- readRDS("rds/seurat_aHom.rds")
isc <- SetIdent(isc, value = isc$cell_type)

# -- Plotting -- 
isc <- SetIdent(object = isc, value = isc$cell_type)
group_order <- c("EE-EC", 
                 "EE-I", 
                 "EE-D", 
                 "EE-K",
                 "Tuft",
                 "Sec Pro",
                 "Pan",
                 "Ent",
                 "Stem")
isc@active.ident <- factor(x = isc@active.ident, levels = group_order)

# -- Genes to use -- 
g_ec <- c("Chga", "Chgb", "Tph1")
g_i <- c("Gcg", "Cck", "Ghrl")
g_d <- c("Sst", "Iapp")
g_k <- c("Gip", "Isl1")
g_tuft <- c("Trpm5", "Sox9", "Dclk1")
g_pan <- c("Defa24", "Lyz1", "Mmp7")
g_ent <- c("Alpi", "Lct", "Fabp2")
g_secprog <- c("Defa17", "Defa-rs1", "Itln1")
g_stem <- c("Mki67", "Olfm4", "Ascl2")
g_all <- c(g_ec, g_i, g_d, g_k, g_tuft, g_secprog, g_pan, g_ent, g_stem)

# -- Plotting -- 
tiff("figures/aHom_heatmap.tiff", units="in", width=10, height=6, res=1000)
DoHeatmap(isc,
          group.colors = unname(col_aHom[group_order]),
          features = g_all, 
          draw.lines = TRUE,
          lines.width = 40,
          size = 3,
          angle = 30)
dev.off()

# -- Cleanup --
rm(isc, group_order, g_ec, g_i, g_d, g_k, g_tuft, 
   g_secprog, g_pan, g_ent, g_stem, g_all)

# ---- Figure 5D - Adult Hom Combined violin ----
# -- Load Libraries --
library(Seurat)
library(ggplot2)
library(cowplot)

# -- Load/Prep Seurat Obj -- 
isc <- readRDS("rds/seurat_aHom.rds")
isc <- SetIdent(object = isc, value = isc$cell_type)
group_order <- c("EE-EC", 
                 "EE-I", 
                 "EE-D", 
                 "EE-K",
                 "Tuft",
                 "Sec Pro",
                 "Pan",
                 "Ent",
                 "Stem")
isc@active.ident <- factor(x = isc@active.ident, levels = group_order)

temp <- read.table("tables/CytoTRACE_plot_table_hom.txt", header = TRUE, stringsAsFactors = FALSE)
cytotrace <- temp$CytoTRACE
names(cytotrace) <- rownames(temp)
cytotrace <- cytotrace[which(grepl("hom", names(cytotrace)))]
isc[['cytotrace']] <- cytotrace

# -- Load genes lists -- 
gene_lists <- list(
  deg_12 = read.csv("tables/de_genes_12.csv", header = TRUE, stringsAsFactors = FALSE)$X,
  decrg_12 = read.csv("tables/de_crg_12.csv", header = TRUE, stringsAsFactors = FALSE)$X
)
gene_lists <- lapply(gene_lists, 
                     function(x) x[which(x %in% rownames(isc))])

gene_lists[[1]] <- gene_lists[[1]][1:50]
gene_lists[[2]] <- gene_lists[[2]][1:50]
write.csv(gene_lists[[1]], "tables/table2_devList.csv")
write.csv(gene_lists[[2]], "tables/table2_crgList.csv")

# -- Score Adult Hom with each gene list -- 
isc <- AddModuleScore(isc, features = gene_lists[1], name = "12_DEG_Score")
isc <- AddModuleScore(isc, features = gene_lists[2], name = "12_CRG_Score")

# -- Generate label portion of chart -- 
names <- ggplot() + 
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 9)) +
  annotate("segment",
           x = rep(0, 10),
           xend = rep(10, 10),
           y = c(0:9),
           yend = c(0:9),
           colour = 'black',
           size = 0.75) +
  annotate("text", 
           x = rep(10, 9),
           y = c(0.5:8.5),
           label = group_order,
           hjust = 1,
           family="Arial") +
  theme(plot.margin = unit(c(9, -10, -14, 0), "pt"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())

# -- Define formatted vln plot function -- 
isc_vln <- function(isc, feature, title, group_order) {
  graph <- VlnPlot(isc, 
                   features = feature, 
                   cols = unname(col_aHom[group_order]),
                   pt.size = 0) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(size = 10, family = "Arial"),
          plot.margin = unit(c(10, 1, 1, 1), "pt")) +
    ggtitle(title) +
    NoLegend()
  
  return(graph)
}

# -- Generate charts -- 
p_1 <- isc_vln(isc, "Mki67", "Mki67", group_order)
p_2 <- isc_vln(isc, "Cdk1", "Cdk1", group_order)
p_3 <- isc_vln(isc, "Lgr5", "Lgr5", group_order)
p_4 <- isc_vln(isc, "Olfm4", "Olfm4", group_order)
p_5 <- isc_vln(isc, "Ascl2", "Ascl2", group_order)
p_6 <- isc_vln(isc, "Chga", "ChgA", group_order)
p_7 <- isc_vln(isc, "Chgb", "ChgB", group_order)
p_8 <- isc_vln(isc, "X12_DEG_Score1", "Dev Score", group_order)
p_9 <- isc_vln(isc, "X12_CRG_Score1", "CRG Score", group_order)
p_10 <- isc_vln(isc, "cytotrace", "Cytotrace", group_order)

# -- Combine charts -- 
tiff("figures/aHom_vln.tiff", units="in", width=11, height=5, res=1000)
plot_grid(names, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10,
          align = "v",
          axis = "b",
          ncol = 11, 
          nrow = 1,
          rel_widths = c(0.8,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
dev.off()

# -- Cleanup --
rm(isc, temp, group_order, names, isc_vln, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, p_10)

# ---- Figure 5E - CytoTRACE ----
# -- Load Libraries --
library(extrafont)
library(scales)
library(RColorBrewer)
library(stringr)

temp <- read.table("../resources/CytoTRACE_results_hom_e17_e12.txt", header = TRUE, stringsAsFactors = FALSE, sep='\t')
temp <- data.frame(
  CytoTRACE = temp$CytoTRACE,
  sample = rownames(temp),
  stringsAsFactors = FALSE
)
temp$sample <- substr(temp$sample, start=20, stop=25)
temp$sample[which(temp$sample == "a.hom")] <- "Adult"
temp$sample <- factor(x = temp$sample, levels = c("E12.5", "E17.5", "Adult"))

# -- Combine charts -- 
pdf("figures/5E_cytotrace.pdf", width=8, height=5)
ggplot(temp, aes(x = CytoTRACE, fill=sample)) +
  geom_density(alpha = 0.4) + 
  scale_fill_manual(values=c("#aba9aa", "#28c941", "#57aeff")) + 
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  theme_classic() + 
  ylab("Density") +
  theme(axis.text = element_text(
    size = 8,
    colour = "black"), 
    axis.title = element_text(
      size = 8,
      colour = "black"), 
    legend.title = element_blank(),
    legend.text = element_text(
      size = 8,
      colour = "black")
  )
dev.off()

# -- Cleanup --
rm(temp)

# Figure 6H - Venn All
library(VennDiagram)
library(RColorBrewer)

outdir <- '../resources/'
markers12 <- read.csv(paste0(outdir,'e12-grp0-markers.csv'))
markers17 <- read.csv(paste0(outdir,'e17-grp0-4-markers.csv'))
markersh <- read.csv(paste0(outdir,'stem-h-markers.csv'))
markersirr <- read.csv(paste0(outdir,'stem-irr-markers.csv'))

e12 <- subset(markers12, logfoldchanges > 2 & pvals_adj < 0.01)$names
e17 <- subset(markers17, logfoldchanges > 2 & pvals_adj < 0.01)$names
ehom <- subset(markersh, logfoldchanges > 2 & pvals_adj < 0.01)$names
eirr <- subset(markersirr, logfoldchanges > 2 & pvals_adj < 0.01)$names

venn.diagram(
  x = list(e12, e17, ehom, eirr),
  category.names = c("E12.5" , "E17.5" , "Adult", "Adult IR"),
  fill = c("#dedede","#157023","#ff0000","#ff7184"),
  alpha=0.75,
  filename = 'logfc-fdr-overlap-color-mid-opacity.png',
  output=TRUE
)



# ---- Figure S7 A - Adult Hom vs Irr Sec Int Vln ---
# -- Load/Prep Seurat Obj -- 
isc <- readRDS("rds/seurat_aHom_aIrr.rds")
isc <- SetIdent(isc, value = isc$cell_type_ee_combined)
group_order <- c("EE H", 
                 "EE IR",
                 "Tuft H",
                 "Tuft IR",
                 "Sec Pro H",
                 "Sec Pro IR",
                 "Pan H", 
                 "Pan IR",
                 "Ent H",
                 "Ent IR",
                 "Sec Pre IR",
                 "Stem H",
                 "Stem IR")
isc@active.ident <- factor(x = isc@active.ident, levels = group_order)
tiff("figures/aHom_aIrr_qcVln.tiff", units="in", width=12, height=5, res=1000)
a <- VlnPlot(isc, features = "nFeature_RNA", cols = unname(col_aHom_aIrr[group_order]), pt.size = 0.25) + 
  NoLegend() +
  ggtitle("Number of Unique Genes per Cell") + 
  theme(axis.title.x = element_blank())
b <- VlnPlot(isc, features = "nCount_RNA", cols = unname(col_aHom_aIrr[group_order]), pt.size = 0.25) + 
  NoLegend()  +
  ggtitle("Number of Molecules per Cell") + 
  theme(axis.title.x = element_blank())
plot_grid(a, b, ncol = 2, nrow = 1)
dev.off()

# ---- Figure 6 - Annotated Adult Hom vs Irr UMAP v2 ----
# -- Load/Prep Seurat Obj -- 
isc <- readRDS("rds/seurat_aHom_aIrr.rds")
isc <- SetIdent(isc, value = isc$cell_type_combined)
group_order <- c("Other H", 
                 "Other IR",
                 "Sec Pre IR",
                 "Stem H",
                 "Stem IR")
isc@active.ident <- factor(x = isc@active.ident, levels = group_order)

col_aHom_aIrr_v2 <- rep("#000000", length(unique(isc$cell_type_combined)))
names(col_aHom_aIrr_v2) <- c("Other H", 
                             "Other IR",
                             "Sec Pre IR",
                             "Stem H",
                             "Stem IR")

col_aHom_aIrr_v2["Other H"] <-          "#616161"
col_aHom_aIrr_v2["Other IR"] <-          "#aaaaaa"
col_aHom_aIrr_v2["Sec Pre IR"] <-    "#800080"
col_aHom_aIrr_v2["Stem H"] <-        "#e30000"
col_aHom_aIrr_v2["Stem IR"] <-        "#ff7184"

# -- Plotting -- 
tiff("figures/aHom_aIrr_umap_annotated.tiff", units="in", width=8, height=6, res=1000)
DimPlot(isc, 
        reduction = "umap", 
        cols = unname(col_aHom_aIrr_v2[group_order]),
        pt.size = 1) + coord_equal()
dev.off()

# ---- Figure S7 - Adult Hom vs Irr Sec Prec UMAP ----
# -- Load Libraries --
#install.packages("viridis")
library(viridis)
# -- Load/Prep Seurat Obj -- 
isc <- readRDS("rds/seurat_aHom_aIrr.rds")
isc <- SetIdent(isc, value = isc$cell_type_combined)
group_order <- c("Other H", 
                 "Other IR",
                 "Sec Pre IR",
                 "Stem H",
                 "Stem IR")
isc@active.ident <- factor(x = isc@active.ident, levels = group_order)

col_aHom_aIrr_v2 <- rep("#000000", length(unique(isc$cell_type_combined)))
names(col_aHom_aIrr_v2) <- c("Other H", 
                             "Other IR",
                             "Sec Pre IR",
                             "Stem H",
                             "Stem IR")

col_aHom_aIrr_v2["Other H"] <-          "#616161"
col_aHom_aIrr_v2["Other IR"] <-          "#aaaaaa"
col_aHom_aIrr_v2["Sec Pre IR"] <-    "#800080"
col_aHom_aIrr_v2["Stem H"] <-        "#e30000"
col_aHom_aIrr_v2["Stem IR"] <-        "#ff7184"

cells <- WhichCells(isc, idents = c("Stem IR", "Sec Pre IR"))

p1 <-  DimPlot(isc, 
               reduction = "umap", 
               cols = unname(col_aHom_aIrr_v2[group_order]),
               pt.size = 0.25) + NoLegend() + coord_equal()

a <- DimPlot(isc, 
             reduction = "umap", 
             cols = unname(col_aHom_aIrr_v2[group_order])[c(3, 5)],
             pt.size = 0.25,
             cells = cells) + NoLegend() +
  coord_equal(xlim=c(-3.25, 0.25), ylim = c(-4.5, -1)) +
  theme(plot.margin = margin(t = 29, r = 0, b = 10, l = -8.5, unit = "pt"),
        axis.title = element_blank(),
        axis.text = element_blank())

myFeatPlot <- function(obj, feat, cells, order) {
  temp <- FeaturePlot(obj, 
                      feature = feat, 
                      pt.size = 0.25,
                      order = order,
                      cells = cells) + 
    scale_colour_gradientn(colors=rev(viridis_pal()(10)),
                           labels=function(x) sprintf("%.1f", x)) +
    coord_equal(xlim=c(-3.25, 0.25), ylim = c(-4.5, -1)) +
    theme(plot.margin = margin(t = 10, 20, 10, 10, unit = "pt"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(hjust = 0.05, size = 14),
          legend.text = element_text(size = 8),
          legend.key.height = unit(0.8, "lines"))
  legend <- get_legend(temp)
  ggdraw(temp + NoLegend()) +
    draw_grob(legend, x = 0.7675, y = -0.17, clip = 'on')
}


b <- myFeatPlot(isc, "Olfm4", cells, TRUE)
c <- myFeatPlot(isc, "Ascl2", cells, TRUE)
d <- myFeatPlot(isc, "Myc", cells, TRUE)
e <- myFeatPlot(isc, "Cdk4", cells, TRUE)
f <- myFeatPlot(isc, "Cdk6", cells, TRUE)

g <- myFeatPlot(isc, "Clu", cells, TRUE)
h <- myFeatPlot(isc, "Anxa1", cells, TRUE)
i <- myFeatPlot(isc, "Mif", cells, TRUE)
j <- myFeatPlot(isc, "Npm1", cells, TRUE)
k <- myFeatPlot(isc, "Cxadr", cells, TRUE)
l <- myFeatPlot(isc, "Cd44", cells, TRUE)

m <- myFeatPlot(isc, "Chga", cells, FALSE)
n <- myFeatPlot(isc, "Cck", cells, FALSE)
o <- myFeatPlot(isc, "Ghrl", cells, FALSE)
p <- myFeatPlot(isc, "Nts", cells, FALSE)
q <- myFeatPlot(isc, "Neurod1", cells, FALSE)
r <- myFeatPlot(isc, "Serpina1c", cells, TRUE)

p2 <- cowplot::plot_grid(
  a, b, c, d, e, f, NULL,
  g, h, i, j, k, l, NULL,
  m, n, o, p, q, r, NULL, 
  nrow = 3,
  ncol = 7,
  rel_widths = rep(c(rep(1, 6), 0.1), 3)
  
)

pdf("figures/Fig6-supC_0-25.pdf", width=18, height=6)
plot_grid(p1, NULL, p2, ncol = 3, nrow = 1, rel_widths = c(1, 0.05, 2.5))
dev.off()


