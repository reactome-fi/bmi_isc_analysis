import scanpy as sc
import os
print(os.getcwd())

out_dir = "bmi_final/"

seurat_file = '../resources/seurat_aHom_aIrr.h5ad'
adata12_file = '../resources/adata_12_frozen.h5ad'
adata17_file = '../resources/adata_17_frozen.h5ad'

adata12 = sc.read_h5ad(adata12_file)
adata17 = sc.read_h5ad(adata17_file)
seurat = sc.read_h5ad(seurat_file)

# Get 12.5 stem markers
sc.tl.rank_genes_groups(adata12, groupby='leiden', groups=['0'])
adata12_dge = sc.get.rank_genes_groups_df(adata12, group=None)
adata12_dge.to_csv(out_dir + 'e12-grp0-markers.csv')

# Get 17.5 stem markers
adata17.obs.leiden.cat.add_categories('0 & 4', inplace=True)
mask = adata17.obs.leiden.isin(['0','4'])
mask.value_counts()
adata17.obs.loc[mask, 'leiden'] = '0 & 4'
sc.tl.rank_genes_groups(adata17, groupby='leiden', groups=['0 & 4'])
adata17_dge = sc.get.rank_genes_groups_df(adata17, group=None)
adata17_dge.to_csv(out_dir + 'e17-grp0-4-markers.csv')

# Get adult/irr stem markers
adult = seurat[seurat.obs.index.str.startswith("Hom_")].copy()
irr = seurat[seurat.obs.index.str.startswith("Irr_")].copy()

sc.tl.rank_genes_groups(adult, groupby='cell_type_combined', groups = ['Stem H'])
stemh_dge = sc.get.rank_genes_groups_df(adult, group='Stem H')
stemh_dge.to_csv(out_dir + 'stem-h-markers.csv')

sc.tl.rank_genes_groups(irr, groupby='cell_type_combined', groups = ['Stem IR'])
stemirr_dge = sc.get.rank_genes_groups_df(irr, group ='Stem IR')
stemirr_dge.to_csv(out_dir + 'stem-irr-markers.csv')
