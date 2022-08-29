import os
from scpy4reactome import pathway_analyzer as pa

e17_raw_counts_dir = '../raw/17_5_gfp/filtered_feature_bc_matrix'
e12_raw_counts_dir = '../raw/12_5_gfp/filtered_feature_bc_matrix'
reactome_gmt_file = "../resources/MouseReactomePathways_Rel_75_122220.gmt"
out_dir = "bmi_final/"

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
os.getcwd()

exec(open('PlotGenerators.py').read())

### Process E17.5 data ###
adata_17 = open_data(e17_raw_counts_dir)
run_umap(adata_17, need_run_variable=True)
adata_17 = adata_17[~adata_17.obs.leiden.isin(['0','14','15','20','22','23'])].copy()
run_umap(adata_17, need_run_variable=True)
adata_17 = adata_17[~adata_17.obs.leiden.isin(['16'])].copy()
run_umap(adata_17, need_run_variable=True)
sc.tl.rank_genes_groups(adata_17, 'leiden', only_positive=True, n_genes=200, use_raw=False)
file = out_dir + 'adata_17_cluster_gene_markers.csv'
df = sc.get.rank_genes_groups_df(adata_17, group=None)
df.to_csv(file)

# Pathway Analysis
anova_results_file = out_dir + 'adata_17_anova_results.csv' 
pathway_results_file = out_dir + 'adata_17_pathway_results.csv'
pa_result = pa.reactome_aucell(adata_17, reactome_gmt_file)
anova_results = pa_result.perform_pathway_anova()
anova_results.to_csv(anova_results_file)
pa_result.adata.obsm['X_aucell'].to_csv(pathway_results_file)
del(pa_result)

# PAGA/CytoTRACE
# output for CytoTRACE run on web
# https://cytotrace.stanford.edu/
cyto_df_17 = adata_17.copy()
cyto_df_17.X = cyto_df_17.layers['counts']
cyto_df_17 = cyto_df_17.to_df()
cyto_df_17.transpose().to_csv(out_dir + 'adata_17_cyto_input_counts.csv')
del(cyto_df_17)

# results from CytoTRACE run on web - output saved and renamed
cyto_web_results_17 = '../resources/CytoTRACE_results_e17_counts.txt'
df17 = pd.read_csv(cyto_web_results_17, index_col=0, delimiter='\t')
df17.index = df17.index.str.replace('.','-')
adata_17.obs['CytoTRACE'] = df17['CytoTRACE']

# Get root cell
root_cell_17 = df17.CytoTRACE.idxmax()
adata_17.uns['iroot'] = np.flatnonzero(adata_17.obs_names==root_cell_17)[0]
del(df17)

# Get PAGA and DPT scores for 17.5 cells only
run_paga(adata_17)
run_dpt(adata_17, root_cell_17)

# Spot check for review...can be removed...generated in figures script
sc.pl.paga(adata_17,
           threshold=0.09,
           color = ['leiden', 'CytoTRACE', 'dpt_pseudotime'],
           pos = adata_17.uns['paga']['pos'],
           cmap = 'viridis',
           plot=False) # If False, do not create the figure, simply compute the layout.

### Process E12.5 data ###
adata_12 = open_data(e12_raw_counts_dir)
# UMAP generation with PAGA initialization
run_umap_via_paga(adata_12)
# Decrease # of leiden clusters from default
sc.tl.leiden(adata_12, resolution=0.3)
sc.tl.rank_genes_groups(adata_12, 'leiden', only_positive=True, n_genes=200, use_raw=False)
file = out_dir + 'adata_12_cluster_gene_markers.csv'
df = sc.get.rank_genes_groups_df(adata_12, group=None)
df.to_csv(file)

# Pathway Analysis
anova_results_file = out_dir + 'adata_12_anova_results.csv' 
pathway_results_file = out_dir + 'adata_12_pathway_results.csv'
pa_result = pa.reactome_aucell(adata_12, reactome_gmt_file)
anova_results = pa_result.perform_pathway_anova()
anova_results.to_csv(anova_results_file)
pa_result.adata.obsm['X_aucell'].to_csv(pathway_results_file)
del(pa_result)

# PAGA/CytoTRACE
# output for CytoTRACE run on web
# https://cytotrace.stanford.edu/
cyto_df_12 = adata_12.copy()
cyto_df_12.X = cyto_df_12.layers['counts']
cyto_df_12 = cyto_df_12.to_df()
cyto_df_12.transpose().to_csv(out_dir + 'adata_12_cyto_input_counts.csv')

# results from CytoTRACE run on web - output saved and renamed
cyto_web_results_12 = '../resources/CytoTRACE_results_e12_counts.txt'
df12 = pd.read_csv(cyto_web_results_12, sep='\t')
df12.index = df12.index.str.replace('.','-')
adata_12.obs['CytoTRACE'] = df12.CytoTRACE

# Get root cell
root_cell_12 = adata_12.obs.CytoTRACE.idxmax()
adata_12.uns['iroot'] = adata_12.obs.index.get_loc(root_cell_12)
del(df12)

# Calc PAGA and DPT scores for 12.5 cells with reduced clusters
sc.tl.paga(adata_12)
sc.tl.dpt(adata_12)
adata_12.uns['paga']['pos'] = _reset_paga_pos(adata_12) # Reset PAGA coords to median UMAP coords
sc.pl.paga(adata_12, pos=adata_12.uns['paga']['pos'], plot=False) #If False, do not create the figure, simply compute the layout.

# # Get PAGA and DPT scores for 12.5 cells
# run_paga(adata_12)
# run_dpt(adata_12, root_cell_12)

sc.pl.paga(adata_12,
           threshold=0.09,
           color = ['leiden', 'CytoTRACE', 'dpt_pseudotime'],
           pos = adata_12.uns['paga']['pos'],
           cmap = 'viridis',
           plot=False) # If False, do not create the figure, simply compute the layout.

#### Build additional E12.5 dataset for Defa24 visualization ###
adata_cat = adata_12.concatenate(adata_17, batch_categories=['E12.5','E17.5'], join='outer')

# Add Sid's cell type annotation
load_cell_annotations(adata_17,CELL_TYPE_ANNOTATION_FILE, 17)
load_cell_annotations(adata_12, CELL_TYPE_ANNOTATION_FILE, 12)


### Build Harmony dataset ###
h17, h12, hm = merge(adata_17, adata_12, batch_categories=batch_categories, batch_correction='harmony')
hm.uns['batch_colors'] = ["#28c941", "#ff8880"]
del(hm.obsm['X_diffmap']) # prevent neighbors calc from using X_diffmap
sc.pp.neighbors(hm, use_rep='X_pca_harmony', random_state=random_state)
sc.tl.leiden(hm, random_state=random_state)
sc.tl.paga(hm)
sc.pl.paga(hm, random_state=random_state, plot=False)
sc.tl.umap(hm, random_state=random_state, init_pos='paga')

### Write AnnDatas to disk ###
adata_17.write(out_dir + 'adata_17.h5ad')
adata_12.write(out_dir + 'adata_12.h5ad')
adata_cat.var = adata_cat.var.astype({
    'mt-E12.5':'bool',
    'highly_variable-E12.5':'bool',
    'mt-E17.5':'bool',
    'highly_variable-E17.5':'bool'})    
adata_cat.write(out_dir + 'adata_12_defa.h5ad')
hm.write(out_dir + 'harmony_dataset.h5ad')
