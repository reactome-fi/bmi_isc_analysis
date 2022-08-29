import os
os.chdir(os.getcwd())
out_dir =  "bmi_figuresl_script/"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
os.chdir(out_dir)
os.getcwd()

exec(open('../PlotGenerators.py').read())

# 
markers12 = pd.read_csv('../../resources/e12-grp0-markers.csv')
markers17 = pd.read_csv('../../resources/e17-grp0-4-markers.csv')
markersh = pd.read_csv('../../resources/stem-h-markers.csv')
markersirr = pd.read_csv('../../resources/stem-irr-markers.csv')

# Processed file paths - set to use frozen adatas to avoid any randomness in UMAPs
seurat_file = '../../resources/seurat_aHom_aIrr.h5ad'
adata12_file = '../../resources/adata_12_frozen.h5ad'
adata17_file = '../../resources/adata_17_frozen.h5ad'
harmony_file = '../../resources/harmony_dataset_frozen.h5ad'

seurat = sc.read_h5ad(seurat_file)
e12 = sc.read_h5ad(adata12_file)
e17 = sc.read_h5ad(adata17_file)
hm = sc.read_h5ad(harmony_file)

# Figure 4
## A
hm.uns['batch_colors'] = ['#28c941', '#a2a2a2',]
plot = sc.pl.umap(hm,color='batch', show=False)
plot = squarify_umap(plot, hm)
plt.savefig('4A2_harmony_square_stem_colored.tiff', bbox_inches='tight', dpi=495, facecolor='white')

## B
adata = sc.read(adata12_file) # for reproducing figures (otherwise UMAP coords vary slightly)
adata.uns['leiden_colors'] = ['#dedede', '#bbbbbb', '#888888', '#555555', '#111111']
plot = sc.pl.umap(adata, color='leiden', legend_loc=4, show=False)
plot = squarify_umap(plot, adata)
plt.savefig('4B_12_5_umap_leiden.tiff', bbox_inches='tight', dpi=600, facecolor='white')

## C
plot = sc.pl.paga(adata, color=['CytoTRACE'], pos=adata.uns['paga']['pos'], fontsize=0, edge_width_scale=0.5, random_state=random_state, show=False)
xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2
plot.xaxis.set_major_locator(plt.MaxNLocator(5))
plot.yaxis.set_major_locator(plt.MaxNLocator(5))
plot.set(adjustable='box', aspect='equal')
plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
plot.set_xbound(xmin-1, xmax-1)
plot.set_ybound(ymin-adjusty, ymax+adjusty)
cbar = plot.figure.properties()['default_bbox_extra_artists'][1]
cbar.set_position((0.646,0.08,0.03,0.8))

fig = plot.figure
for i in range(5):
    fig.properties()['default_bbox_extra_artists'][4].set_visible(False) 
dummy = plt.figure()
new_manager = dummy.canvas.manager
new_manager.canvas.figure = fig
fig.set_canvas(new_manager.canvas)
plt.savefig('4C_12_5_paga_cyto.tiff', bbox_inches='tight', dpi=600, facecolor='white')

## D
adata = sc.read_h5ad(adata17_file)
testplot = [sc.pl.umap(adata, color='leiden', legend_loc=4, show=False)]
xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2
for plot in testplot:
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
fig = testplot[0].figure
dummy = plt.figure()
new_manager = dummy.canvas.manager
new_manager.canvas.figure = fig
fig.set_canvas(new_manager.canvas)
plt.savefig('4D_17_5_leiden.tiff', bbox_inches='tight', dpi=600, facecolor='white')

## E
plot = sc.pl.paga(adata,
           #threshold=0.09,
           color = ['cyto_web17'],#CytoTRACE
           pos = adata.uns['paga']['pos'],
           cmap = 'viridis',
           root = adata.uns['iroot'],
           #fontsize=0,
           edge_width_scale=0.5,
           random_state=random_state,
           show=False)
xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2
plot.xaxis.set_major_locator(plt.MaxNLocator(6))
plot.yaxis.set_major_locator(plt.MaxNLocator(5))
plot.set(adjustable='box', aspect='equal')
plot.set_xbound(xmin-1, xmax+1)
plot.set_ylim(plot.get_ylim()[0] - adjusty-1, plot.get_ylim()[1] + adjusty+1)
cbar = plot.figure.properties()['default_bbox_extra_artists'][1]
cbar.set_position((0.646,0.08,0.03,0.8))
fig = plot.figure

for i in range(18):
    fig.properties()['default_bbox_extra_artists'][4].set_visible(False) # use idx 3 for removing leiden text
dummy = plt.figure()
new_manager = dummy.canvas.manager
new_manager.canvas.figure = plot.figure
fig.set_canvas(new_manager.canvas)
plt.savefig('4E_17_5_cyto.tiff', bbox_inches='tight', dpi=600, facecolor='white')

plot = sc.pl.paga(adata,
           #threshold=0.09,
           color = ['dpt_pseudotime'],
           pos = adata.uns['paga']['pos'],
           cmap = 'viridis',
           root = adata.uns['iroot'],
           #fontsize=0,
           edge_width_scale=0.5,
           random_state=random_state,
           show=False)
xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2
plot.xaxis.set_major_locator(plt.MaxNLocator(6))
plot.yaxis.set_major_locator(plt.MaxNLocator(5))
plot.set(adjustable='box', aspect='equal')
plot.set_xbound(xmin-1, xmax+1)
plot.set_ylim(plot.get_ylim()[0] - adjusty-1, plot.get_ylim()[1] + adjusty+1)
cbar = plot.figure.properties()['default_bbox_extra_artists'][1]
cbar.set_position((0.646,0.08,0.03,0.8))
fig = plot.figure
for i in range(18):
    fig.properties()['default_bbox_extra_artists'][4].set_visible(False) # use idx 3 for removing leiden text

dummy = plt.figure()
new_manager = dummy.canvas.manager
new_manager.canvas.figure = plot.figure
fig.set_canvas(new_manager.canvas)
plt.savefig('4E_17_5_dpt.tiff', bbox_inches='tight', dpi=600, facecolor='white')

## F & J
adata17 = sc.read_h5ad(adata17_file)
adata12 = sc.read_h5ad(adata12_file)

batch_categories = ['E12.5', 'E17.5']
adata_cat = adata12.concatenate(adata17, batch_categories=batch_categories, join='outer')

genes = ['Pou5f1', 'Alpi', 'Chga', 'Dkk3','Sfrp5','Wnt3', 'Wnt5a']
max_expr = {}
for i in genes:
    max_expr[i] = max([adata17[:,i].X.max(), adata12[:,i].X.max()])

pou_max = max([adata17[:,'Pou5f1'].X.max(), adata12[:,'Pou5f1'].X.max()])
alpi_max = max([adata17[:,'Alpi'].X.max(), adata12[:,'Alpi'].X.max()])
defa_max = max([adata17[:,'Defa24'].X.max()]) #, adata12[:,'Defa24'].X.max()]
chga_max = max([adata17[:,'Chga'].X.max(), adata12[:,'Chga'].X.max()])

### 17.5
ticks = [-10, 0, 10, 20]
xmin, xmax = min(adata17.obsm['X_umap'][:,0]), max(adata17.obsm['X_umap'][:,0])
ymin, ymax = min(adata17.obsm['X_umap'][:,1]), max(adata17.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2
genes = ['Defa24']
for i in genes:
    plot = sc.pl.umap(adata17, color=i, cmap='viridis_r', vmax=defa_max, show=False)
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
    file_name = '4_17_5_' + i + '.tiff'
    plt.savefig(file_name, bbox_inches='tight', dpi=600, facecolor='white')

genes = ['Pou5f1', 'Alpi', 'Chga', 'Dkk3','Sfrp5','Wnt3', 'Wnt5a']
for i in genes:
    plot = sc.pl.umap(adata17, color=i, cmap='viridis_r', vmax=max_expr[i], show=False)
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
    file_name = '4_17_5_' + i + '.tiff'
    plt.savefig(file_name, bbox_inches='tight', dpi=600, facecolor='white')

### 12.5
ticks12 = [-4, 0, 4, 8, 12]
xmin12, xmax12 = min(adata12.obsm['X_umap'][:,0]), max(adata12.obsm['X_umap'][:,0])
ymin12, ymax12 = min(adata12.obsm['X_umap'][:,1]), max(adata12.obsm['X_umap'][:,1])
xrange12 = xmax12 - xmin12
yrange12 = ymax12 - ymin12
adjusty12 = (xrange12 - yrange12)/2

genes = ['Pou5f1', 'Alpi', 'Chga', 'Dkk3','Sfrp5','Wnt3', 'Wnt5a']
for i in genes:
    plot = sc.pl.umap(adata12, color=i, cmap='viridis_r', vmax=max_expr[i], show=False)
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    plot.set_ylim(plot.get_ylim()[0] - adjusty12, plot.get_ylim()[1] + adjusty12)
    file_name = '4_12_5_' + i + '.tiff'
    plt.savefig(file_name, bbox_inches='tight', dpi=600, facecolor='white')

batch_categories = ['E12.5', 'E17.5']
plot = sc.pl.umap(adata_cat[adata_cat.obs.batch=='E12.5'], color='Defa24', cmap='viridis_r', vmax=defa_max, show=False)
plot.set_xticks(ticks12)
plot.set_yticks(ticks12)
plot.set(adjustable='box', aspect='equal')
plt.savefig('4_12_5_Defa.tiff', bbox_inches='tight', dpi=600, facecolor='white')

# Fig 4 & Sup Fig 5
## Pathways

pathway_list = [
    'mRNA Splicing', "Collagen degradation", 'Cobalamin (Cbl, vitamin B12) transport and metabolism', 'Digestion of dietary lipid',
    'Fatty acid metabolism', 'Growth hormone receptor signaling', 'HDMs demethylate histones', 'RUNX3 regulates YAP1-mediated transcription']

for i in pathway_list:
    adata17.obs[i] = adata17.obsm['X_aucell'][i]

pathway_list12 = ['mRNA Splicing', 'Digestion of dietary lipid','Growth hormone receptor signaling', 'HDMs demethylate histones', 'RUNX3 regulates YAP1-mediated transcription']

for i in pathway_list12:
    adata12.obs[i] = adata12.obsm['X_aucell'][i]

# a number of pathways not scored in adata12 (<80% genes present), set to 0
adata12.obs['Collagen degradation'] = 0
adata12.obs['Cobalamin (Cbl, vitamin B12) transport and metabolism'] = 0
adata12.obs['Fatty acid metabolism'] = 0

pathway_max = {}
for i in pathway_list:
    pathway_max[i] = max_expr[i] = max([adata17.obs[i].max(), adata12.obs[i].max()])
# set pathway max expr scale to 12.5 data
pathway_max['RUNX3 regulates YAP1-mediated transcription'] = adata12.obs[i].max()
pathway_max['Oncogene Induced Senescence'] = adata12.obs[i].max()
pathway_max['Chromatin modifying enzymes'] = adata12.obs[i].max()

for i in pathway_list:
    plot = sc.pl.umap(adata17, color=i, cmap='viridis_r', vmax=pathway_max[i], show=False)
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
    file_name = '4G_17_5_' + str(i)[:5] + '.tiff'
    plt.savefig(file_name, bbox_inches='tight', dpi=600, facecolor='white')

for i in pathway_list:
    plot = sc.pl.umap(adata12, color=i, cmap='viridis_r', vmax=pathway_max[i], show=False)
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    plot.set_ylim(plot.get_ylim()[0] - adjusty12, plot.get_ylim()[1] + adjusty12)
    file_name = '4G_12_5_' + str(i)[:5]  + '.tiff'
    plt.savefig(file_name, bbox_inches='tight', dpi=600, facecolor='white')

## B & C
adata = sc.read(adata12_file) # for reproducing figures (otherwise UMAP coords vary slightly)
plot = sc.pl.paga(adata, color=['dpt_pseudotime'], pos=adata.uns['paga']['pos'], fontsize=0, edge_width_scale=0.5, random_state=random_state, show=False)
xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2

plot.xaxis.set_major_locator(plt.MaxNLocator(5))
plot.yaxis.set_major_locator(plt.MaxNLocator(5))
plot.set(adjustable='box', aspect='equal')
plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
plot.set_xbound(xmin-1, xmax-1)
plot.set_ybound(ymin-adjusty, ymax+adjusty)
cbar = plot.figure.properties()['default_bbox_extra_artists'][1]
cbar.set_position((0.646,0.08,0.03,0.8))

fig = plot.figure
for i in range(5):
    fig.properties()['default_bbox_extra_artists'][4].set_visible(False) 
dummy = plt.figure()
new_manager = dummy.canvas.manager
new_manager.canvas.figure = fig
fig.set_canvas(new_manager.canvas)
plt.savefig('5B_12_5_paga_dpt.tiff', bbox_inches='tight', dpi=600, facecolor='white')

## D & E
adata = sc.read(adata17_file) # for reproducing figures (otherwise UMAP coords vary slightly)
plot = sc.pl.paga(adata, color=['dpt_pseudotime'], pos=adata.uns['paga']['pos'], fontsize=0, edge_width_scale=0.5, random_state=random_state, show=False)

xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
xrange = xmax - xmin
yrange = ymax - ymin
adjusty = (xrange - yrange)/2

plot.xaxis.set_major_locator(plt.MaxNLocator(5))
plot.yaxis.set_major_locator(plt.MaxNLocator(5))
plot.set(adjustable='box', aspect='equal')
plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
plot.set_xbound(xmin-1, xmax-1)
plot.set_ybound(ymin-adjusty, ymax+adjusty)
cbar = plot.figure.properties()['default_bbox_extra_artists'][1]
cbar.set_position((0.646,0.08,0.03,0.8))

fig = plot.figure
for i in range(5):
    fig.properties()['default_bbox_extra_artists'][4].set_visible(False) 
dummy = plt.figure()
new_manager = dummy.canvas.manager
new_manager.canvas.figure = fig
fig.set_canvas(new_manager.canvas)
plt.savefig('5B_17_5_paga_dpt.tiff', bbox_inches='tight', dpi=600, facecolor='white')

## G
from matplotlib_venn import venn2

marks12 = set(markers12[(markers12.logfoldchanges > 2) & (markers12.pvals_adj < 0.01)]['names'])
marks17 = set(markers17[(markers17.logfoldchanges > 2) & (markers17.pvals_adj < 0.01)]['names'])
markshom = set(markersh[(markersh.logfoldchanges > 2) & (markersh.pvals_adj < 0.01)]['names'])
marksirr = set(markersirr[(markersirr.logfoldchanges > 2) & (markersirr.pvals_adj < 0.01)]['names'])

col = {
    'e12':"#dedede",
    'e17':"#157023",
    'hom': "#ff0000",
    'irr':"#ff7184"}

venn2([marks12, marks17], set_labels = ('E12.5', 'E17.5'), set_colors=[col['e12'], col['e17']])
plt.title("DEGs: logFC > 2 & FDR < 0.01")
plt.savefig('logFC-fdr-e12e17', dpi=600)

## H
e12_ids = adata12[adata12.obs.leiden=='0'].obs_names.str[:] + '-E12.5'
e17_ids = adata17[adata17.obs.leiden.isin(['0','4'])].obs_names 

hm.obs.batch.cat.add_categories('E12.5 - Stem', inplace=True)
hm.obs.batch.cat.add_categories('E17.5 - Stem', inplace=True)
hm.obs.loc[e12_ids, 'batch'] = 'E12.5 - Stem'
hm.obs.loc[e17_ids, 'batch'] = 'E17.5 - Stem'
hm.uns['batch_colors'] = ['#a2a2a2', '#dedede', '#28c941', '#157023']
hm.obs.batch.cat.reorder_categories(['E12.5', 'E12.5 - Stem', 'E17.5', 'E17.5 - Stem'], inplace=True)

plot = sc.pl.umap(hm,color='batch', show=False)
plot = squarify_umap(plot, hm)
plt.savefig('4A2_harmony_square_stem_colored.tiff', bbox_inches='tight', dpi=495, facecolor='white')

# Figure 5
## C
venn2([marks12, markshom], set_labels = ('E12.5', 'Adult'), set_colors=[col['e12'], col['hom']])
plt.title("DEGs: logFC > 2 & FDR < 0.01")
plt.savefig('5C_E12hom', dpi=600)

venn2([marks17, markshom], set_labels = ('E17.5', 'Adult'), set_colors=[col['e17'], col['hom']])
plt.title("DEGs: logFC > 2 & FDR < 0.01")
plt.savefig('5C_E17hom', dpi=600)

# Figure 6
## E
venn2([markshom, marksirr], set_labels = ('Adult', 'Adult IR'), set_colors=[col['hom'], col['irr']])
plt.title("DEGs: logFC > 2 & FDR < 0.01")
plt.savefig('6E_HomvIrr', dpi=600)


