"""
This script is used to generate some plots for Missy's manuscript related to the mouse intestinal stem cell development
project
"""
import os
import os.path as path

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

DPI = 720 
# Need a customized map
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt

# To be used as marker genes for scoring clusters
MOUSE_CELL_MARKERS = './resources/Mouse_cell_markers.txt'

colors_custom = ['#cfcfcf', "#3300FF"] # Light grey to blue
cmap_custom = LinearSegmentedColormap.from_list("mycmap", colors_custom)
cmap_custom = 'viridis_r'

# Make sure this points to the right file
CELL_TYPE_ANNOTATION_FILE = '../resources/12_17_cell_annotations.csv'
batch_categories = ['E17.5', 'E12.5']
random_state = 17051256

def load_cell_annotations(
    annData: sc.AnnData,
    cell_annot_file: str = CELL_TYPE_ANNOTATION_FILE,
    timepoint: int = None)->None:
    # append cell annotations to adatas, need to alter Seurat cell ids and select by batch
    
    cell_annot_df = pd.read_csv(cell_annot_file)
    cell_annot_df['cell'] = cell_annot_df.cell.str.split('_').apply(lambda x: ''.join((x[2],'-1')))
    if timepoint is None:
        print('need timepoint')
    else:
        df = cell_annot_df[cell_annot_df.timepoint==timepoint].set_index('cell')
        annData.obs['cell_type'] = df['annotated_clusters']

# The following functions are modified from code hosted at fi_sc_analysis
def _reset_paga_pos(adata):
    """
    The following code is based on scatterplots.py to calculate paga positions directly from
    the umap. The original paga's pos is not matched with umap because of different layout algorithm
    used.
    :param adata:
    :return: the median coordinates of all cells in individual clusters.
    """
    key = 'leiden'
    umap_key = "X_umap"
    clusters = adata.obs[key].cat.categories
    cluster_pos = np.zeros((len(clusters), 2))
    for icluster, cluster in enumerate(clusters) :
        selected = adata.obs[key] == cluster
        selected_umap = adata.obsm[umap_key][selected, :]
        x_pos, y_pos = np.median(selected_umap, axis=0)
        cluster_pos[icluster] = [x_pos, y_pos]
    return cluster_pos


def run_paga(adata:sc.AnnData):
    sc.pp.neighbors(adata, random_state=random_state)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False, random_state=random_state, show=False)
    adata.uns['paga']['pos'] = _reset_paga_pos(adata)


def run_umap_via_paga(adata:sc.AnnData, use_rep = None):
    """
    Follow the procedure outlined in https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
    by running paga first as the init_pos.
    :param adata: assumed this data has been cleaned and normalized.
    :param use_rep: for the harmony batch corrected data, X_pca_harmony should be used.
    :return:
    """
    if use_rep is None:
        sc.pp.highly_variable_genes(adata)
        sc.pp.pca(adata, random_state=random_state)
        use_rep = 'X_pca'
    sc.pp.neighbors(adata, use_rep=use_rep, random_state=random_state)
    sc.tl.leiden(adata, random_state=random_state)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False, random_state=random_state, show=False)
    sc.tl.umap(adata, random_state=random_state, init_pos='paga')
    adata.uns['paga']['pos'] = _reset_paga_pos(adata)


def plot_paga(anndata: sc.AnnData):
    # This is the root cell selected based on the maximum cytotrace value
    # This cell is also used for dpt inference
    root_cell = 'GTGCACGCACCCAAGC-1-1-e17_5'
    root_index = np.flatnonzero(anndata.obs_names == root_cell)[0]
    sc.pl.paga(anndata,
               threshold=0.45,
               color = ['leiden', 'cytotrace', 'dpt_pseudotime'],
               pos = anndata.uns['paga']['pos'],
               cmap = cmap_custom,
               root = root_index) # Don't make any difference by assigning root


def run_dpt(adata:sc.AnnData,
            root_cell: str):
    # Get the number index of the root_cell
    root_index = np.flatnonzero(adata.obs_names == root_cell)[0]
    adata.uns['iroot'] = root_index
    sc.tl.dpt(adata)


def open_data(file: str,
              preprocess: bool = True) -> sc.AnnData:
    """
    Open data the data and then pre-process it
    :param file:
    :param preprocess: true for preprocess
    :return:
    """
    adata = None
    if path.isdir(file):
        adata = sc.read_10x_mtx(file, var_names='gene_symbols', cache=True)
    else: # This should be a file. Assume it is a h5ad
        adata = sc.read_h5ad(file)
    if adata is None:
        raise ValueError('Cannot open file {}.'.format(str))
    if not preprocess:
        return adata
    # Copy from the scanpy_wrapper and from Chris' code
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # Have not done any filtering based on total gene numbers
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    adata.layers['counts'] = adata.X
    sc.pp.normalize_total(adata, 1E+4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    return adata


def run_umap(adata,
             need_run_variable: bool=False):
    if need_run_variable:
        sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata, random_state=random_state)
    sc.pp.neighbors(adata, random_state=random_state)
    sc.tl.leiden(adata, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    return adata


def merge(source: sc.AnnData,
          target: sc.AnnData,
          batch_categories: str = None,
          batch_correction: str = None):
    """
    Merge the source dataset into the target dataset
    :param source:
    :param target:
    :param batch_categories:
    :param batch_correction:
    :return:
    """
    merged = source.concatenate(target, batch_categories=batch_categories)
    if batch_correction is not None:
        sc.pp.pca(merged) # Need to conduct pca first for both batch correction
        if batch_correction == 'harmony':
            sce.pp.harmony_integrate(merged, 'batch')
        elif batch_correction == 'bbknn':
            # See https://github.com/theislab/scanpy/issues/514
            import bbknn
            bbknn.bbknn(merged, 'batch')
        else:
            raise ValueError("{} is not supported in batch correction.".format(batch_correction))
    # Want to split them into the original data for doing whatever it is needed
    if batch_categories is None:
        batch_categories = ['0', '1']
    src_data = merged[merged.obs['batch'] == batch_categories[0]]
    target_data = merged[merged.obs['batch'] == batch_categories[1]]
    return src_data, target_data, merged


def squarify_umap(plot, adata, colorbar=False, is_seurat=False):
    # modify default scanpy figures for manuscript
    xmin, xmax = min(adata.obsm['X_umap'][:,0]), max(adata.obsm['X_umap'][:,0])
    ymin, ymax = min(adata.obsm['X_umap'][:,1]), max(adata.obsm['X_umap'][:,1])
    xrange = xmax - xmin
    yrange = ymax - ymin
    adjusty = (xrange - yrange)/2
    plot.xaxis.set_major_locator(plt.MaxNLocator(6))
    plot.yaxis.set_major_locator(plt.MaxNLocator(6))
    plot.set(adjustable='box', aspect='equal')
    if is_seurat:
        adjusty = (yrange - xrange)/2
        plot.set_xlim(plot.get_xlim()[0] - adjusty, plot.get_xlim()[1] + adjusty)
    else:
        if adjusty > 0:
            plot.set_ylim(plot.get_ylim()[0] - adjusty, plot.get_ylim()[1] + adjusty)
        else:
            plot.set_xlim(plot.get_xlim()[0] + adjusty, plot.get_xlim()[1] - adjusty)
    return plot    


