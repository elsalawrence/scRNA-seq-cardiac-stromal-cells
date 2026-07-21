import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import os
import bbknn as bk
import scvelo as scv
import argparse
import anndata
sc.settings.set_figure_params(dpi=300,fontsize=10)


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file_name', type=str, help='File name for an h5ad AnnData object')
parser.add_argument('-s', '--save_file', type=str, help='Name of file you want to save')
args = parser.parse_args()
obj_name = '_' + args.file_name.rstrip('.h5ad')

adata_raw = sc.read_h5ad(args.file_name)


adata_raw.X.min()


hvg_adatalognorm = adata_raw.copy()


hvg_adatalognorm.obs['n_counts'] = hvg_adatalognorm.X.sum(axis=1).A1


sc.pl.violin(hvg_adatalognorm, ['n_genes', 'n_counts', 'percent_mito'], jitter=0, multi_panel=True, stripplot = False, save=obj_name)

sc.pl.highest_expr_genes(hvg_adatalognorm, n_top=20, save=obj_name)

sc.pp.filter_cells(hvg_adatalognorm, min_genes=200)
sc.pp.filter_genes(hvg_adatalognorm, min_cells=50)

hvg_adatalognorm.var_names[hvg_adatalognorm.var_names.str.startswith('mt-')]

mito_genes = hvg_adatalognorm.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
hvg_adatalognorm.obs['percent_mito'] = np.sum(
    hvg_adatalognorm[:, mito_genes].X, axis=1).A1 / np.sum(hvg_adatalognorm.X, axis=1).A1


sc.pl.violin(hvg_adatalognorm, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0, multi_panel=True, stripplot = False, save=obj_name+'afterfilter')

sc.pl.scatter(hvg_adatalognorm, x='n_counts', y='percent_mito', save=obj_name+'mito')
sc.pl.scatter(hvg_adatalognorm, x='n_counts', y='n_genes', save=obj_name+'ngenes')

sc.logging.print_versions()

sc.pp.normalize_total(hvg_adatalognorm, target_sum=1e4)

sc.pp.log1p(hvg_adatalognorm)

hvg_adatalognorm.raw = hvg_adatalognorm.copy()

sc.pp.highly_variable_genes(hvg_adatalognorm, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(hvg_adatalognorm)

hvg_adatalognorm = hvg_adatalognorm[:, hvg_adatalognorm.var.highly_variable]

sc.pp.regress_out(hvg_adatalognorm, ['n_counts', 'percent_mito'])


sc.pp.scale(hvg_adatalognorm, max_value=10)


sc.tl.pca(hvg_adatalognorm, svd_solver='arpack')


sc.pl.pca(hvg_adatalognorm, color='Wif1', save=obj_name+'wif1')


sc.pl.pca_variance_ratio(hvg_adatalognorm, log=True, save=obj_name)


hvg_adatalognorm.write(args.save_file)