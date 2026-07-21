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


adatalognorm = adata_raw.copy()


adatalognorm.obs['n_counts'] = adatalognorm.X.sum(axis=1).A1


sc.pl.violin(adatalognorm, ['n_genes', 'n_counts', 'percent_mito'], jitter=0, multi_panel=True, stripplot = False, save=obj_name)

sc.pl.highest_expr_genes(adatalognorm, n_top=20, save=obj_name)

sc.pp.filter_cells(adatalognorm, min_genes=200)
sc.pp.filter_genes(adatalognorm, min_cells=50)

adatalognorm.var_names[adatalognorm.var_names.str.startswith('mt-')]

mito_genes = adatalognorm.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adatalognorm.obs['percent_mito'] = np.sum(
    adatalognorm[:, mito_genes].X, axis=1).A1 / np.sum(adatalognorm.X, axis=1).A1


sc.pl.violin(adatalognorm, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0, multi_panel=True, stripplot = False, save=obj_name+'afterfilter')

sc.pl.scatter(adatalognorm, x='n_counts', y='percent_mito', save=obj_name+'mito')
sc.pl.scatter(adatalognorm, x='n_counts', y='n_genes', save=obj_name+'ngenes')

sc.logging.print_versions()

sc.pp.normalize_total(adatalognorm, target_sum=1e4)

sc.pp.log1p(adatalognorm)

adatalognorm.raw = adatalognorm.copy()


sc.pp.regress_out(adatalognorm, ['n_counts', 'percent_mito'])


sc.pp.scale(adatalognorm, max_value=10)


sc.tl.pca(adatalognorm, svd_solver='arpack')


sc.pl.pca(adatalognorm, color='Wif1', save=obj_name+'wif1')


sc.pl.pca_variance_ratio(adatalognorm, log=True, save=obj_name)


adatalognorm.write(args.save_file)