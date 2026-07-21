# Single-cell RNA-sequencing of cardiac stromal cells following ischemic injury

**MSc thesis — Elsa Lawrence**
MSc Genes, Drugs and Stem Cells: Novel Therapies, Imperial College London (2020)
Cardiac Function Research Group, National Heart and Lung Institute (NHLI), Hammersmith Campus

> *Single cell RNA-sequencing captures phenotypic changes in the repertoire of cardiac stromal cells following ischemic injury*

This repository shares my master's thesis and a selection of the analysis code behind it. It is intended to give an overview of the project and a look at how the single-cell analysis was done — not to be a fully reproducible pipeline, since the underlying data is not redistributed here (see [Data](#data)).

---

## Abstract

Myocardial infarction (MI) is a leading cause of death worldwide. Cardiac fibroblasts play a leading role in wound healing in response to MI following cardiomyocyte death, by promoting maladaptive remodelling which contributes to loss of cardiac function. Reperfusion therapies, currently the first choice of treatment, restore blood flow but are insufficient to restore the lost cardiac muscle tissue and function. A better understanding of the cardiac fibroblast cell landscape is crucial for the development of strategies to induce cardiac repair. Single-cell RNA-sequencing (scRNA-seq) provides a new tool to explore the cardiac stromal compartment and understand the cellular drivers of adverse cardiac remodelling with higher resolution and precision.

Here, we computationally analysed a mix of newly generated and publicly available scRNA-seq data from unoperated, sham-operated, and operated mouse hearts collected at various timepoints following experimental MI. We used low-resolution Leiden clustering to identify and annotate **10 major cell types** within a whole dataset of **75,145 cells**. Further sub-clustering identified **7 subpopulations** within the cardiac stromal fibroblast compartment (**42,181 cells**) and **12 subpopulations** within the macrophage compartment. Within fibroblasts, specificity of *Wif1* expression was confirmed globally. Interactions were computationally predicted between fibroblasts and macrophages using CellPhoneDB, identifying a key receptor–ligand pair, **CD74/MIF**, which might be a promising target for cardioprotective mechanisms. Overall, this analysis gives a deeper insight into the heterogeneity of cardiac stromal fibroblast cells and their dynamic changes in gene expression and function post-MI.

## Key findings

- Annotation of **10 major cell types** across 75,145 cells from unoperated, sham-operated and post-MI mouse hearts.
- **7 fibroblast subpopulations** and **12 macrophage subpopulations** resolved by sub-clustering.
- Confirmation of *Wif1* as a fibroblast-specific marker.
- CellPhoneDB prediction of fibroblast–macrophage crosstalk, highlighting the **CD74/MIF** receptor–ligand pair.

## Methods overview

The analysis was carried out in Python using the [Scanpy](https://scanpy.readthedocs.io/) ecosystem, with the following main steps:

1. **Loading & merging** 10x Genomics count matrices from newly generated and public datasets.
2. **Quality control & preprocessing** — doublet detection (Scrublet), cell/gene filtering, mitochondrial-content filtering, normalisation, log-transformation, regression and scaling.
3. **Batch integration** with BBKNN.
4. **Dimensionality reduction & clustering** — PCA, UMAP and Leiden clustering at low resolution for major cell types, then higher resolution for sub-clustering.
5. **Marker gene / differential expression** analysis (`rank_genes_groups`).
6. **Cell–cell interaction prediction** between fibroblasts and macrophages using [CellPhoneDB](https://www.cellphonedb.org/).

**Tools:** `scanpy`, `anndata`, `scrublet`, `bbknn`, `scvelo`, `scipy`, `numpy`, `pandas`, `matplotlib`, and CellPhoneDB.

## Repository structure

```
.
├── README.md
├── LICENSE
├── environment.yml            # conda environment for the Python analysis
├── docs/
│   └── Lawrence_MSc_thesis_2020.pdf
└── code/
    └── python/
        ├── preprocessing.py                  # QC, filtering, normalisation, PCA (CLI script)
        ├── preprocessing_hvg.py              # preprocessing variant using highly variable genes
        ├── 00_getting_raw_data.ipynb         # loading data and building the AnnData object
        ├── 01_load_merge_annotate.ipynb      # merging datasets, annotating cells, subsetting fibroblasts
        ├── 02_fibroblast_subclustering.ipynb # fibroblast sub-clustering, markers, DEG analysis
        ├── 03_macrophage_subclustering.ipynb # macrophage sub-clustering, markers, SCCAF assessment
        ├── 04_cellphonedb_interactions.ipynb # fibroblast–macrophage interaction analysis (CellPhoneDB output)
        └── 05_cluster_assessment.ipynb       # choosing cluster resolution: elbow/SSE, silhouette, CH, gap statistic, SCCAF; correlation matrices
```

> **Note on scope.** This is a representative selection of the code, chosen to show the workflow. Personal HPC file paths have been replaced with generic `data/` placeholders, so the notebooks are for reading rather than direct execution.
>
> **External dependencies in the notebooks.** A couple of the analysis notebooks call tools that aren't part of this repository: [SCCAF](https://github.com/SCCAF/sccaf) (for cluster self-projection accuracy, used in notebooks 03 and 05) and an in-house `ClusterSupport` module for the gap-statistic step (notebook 05). Install SCCAF separately if you want to run those cells; the `ClusterSupport` import can be skipped.

## Data

The published figures and full results are in the thesis PDF (`docs/`). The raw sequencing data is **not** included in this repository:

- **Public data** — mouse cardiac interstitial cells from sham-operated and 3/7-day post-MI hearts are available from ArrayExpress (EMBL-EBI) under accession **[E-MTAB-7376](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7376)** (Farbehi *et al.*, 2019).
- **Lab-generated data** — the PDGFRα-fated scRNA-seq datasets were generated by and belong to the Noseda lab and are not redistributed here.

## Citation

Lawrence, E. (2020). *Single cell RNA-sequencing captures phenotypic changes in the repertoire of cardiac stromal cells following ischemic injury.* MSc thesis, Imperial College London.

## Acknowledgements

Supervised by Dr Michela Noseda, Dr Antonio M. A. Miranda and Prof. Michael Schneider (Cardiac Function Research Group, NHLI, Imperial College London), with programming guidance from Dr Michael Lee, Dr Alik Huseynov and Matt Greenig, and the wider Noseda lab.

## License

Code in this repository is released under the [MIT License](LICENSE). The thesis document (`docs/`) is © Elsa Lawrence, 2020, and is shared for reference; please cite it rather than reproducing it.
