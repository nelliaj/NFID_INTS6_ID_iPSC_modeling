import hnoca.snapseed as snap
from hnoca.snapseed.utils import read_yaml
import scanpy as sc
import harmonypy as hm
import pandas as pd

# Read in the marker genes (imported from "https://github.com/theislab/neural_organoid_atlas/blob/main/supplemental_files/Data_S1_snapseed_markers.yaml")

marker_genes = read_yaml("deconvolution/Data_S1_snapseed_markers_updated.yaml")
adata=sc.read_h5ad("deconvolution/merged_seurat_coculture_subset_v3.h5ad")

# Filter and preprocessing D28 Ngn2 single-cell dataset from Pietiläinet et al. 2023, Cell Reports: https://doi.org/10.1016/j.celrep.2022.111988

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt',], percent_top=None, log1p=True, inplace=True)

adata.layers["counts"]=adata.X.copy()

#Processing
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

adata.layers["lognorm"]=adata.X.copy()

#Find highly variable genes
sc.pp.highly_variable_genes(adata)
#adata = adata[:, adata.var.highly_variable]

adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver='arpack')


# PCA matrix: cells × PCs
X = adata_hvg.obsm['X_pca']  
meta = pd.DataFrame({'orig.ident': adata_hvg.obs['orig.ident'].astype(str).values})
vars_use = ['orig.ident']
ho = hm.run_harmony(X, meta, vars_use=vars_use, theta=2, nclust=100, max_iter_harmony=10)

adata_hvg.obsm['X_pca_harmony'] = ho.Z_corr
print(adata_hvg.obsm['X_pca_harmony'].shape)

sc.pp.neighbors(adata_hvg, use_rep='X_pca_harmony')
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=1.1)


all_same = all(a == b for a, b in zip(adata.obs_names, adata_hvg.obs_names))
print(all_same)
adata.obs["leiden"]=adata_hvg.obs["leiden"]

# Or for more complex hierarchies
results_snap=snap.annotate_hierarchy(
  adata,
  marker_genes,
  group_name="leiden",
  layer="lognorm",
)

tab_assignments = results_snap["assignments"]

tab_assignments.to_csv("deconvolution/assignments_snapseed_pau_040426.txt", sep="\t", index=True)  # index=False prevents writing row numbers


tab_metrics_lv1 = results_snap["metrics"]["level_1"]
tab_metrics_lv2 = results_snap["metrics"]["level_2"]
tab_metrics_lv3 = results_snap["metrics"]["level_3"]
tab_metrics_lv4 = results_snap["metrics"]["level_4"]
tab_metrics_lv5 = results_snap["metrics"]["level_5"]


tab_metrics_lv1.to_csv("deconvolution/metrics_level1_snapseed_pau_040426.txt", sep="\t", index=False)
tab_metrics_lv2.to_csv("deconvolution/metrics_level2_snapseed_pau_040426.txt", sep="\t", index=False)
tab_metrics_lv3.to_csv("deconvolution/metrics_level3_snapseed_pau_040426.txt", sep="\t", index=False)
tab_metrics_lv4.to_csv("deconvolution/metrics_level4_snapseed_pau_040426.txt", sep="\t", index=False)
tab_metrics_lv5.to_csv("deconvolution/metrics_level5_snapseed_pau_040426.txt", sep="\t", index=False)

adata_hvg.obs.to_csv("deconvolution/metadata_scanpy_pau_040426.txt", sep="\t", index=True)

