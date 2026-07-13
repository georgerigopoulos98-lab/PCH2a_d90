# Import Libraries
import scanpy as sc
import matplotlib.pyplot as plt
import harmonypy as hm
import pandas as pd
import numpy as np

# Data retrieval
adata = sc.read_h5ad(r"\\172.23.94.116\AG Mayer$\Students\George\PCH2a_d90\data\raw_data\integrated_data_scvi_d90.h5ad")
adata

#they are under Featurename instead of names???? (old info from me)

#AnnDara object metadata:

#obs (observation (cell) annotations, holds per-cell metadata, like sample_id, batch, total_counts etc.):
#'barcode', 'sample', 'sample_id', 'sequencing_id', 'experiment', 'experiment_lane_numbers', 'lane', 'lane_index', 'batch_lane', 'batch', 'batch_id', 'condition', 'cell_line', 'timepoint', 'protocol', 'analysis_group', 'organism_ontology_term_id', 'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id', 'development_stage_ontology_term_id', 'suspension_type', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'outlier_mad', 'manual_outlier', 'doublet_score', 'predicted_doublet', 'doublet_status', 'Score_ER_stress', 'Score_Glycolysis', 'Score_HeatShock', 'Score_OxidativeStress', 'Score_Apoptosis', 'Score_Endoderm', 'Score_Mesoderm', 'Score_Neurogenesis', 'Score_Gliogenesis', 'Score_NSDevelopment', 'Score_ER_stress_robust_z', 'Score_Glycolysis_robust_z', 'Score_HeatShock_robust_z', 'Score_OxidativeStress_robust_z', 'Score_Apoptosis_robust_z', 'Score_Endoderm_robust_z', 'Score_Mesoderm_robust_z', 'lightweight_stress_score', 'off_target_score', 'stressed_cells', 'off_target_tissue', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1_0', 'leiden_res1_5', 'sepp_author_cell_type', 'sepp_author_cell_type_score', 'sepp_precisest_label', 'sepp_precisest_label_score', 'sepp_subtype', 'sepp_subtype_score', 'sepp_development_stage', 'sepp_development_stage_score', 'sepp_author_stage', 'sepp_author_stage_score', 'sepp_dev_state', 'sepp_dev_state_score', 'aldinger_fig_cell_type', 'aldinger_fig_cell_type_score', 'aldinger_figure_clusters', 'aldinger_figure_clusters_score', 'aldinger_age', 'aldinger_age_score', 'hnoca_level_1_pca', 'hnoca_level_1_pca_score', 'hnoca_level_2_pca', 'hnoca_level_2_pca_score', 'hnoca_level_3_pca', 'hnoca_level_3_pca_score', 'hnoca_level_4_pca', 'hnoca_level_4_pca_score', 'quadrato_final_clusters', 'quadrato_final_clusters_score', 'nano_Class', 'nano_Class_score', 'nano_State', 'nano_State_score', 'nano_Type_v1', 'nano_Type_v1_score', 'nano_Subtype_v1', 'nano_Subtype_v1_score', 'leiden_scVI_res0_25', 'leiden_scVI_res0_5', 'leiden_scVI_res0_75', 'leiden_scVI_res1_0', 'leiden_scVI_res1_5'
adata.obs.head() #Or specific ones:
adata.obs["barcode"]

#var (variable (gene) annotations, holds per-gene metadata, like gene_ids, highly_variable_rank, etc.):
#'gene_ids', 'feature_types', 'genome', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'highly_variable_scvi_timepoint'
adata.var.head() #Or specific ones:
adata.var["feature_types"]

#uns (unstructured annotations, holds everything else that doesn't fit into obs or var, like sample_colors, batch_colors, etc.):
#'aldinger_age_colors', 'aldinger_fig_cell_type_colors', 'aldinger_figure_clusters_colors', 'analysis_group_colors', 'assay_ontology_term_id_colors', 'batch_colors', 'batch_id_colors', 'batch_lane_colors', 'cell_line_colors', 'condition_colors', 'development_stage_ontology_term_id_colors', 'doublet_status_colors', 'experiment_colors', 'experiment_lane_numbers_colors', 'hnoca_level_1_pca_colors', 'hnoca_level_2_pca_colors', 'hnoca_level_3_pca_colors', 'hnoca_level_4_pca_colors', 'hvg', 'lane_colors', 'lane_index_colors', 'leiden_res0_25', 'leiden_res0_25_colors', 'leiden_res0_5', 'leiden_res0_5_colors', 'leiden_res1_0', 'leiden_res1_0_colors', 'leiden_res1_5', 'leiden_res1_5_colors', 'leiden_scVI_res0_25', 'leiden_scVI_res0_25_colors', 'leiden_scVI_res0_5', 'leiden_scVI_res0_5_colors', 'leiden_scVI_res0_75', 'leiden_scVI_res0_75_colors', 'leiden_scVI_res1_0', 'leiden_scVI_res1_0_colors', 'leiden_scVI_res1_5', 'leiden_scVI_res1_5_colors', 'log1p', 'manual_outlier_colors', 'nano_Class_colors', 'nano_State_colors', 'nano_Subtype_v1_colors', 'nano_Type_v1_colors', 'neighbors', 'neighbors_scVI', 'off_target_tissue_colors', 'organism_ontology_term_id_colors', 'outlier_mad_colors', 'pca', 'pca_scVI', 'predicted_doublet_colors', 'protocol_colors', 'quadrato_final_clusters_colors', 'sample_colors', 'sample_id_colors', 'scvi_timepoint_integration', 'sepp_author_cell_type_colors', 'sepp_author_stage_colors', 'sepp_dev_state_colors', 'sepp_development_stage_colors', 'sepp_precisest_label_colors', 'sepp_subtype_colors', 'sequencing_id_colors', 'stressed_cells_colors', 'suspension_type_colors', 'timepoint_colors', 'tissue_ontology_term_id_colors', 'tissue_type_colors', 'umap'
adata.uns.keys() #Or specific ones:
adata.uns["timepoint_colors"]

#obsm (observation multi-dimensional per-cell data, holds stuff like X_umap, etc.): 
#'X_pca', 'X_pca_scVI', 'X_scVI', 'X_scVI_timepoint', 'X_umap', 'X_umap_scVI', 'X_umap_scVI_timepoint'
adata.obsm.keys() #Or specific ones:
adata.obsm["X_scVI"]

#varm (variable multi-dimensional per-gene data, holds stuff like PCs): 
#'PCs'
adata.varm.keys() #Or specific ones:
adata.varm["PCs"]

#layers (Alternative version of the main data matrix, preserves multiple parallel versions of the expression matrix without overwriting anything): 
#'counts', 'log1p_norm'
adata.layers.keys() #Or specific ones:
adata.layers["counts"]

#obsp (observation pairwise matrices (cell-cell relationships), holds stuff like connectivities and distances): 
#'connectivities', 'distances', 'neighbors_scVI_connectivities', 'neighbors_scVI_distances'
adata.obsp.keys() #Or specific ones:
adata.obsp["connectivities"]


#Lizia script:
#scVI (single-cell Variation Interference) (batch effect correction):
def activate_scvi_embeddings(adata):
    if "X_umap_scVI" not in adata.obsm:
        if "X_umap" in adata.obsm:
            adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy()
        else:
            raise KeyError("Integrated object has neither X_umap_scVI nor X_umap. Re-run annotation_integration_scvi.ipynb.")
    adata.obsm["X_umap"] = adata.obsm["X_umap_scVI"].copy()
 
    if "X_pca_scVI" not in adata.obsm and "X_pca" in adata.obsm:
        adata.obsm["X_pca_scVI"] = adata.obsm["X_pca"].copy()
    if "X_pca_scVI" in adata.obsm:
        adata.obsm["X_pca"] = adata.obsm["X_pca_scVI"].copy()
    return "umap_scVI"

EMBEDDING_BASIS = activate_scvi_embeddings(adata)



# Display Annotations

# Check relevant keys in adata

print("adata.obs keys:", adata.obs_keys())
print("adata.var keys:", adata.var_keys())
print("adata.uns keys:", adata.uns_keys())

adata.obs['suspension_type'].head() # Suspensio_type: All are nucleus.
adata.obs['sex'].head() # All are female. 
adata.obs['tissue'].head() # All are primary motor cortex samples.
adata.obs['development_stage'].head() # Most likely just 2 patients, samples should be from post-mortem adult brains.

adata.obs['assay_ontology_term_id'].head() # (2, object): ['EFO:0009922', 'EFO:0030059']
adata.uns['donor_id_colors'] #array(['#1f77b4', '#ff7f0e'], dtype=object)


# Visualize annotations

sc.pl.umap(adata, color='PrimaryAnnotation', title='Primary Annotations')
sc.pl.umap(adata, color='ThirdAnnotation', title='Third Annotations')
sc.pl.umap(adata, color='cell_type', title='Cell Types')

sc.pl.umap(adata, color='donor_id', title='ALS vs Healthy') # Can also use 'disease' instead of 'donor_id'
sc.pl.umap(adata, color='assay', title='Assay Type') # Different assay was used for each donor type. The ALS ones were analyzed with 10x Multiome which also includes ATAC-seq.
sc.pl.umap(adata, color='development_stage', title='Development stage') # 80 year old with ALS and 73 year old healthy control.



#Background:
#Cell types are further divided in cell states (subtypes).
#Manual annotation or via DGE
#Automated annotation 

#DGE between disease and control
#GO Enrichement analysis/ GSEA
