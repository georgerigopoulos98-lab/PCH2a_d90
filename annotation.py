# Import Libraries
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import harmonypy as hm
import pandas as pd
import numpy as np

# Data retrieval
adata = sc.read_h5ad(r"\\172.23.94.116\AG Mayer$\Students\George\PCH2a_d90\data\raw_data\integrated_data_scvi_d90.h5ad")
adata

#Lizia script:
#scVI (single-cell Variation Interference):
def activate_scvi_embeddings(adata):
    if "X_umap_scVI" not in adata.obsm:
        if "X_umap" in adata.obsm:
            adata.obsm["X_umap_scVI"] = adata.obsm["X_umap"].copy() #If there's no X_umap_scVI, then copy existing X_umap into X_umap_scVI
        else:
            raise KeyError("Integrated object has neither X_umap_scVI nor X_umap. Re-run annotation_integration_scvi.ipynb.")
    adata.obsm["X_umap"] = adata.obsm["X_umap_scVI"].copy() #Make the scVI UMAP the active UMAP
 
    if "X_pca_scVI" not in adata.obsm and "X_pca" in adata.obsm:
        adata.obsm["X_pca_scVI"] = adata.obsm["X_pca"].copy()
    if "X_pca_scVI" in adata.obsm:
        adata.obsm["X_pca"] = adata.obsm["X_pca_scVI"].copy() #Likewise, active PCA becomes the scVI PCA
    return "umap_scVI"

EMBEDDING_BASIS = activate_scvi_embeddings(adata) #so, EMBEDDING_BASIS == "umap_scVI"

#Plot:
sc.pl.embedding(adata, basis=EMBEDDING_BASIS) #sc.pl. is scanpy plot, embedding is a plotting function for any embedding stored in adata.obsm



#AnnDara object metadata:

#obs (observation (cell) annotations, holds per-cell metadata, like sample_id, batch, total_counts etc.):
#'barcode', 'sample', 'sample_id', 'sequencing_id', 'experiment', 'experiment_lane_numbers', 'lane', 'lane_index', 'batch_lane', 'batch', 'batch_id', 'condition', 'cell_line', 'timepoint', 'protocol', 'analysis_group', 'organism_ontology_term_id', 'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id', 'development_stage_ontology_term_id', 'suspension_type', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'outlier_mad', 'manual_outlier', 'doublet_score', 'predicted_doublet', 'doublet_status', 'Score_ER_stress', 'Score_Glycolysis', 'Score_HeatShock', 'Score_OxidativeStress', 'Score_Apoptosis', 'Score_Endoderm', 'Score_Mesoderm', 'Score_Neurogenesis', 'Score_Gliogenesis', 'Score_NSDevelopment', 'Score_ER_stress_robust_z', 'Score_Glycolysis_robust_z', 'Score_HeatShock_robust_z', 'Score_OxidativeStress_robust_z', 'Score_Apoptosis_robust_z', 'Score_Endoderm_robust_z', 'Score_Mesoderm_robust_z', 'lightweight_stress_score', 'off_target_score', 'stressed_cells', 'off_target_tissue', 'leiden_res0_25', 'leiden_res0_5', 'leiden_res1_0', 'leiden_res1_5', 'sepp_author_cell_type', 'sepp_author_cell_type_score', 'sepp_precisest_label', 'sepp_precisest_label_score', 'sepp_subtype', 'sepp_subtype_score', 'sepp_development_stage', 'sepp_development_stage_score', 'sepp_author_stage', 'sepp_author_stage_score', 'sepp_dev_state', 'sepp_dev_state_score', 'aldinger_fig_cell_type', 'aldinger_fig_cell_type_score', 'aldinger_figure_clusters', 'aldinger_figure_clusters_score', 'aldinger_age', 'aldinger_age_score', 'hnoca_level_1_pca', 'hnoca_level_1_pca_score', 'hnoca_level_2_pca', 'hnoca_level_2_pca_score', 'hnoca_level_3_pca', 'hnoca_level_3_pca_score', 'hnoca_level_4_pca', 'hnoca_level_4_pca_score', 'quadrato_final_clusters', 'quadrato_final_clusters_score', 'nano_Class', 'nano_Class_score', 'nano_State', 'nano_State_score', 'nano_Type_v1', 'nano_Type_v1_score', 'nano_Subtype_v1', 'nano_Subtype_v1_score', 'leiden_scVI_res0_25', 'leiden_scVI_res0_5', 'leiden_scVI_res0_75', 'leiden_scVI_res1_0', 'leiden_scVI_res1_5'
adata.obs.head() #Or specific ones:
adata.obs["tissue_type"]

#2 samples: S17818Nr5 and S17818Nr6.
#2 conditions: CTRL and PCH. (all D90 and QB protocol)

#var (variable (gene) annotations, holds per-gene metadata, like gene_ids, highly_variable_rank, etc.):
#'gene_ids', 'feature_types', 'genome', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'highly_variable_scvi_timepoint'
adata.var.head() #Or specific ones:
adata.var["mt"]

#uns (unstructured annotations, holds everything else that doesn't fit into obs or var, like sample_colors, batch_colors, etc.):
#'aldinger_age_colors', 'aldinger_fig_cell_type_colors', 'aldinger_figure_clusters_colors', 'analysis_group_colors', 'assay_ontology_term_id_colors', 'batch_colors', 'batch_id_colors', 'batch_lane_colors', 'cell_line_colors', 'condition_colors', 'development_stage_ontology_term_id_colors', 'doublet_status_colors', 'experiment_colors', 'experiment_lane_numbers_colors', 'hnoca_level_1_pca_colors', 'hnoca_level_2_pca_colors', 'hnoca_level_3_pca_colors', 'hnoca_level_4_pca_colors', 'hvg', 'lane_colors', 'lane_index_colors', 'leiden_res0_25', 'leiden_res0_25_colors', 'leiden_res0_5', 'leiden_res0_5_colors', 'leiden_res1_0', 'leiden_res1_0_colors', 'leiden_res1_5', 'leiden_res1_5_colors', 'leiden_scVI_res0_25', 'leiden_scVI_res0_25_colors', 'leiden_scVI_res0_5', 'leiden_scVI_res0_5_colors', 'leiden_scVI_res0_75', 'leiden_scVI_res0_75_colors', 'leiden_scVI_res1_0', 'leiden_scVI_res1_0_colors', 'leiden_scVI_res1_5', 'leiden_scVI_res1_5_colors', 'log1p', 'manual_outlier_colors', 'nano_Class_colors', 'nano_State_colors', 'nano_Subtype_v1_colors', 'nano_Type_v1_colors', 'neighbors', 'neighbors_scVI', 'off_target_tissue_colors', 'organism_ontology_term_id_colors', 'outlier_mad_colors', 'pca', 'pca_scVI', 'predicted_doublet_colors', 'protocol_colors', 'quadrato_final_clusters_colors', 'sample_colors', 'sample_id_colors', 'scvi_timepoint_integration', 'sepp_author_cell_type_colors', 'sepp_author_stage_colors', 'sepp_dev_state_colors', 'sepp_development_stage_colors', 'sepp_precisest_label_colors', 'sepp_subtype_colors', 'sequencing_id_colors', 'stressed_cells_colors', 'suspension_type_colors', 'timepoint_colors', 'tissue_ontology_term_id_colors', 'tissue_type_colors', 'umap'
adata.uns.keys() #Or specific ones:
adata.uns["pca"]

#obsm (observation multi-dimensional per-cell data, holds stuff like X_umap, etc.): 
#'X_pca', 'X_pca_scVI', 'X_scVI', 'X_scVI_timepoint', 'X_umap', 'X_umap_scVI', 'X_umap_scVI_timepoint'
adata.obsm.keys() #Or specific ones:
adata.obsm["X_umap"]

#varm (variable multi-dimensional per-gene data, holds stuff like PCs): 
#'PCs'
adata.varm.keys() #Or specific ones:
adata.varm["PCs"]

#layers (Alternative version of the main data matrix, preserves multiple parallel versions of the expression matrix without overwriting anything): 
#'counts', 'log1p_norm'
adata.layers.keys() #Or specific ones:
adata.layers["log1p_norm"]

#obsp (observation pairwise matrices (cell-cell relationships), holds stuff like connectivities and distances): 
#'connectivities', 'distances', 'neighbors_scVI_connectivities', 'neighbors_scVI_distances'
adata.obsp.keys() #Or specific ones:
adata.obsp["distances"]


#Check these annotations out:
adata.obs["quadrato_final_clusters"]
adata.obs["quadrato_final_clusters_score"]
adata.obs["nano_Class"]
adata.obs["nano_Class_score"]
adata.obs["nano_State"]
adata.obs["nano_State_score"]
adata.obs["nano_Type_v1"]
adata.obs["nano_Type_v1_score"]
adata.obs["nano_Subtype_v1"]
adata.obs["nano_Subtype_v1_score"]

#And these:
adata.obs["sepp_author_cell_type"]
adata.obs["sepp_author_cell_type_score"]
adata.obs["sepp_precisest_label"]
adata.obs["sepp_precisest_label_score"]
adata.obs["sepp_subtype"]
adata.obs["sepp_subtype_score"]
adata.obs["sepp_development_stage"]
adata.obs["sepp_development_stage_score"]
adata.obs["sepp_author_stage"]
adata.obs["sepp_author_stage_score"]
adata.obs["sepp_dev_state"]
adata.obs["sepp_dev_state_score"]
adata.obs["aldinger_fig_cell_type"]
adata.obs["aldinger_fig_cell_type_score"]
adata.obs["aldinger_figure_clusters"]
adata.obs["aldinger_figure_clusters_score"]
adata.obs["aldinger_age"]
adata.obs["aldinger_age_score"]


sc.pl.umap(
    adata, 
    color='aldinger_age_score',
    legend_fontsize=14,       # bump this up (default is ~small)
    legend_fontweight='bold',
    legend_loc='right margin', # or 'on data' if you want labels on clusters
    title='Aldinger Age Score'
)



#Do leiden clustering for manual annotation.
#Do compositional analysis (stacked barplot) between 2 conditions after finalizing annotations.
#Barplot example from previous workflow:
# Visualize cell type amounts between patients

# Count cell types per donor
ct_counts = (
    adata.obs
    .groupby(["donor_id", "cell_type"])
    .size()
    .unstack(fill_value=0)
)

# Convert counts to proportions per donor
ct_proportions = ct_counts.div(ct_counts.sum(axis=1), axis=0)

# Plot stacked barplot of proportions
ax = ct_proportions.plot(
    kind="bar",
    stacked=True,
    figsize=(10, 6),
    edgecolor="black"
)

plt.title("Cell type composition per donor")
plt.ylabel("Proportion of cells")
plt.xlabel("Donor ID")
plt.ylim(0, 1)  # Proportions go from 0 to 1
plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()



#Background:
#Cell types are further divided in cell states (subtypes).
#Manual annotation or via DGE
#Automated annotation 

#DGE between disease and control
#GO Enrichement analysis/ GSEA



#Questions:
#how many samples do I originally have in total, just 2? 
#same culture conditions, just sifferent cell lines?/ how many organoids sequenced?
#existing annotations?
