# Import Libraries
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Data retrieval
adata = sc.read_h5ad(#Path here)
adata

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
adata.obs["total_counts"]

#Cell line: BIONi010-C HOM and WT
#2 samples: S17818Nr5 and S17818Nr6.
#2 conditions: CTRL and PCH. (all D90 and QB protocol (quadrato))

#var (variable (gene) annotations, holds per-gene metadata, like gene_ids, highly_variable_rank, etc.):
#'gene_ids', 'feature_types', 'genome', 'mt', 'ribo', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'highly_variable_scvi_timepoint'
adata.var.head() #Or specific ones:
adata.var["total_counts"]

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


# Plot conditions:
sc.pl.umap(
    adata, 
    color='condition',
    size=20,
    legend_fontsize=14,       # bump this up (default is ~small)
    legend_fontweight='bold',
    legend_loc='right margin', # or 'on data' if you want labels on clusters
    title='Conditions'
)

#Plot conditions separately:

conditions = adata.obs['condition'].unique()

fig, axes = plt.subplots(
    1,
    len(conditions),
    figsize=(7 * len(conditions), 6)
)

for ax, condition in zip(axes, conditions):
    subset = adata[adata.obs['condition'] == condition].copy()

    sc.pl.umap(
        subset,
        color='condition',
        ax=ax,
        size=20,
        title=condition,
        show=False,
        legend_loc='none'
    )

    ax.set_title(
        condition,
        fontsize=20,
        fontweight='bold'
    )

plt.tight_layout()
plt.show()






#Clustering

#Compare leiden resolutions:
sc.pl.embedding(
    adata,
    basis="X_umap_scVI",
    color=[
        "leiden_scVI_res0_25",
        "leiden_scVI_res0_5",
        "leiden_scVI_res0_75",
        "leiden_scVI_res1_0",
        "leiden_scVI_res1_5"
    ],
    ncols=2,
    legend_loc="on data",
    frameon=False
)

#Use resolution 0.5 SCVI (n_neighbours 20 was used for these).

#Use known cell type markers for manual annotation:

sc.pl.umap(
    adata,
    color=[
        "PAX6",
        "NFIA",
        "NFIB",
        "FOXG1"
    ]
)
                                                                                   
#Check highly expressed markers for each cluster:

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_scVI_res0_5",
    method="wilcoxon",
    key_added="rank_genes_scVI_res0_5"
)

sc.pl.rank_genes_groups(
    adata,
    key="rank_genes_scVI_res0_5",
    n_genes=10,
    sharey=False
)





#Manual annotations:

#0.5 Clusters:
#0: Radial Glia (VIM) to Astrocytes from left to right (right side mature) (CLU)
#1: Radial Glia (VIM) to Astrocytes from left to right (right side mature) (CLU)
#2: Neurons
#3: GABAergic (inhibitory) neurons (EBF3)
#4: GABAergic interneurons (FOXP2) (NXPH2)
#5: Purkinje cells (lower half) (SKOR2) (GRID2)
#6: Oligodendrocyte precursor cells (OLIG1) (OLIG2) and top right specificcaly dividing! (TOP2A)
#7: Cerebellar granule neurons (Glutamatergic) (FGF12)
#8: Glutamatergic (excitatory) neurons bottom half (NEUROD6) (SLC17A6) AND Gabaergic (inhibitory) neurons top half (NR2F2) (SOX14)
#9: GABAergic interneurons (FOXP2)
#10:
#11: 
#12: GABAergic interneurons (CALB2) and Purkinje cells (upper half specifically) (SKOR2)
#13:
#14:
#15: Serotonergic neurons (TPH2)
#16: Fibroblasts (COL3A1)




#All cluster right: radial glia / progenitors (SOX2)

#All cluster left side: Neurons

#Nano phases:
#ALL the clusters on the right: Non-dividing
#In between: Dividing
#ALL in left: Post mitotic

#All between 17-20PCW according to Aldinger and Sepp

sc.pl.embedding(
    adata,
    basis="X_umap_scVI",
    color=[
        "leiden_scVI_res1_5"
    ],
    ncols=2,
    legend_loc="on data",
    frameon=False
)

#Name clusters manually (use 1.5 resolution for more specificity):

#General annotations: eCN (excitatory cerebellar nuclei), iCN/PC (inhibitory cerebellar nuclei/ Purkinje cells), oligodendrocyte precursor cell (OPC) and Progenitor/ Astrocytes.

#General clusters:
cluster_to_celltype = {
    '0': 'Progenitor/ Astrocytes',
    '1': 'Progenitor/ Astrocytes',
    '2': 'iCN/PC',
    '3': 'OPC',
    '4': 'iCN/PC',
    '5': 'Progenitor/ Astrocytes',
    '6': 'iCN/PC',
    '7': 'Progenitor/ Astrocytes',
    '8': 'iCN/PC',
    '9': 'iCN/PC',
    '10': 'iCN/PC',
    '11': 'iCN/PC',
    '12': 'iCN/PC',
    '13': 'iCN/PC',
    '14': 'iCN/PC',
    '15': 'iCN/PC',
    '16': 'iCN/PC',
    '17': 'iCN/PC',
    '18': 'iCN/PC',
    '19': 'iCN/PC',
    '20': 'iCN/PC',
    '21': 'iCN/PC',
    '22': 'OPC - dividing',
    '23': 'iCN/PC',
    '24': 'iCN/PC',
    '25': 'iCN/PC',
    '26': 'iCN/PC',
    '27': 'Serotonergic neurons',
    '28': 'Progenitor/ Astrocytes',
    '29': 'Fibroblasts',
    '30': 'iCN/PC'
}

#Add to metadata
adata.obs['cell_type'] = (
    adata.obs['leiden_scVI_res1_5']
    .map(cluster_to_celltype)
)

#Specific subtypes:
cluster_to_subtype = {
    '0': 'Radial Glia (VIM) to Astrocytes (CLU) (left to right)',
    '1': 'Radial Glia (VIM) to Astrocytes (CLU) (left to right)',
    '2': 'GABAergic interneurons (FOXP2/NXPH2)',
    '3': 'Oligodendrocyte precursor cells (OLIG1/OLIG2)',
    '4': 'Purkinje cells (SKOR2/GRID2)',
    '5': 'Radial Glia (VIM) to Astrocytes (CLU) (left to right)',
    '6': 'Ventricular zone (VZ) neuroblasts (NFIA/NFIB)',
    '7': 'Radial Glia (VIM) to Astrocytes (CLU) (left to right)',
    '8': 'iCN/PC (GAD1/GAD2)',
    '9': 'iCN/PC (GAD1/GAD2)',
    '10': 'iCN/PC (GAD1/GAD2)',
    '11': 'Ventricular zone (VZ) neuroblasts (NFIA/NFIB)',
    '12': 'Glutamatergic neurons bottom half (NEUROD6/SLC17A6) /\nGABAergic neurons top half (NR2F2/SOX14)',
    '13': 'iCN/PC (GAD1/GAD2)',
    '14': 'iCN/PC (GAD1/GAD2)',
    '15': 'iCN/PC (GAD1/GAD2)',
    '16': 'Ventricular zone (VZ) neuroblasts (NFIA/NFIB)',
    '17': 'Ventricular zone (VZ) neuroblasts (NFIA/NFIB)',
    '18': 'GABAergic interneurons (FOXP2/NXPH2)',
    '19': 'GABAergic interneurons (FOXP2)',
    '20': 'iCN/PC (GAD1/GAD2)',
    '21': 'iCN/PC (GAD1/GAD2)',
    '22': 'Oligodendrocyte precursor cells - divinding (TOP2A)',
    '23': 'GABAergic interneurons (CALB2)',
    '24': 'iCN/PC (GAD1/GAD2)',
    '25': 'iCN/PC (GAD1/GAD2)',
    '26': 'Purkinje cells (SKOR2)',
    '27': 'Serotonergic neurons (TPH2)',
    '28': 'Radial Glia (VIM) to Astrocytes (CLU) (left to right)',
    '29': 'Fibroblasts (COL3A1)',
    '30': 'iCN/PC (GAD1/GAD2)'
}

#Add to metadata
adata.obs['cell_subtype'] = (
    adata.obs['leiden_scVI_res1_5']
    .map(cluster_to_subtype)
)

#Visualize final annotations
sc.pl.umap(
    adata,
    color='cell_subtype',
    legend_fontsize=14,
    legend_fontweight='bold',
    legend_loc='right margin',
    size=20,
    title='Subtypes'
)

sc.pl.umap(
    adata,
    color='cell_type',
    legend_fontsize=14,
    legend_fontweight='bold',
    legend_loc='right margin',
    size=20,
    title='General annotations'
)


#Cluster 16 I previously annotated as: Cerebellar granule neurons (Glutamatergic) (FGF12) OR off-target forebrain FOXG1(?).
#Also clusters 6 and 11 as: GABAergic neurons (EBF3)


#Check amount of cells in each annotated cluster:

adata.obs["cell_type"].value_counts()
adata.obs["cell_subtype"].value_counts()

#And between each condition:

pd.crosstab(
    adata.obs["cell_type"],
    adata.obs["condition"]
)

pd.crosstab(
    adata.obs["cell_subtype"],
    adata.obs["condition"]
)


#Compositional analysis (stacked barplot) between 2 conditions:
# Count cell types per donor
ct_counts = (
    adata.obs
    .groupby(["condition", "cell_type"])
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

plt.title("Cell type composition")
plt.ylabel("Proportion of cells")
plt.xlabel("Condition")
plt.ylim(0, 1)  # Proportions go from 0 to 1
plt.legend(fontsize=14, title_fontsize=14, title="Cell type", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()


#Same for subtypes:
ct_counts = (
    adata.obs
    .groupby(["condition", "cell_subtype"])
    .size()
    .unstack(fill_value=0)
)

ct_proportions = ct_counts.div(ct_counts.sum(axis=1), axis=0)

ax = ct_proportions.plot(
    kind="bar",
    stacked=True,
    figsize=(10, 6),
    edgecolor="black"
)

plt.title("Cell subtype composition")
plt.ylabel("Proportion of cells")
plt.xlabel("Condition")
plt.ylim(0, 1)  # Proportions go from 0 to 1
plt.legend(fontsize=12, title_fontsize=12, title="Cell subtype", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()


#Save annotated dataset and do pseudobulk DGE in R:
adata.write("data/processed_data/PCH2a_d90_annotated.h5ad")



