# shell: pip install git+https://github.com/satijalab/panhumanpy.git

#python:
import panhumanpy as ph
import anndata as ad
import scanpy as sc

# Load my data from annotations.py

# Create Azimuth object
azimuth = ph.AzimuthNN(adata)

#Check predicitions:
cell_metadata = azimuth.cells_meta

#Inspect them:
print(cell_metadata.columns)
print(cell_metadata.head())

#Add them to my data:
adata.obs = adata.obs.join(cell_metadata)
#or
adata.obs["cell_type"] = cell_metadata["predicted.celltype.l2"]

#Visualize annotations:
sc.pl.umap(
    adata, 
    color='azimuth_broad',
    legend_fontsize=14,       # bump this up (default is ~small)
    legend_fontweight='bold',
    legend_loc='right margin', # or 'on data' if you want labels on clusters
    title='Broad Annotations'
)

sc.pl.umap(
    adata, 
    color='azimuth_medium',
    legend_fontsize=14,       # bump this up (default is ~small)
    legend_fontweight='bold',
    legend_loc='right margin', # or 'on data' if you want labels on clusters
    title='Medium Annotations'
)

#More detailed ones:
sc.pl.umap(adata, color='final_level_labels', title='Medium Annotations')

#Confidence score:
sc.pl.umap(adata, color='final_level_confidence', title='Confidence score')




