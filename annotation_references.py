#shell: pip install scvi-tools

import scvi
import scanpy as sc

ref = sc.read_h5ad("data/references/Sepp_cerebellum.h5ad")
query = sc.read_h5ad(r"\\172.23.94.116\AG Mayer$\Students\George\PCH2a_d90\data\raw_data\integrated_data_scvi_d90.h5ad")

#Combine
ref.obs["dataset"] = "reference"
query.obs["dataset"] = "query"
query.obs["cell_type"] = "Unknown"

adata = ad.concat(
    [ref, query],
    label="dataset",
    keys=["reference", "query"]
)

#Prepare
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key="dataset"
)

#Train
model = scvi.model.SCVI(adata)
model.train()

#Extract latent space
adata.obsm["X_scVI"] = model.get_latent_representation()

#Train SCANVI using reference labels:
scanvi = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key="cell_type",
    unlabeled_category="Unknown"
)

scanvi.train()

#Predict:
predicted = scanvi.predict(
    adata[adata.obs.dataset=="query"]
)

query.obs["cell_type"] = predicted


sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)

#UMAP coordinates are store in adata.obsm["X_umap"]
#Visualize:
sc.pl.umap(adata, color="leiden") #Leiden clusters
sc.pl.umap(adata, color="dataset") #Ref vs query
sc.pl.umap(adata, color="cell_type") #cell type
sc.pl.umap(adata, color=["dataset", "cell_type"])

#Save UMAP:
sc.pl.umap(adata, color="cell_type", save="results/figures/_celltypes.png")

#Save combined dataset:
adata.write("data/processed_data/combined_Sepp_cerebellum.h5ad")

