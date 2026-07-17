#Compare Memento to DGE scanpy (maybe also others in R)
# pip install memento-de

import memento
from statsmodels.stats.multitest import multipletests

#Check my Anndata:
adata.obs["sample"] # S17818Nr6 (CTRL) vs S17818Nr5 (PCH)
adata.obs["condition"] # CTRL vs PCH
adata.layers["counts"] # raw counts in layers
adata.obs["azimuth_broad"] # Cell annotations

# Restore raw counts
adata.X = adata.layers["counts"].copy()

# Parameters
condition_col = "condition"
celltype_col = "azimuth_broad"

control = "CTRL"
case = "PCH"

capture_rate = 0.15      # 10x v2
num_boot = 5000
num_cpus = 8

# DEG thresholds
logfc_thresh = 1.0     # 2-fold change
fdr_thresh = 0.05

# Differential expression
results_all = {}

celltypes = sorted(adata.obs[celltype_col].unique())

for ct in celltypes:

    print(f"Running {ct}")

    subset = adata[
        adata.obs[celltype_col] == ct
    ].copy()

    # Skip tiny populations
    if subset.n_obs < 100:
        print("Skipping (<100 cells)")
        continue

    # Encode condition
    subset.obs["stim"] = (
        subset.obs[condition_col] == case
    ).astype(int)

    # Need cells from both groups
    if subset.obs["stim"].nunique() != 2:
        print("Skipping (only one condition present)")
        continue

    # Run Memento
    res = memento.binary_test_1d(
        adata=subset,
        treatment_col="stim",
        capture_rate=capture_rate,
        num_boot=num_boot,
        num_cpus=num_cpus
    )

    # Multiple testing correction
    res["FDR"] = multipletests(
        res["de_pval"],
        method="fdr_bh"
    )[1]

    res["celltype"] = ct

    # Apply DEG thresholds

    res["significant"] = (
        (res["FDR"] < fdr_thresh) &
        (abs(res["de_coef"]) > logfc_thresh)
    )


    res["direction"] = "Not significant"

    res.loc[
        (res["significant"]) &
        (res["de_coef"] > logfc_thresh),
        "direction"
    ] = "PCH up"


    res.loc[
        (res["significant"]) &
        (res["de_coef"] < -logfc_thresh),
        "direction"
    ] = "CTRL up"

    # Volcano colors

    res["color"] = "gray"

    res.loc[
        res["direction"] == "PCH up",
        "color"
    ] = "red"

    res.loc[
        res["direction"] == "CTRL up",
        "color"
    ] = "blue"

    # Save
    outfile = (
        f"memento_{ct.replace(' ','_')}_DE.csv"
    )

    res.to_csv(outfile, index=False)

    results_all[ct] = res

print("Done.")

#Combine all results:
combined = pd.concat(
    results_all.values(),
    ignore_index=True
)

combined.to_csv(
    "Memento_All_Celltypes.csv",
    index=False
)

#Significant genes:
sig = combined[
    combined["significant"]
].sort_values(
    ["celltype","FDR"]
)

sig.to_csv(
    "Memento_significant_DEGs.csv",
    index=False
)


print(
    f"Total significant DEGs: {sig.shape[0]}"
)

#Visualize results
# >0 means gene is upregulated in PCH relative to CTRL
# <0 means gene is downregulated in PCH
#Volcano plot:
#Label top genes:
for ct, res in results_all.items():

    plt.figure(figsize=(8,6))


    plt.scatter(
        res["de_coef"],
        -np.log10(res["FDR"] + 1e-300),
        c=res["color"],
        alpha=0.5
    )


    # Threshold lines

    plt.axhline(
        -np.log10(fdr_thresh),
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axvline(
        logfc_thresh,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axvline(
        -logfc_thresh,
        linestyle="--",
        color="black",
        linewidth=1
    )


    # Label top genes

    top = (
        res[
            res["significant"]
        ]
        .sort_values("FDR")
        .head(15)
    )


    for _, row in top.iterrows():

        plt.text(
            row["de_coef"],
            -np.log10(row["FDR"] + 1e-300),
            row["gene"],
            fontsize=8
        )


    plt.xlabel(
        "Log2 fold change (PCH vs CTRL)"
    )

    plt.ylabel(
        "-log10(FDR)"
    )

    plt.title(
        f"Memento DE: {ct} (PCH vs CTRL)"
    )


    plt.tight_layout()

    plt.savefig(
        f"Volcano_Memento_{ct.replace(' ','_')}.png",
        dpi=300
    )

    plt.show()


#Top DE genes barplots:
for ct, res in results_all.items():

    top = (
        res[
            res["significant"]
        ]
        .sort_values("FDR")
        .head(20)
        .sort_values("de_coef")
    )


    if top.empty:
        continue


    plt.figure(figsize=(8,6))


    plt.barh(
        top["gene"],
        top["de_coef"]
    )


    plt.xlabel(
        "Log2 fold change (PCH vs CTRL)"
    )

    plt.title(
        f"Top Memento DE genes: {ct}"
    )

    plt.tight_layout()

    plt.show()




#STANDARD SCANPY WORKFLOW FOR DGE:
# Show differentially expressed markers per patient

# Make sure donor_id is categorical
adata.obs['condition'] = adata.obs['condition'].astype('category')

# Perform differential expression analysis
# Compare using Wilcoxon rank-sum test
sc.tl.rank_genes_groups(
    adata,
    groupby='condition',
    groups=['PCH'],
    reference='CTRL',
    method='wilcoxon'
)

# Get DE results as a pandas dataframe
de_results = sc.get.rank_genes_groups_df(adata, group='PCH')

# Add marker names from adata.var['feature_name']
de_results['marker_name'] = de_results['names'].map(
    lambda x: adata.var.loc[x, 'feature_name'] if x in adata.var.index else x
)

# Calculate -log10(p-value)
de_results['neg_log10_pval'] = -np.log10(de_results['pvals'] + 1e-300)

# Thresholds
logfc_thresh = 1.0    #should be 1! not less, also don't do it across all cell types since different proportions can affect DEG.
pval_thresh = 0.05     #always 0.05 or less!

# Color by up/down regulation
de_results["neg_log10_padj"] = -np.log10(de_results["pvals_adj"] + 1e-300)

def volcano_color(row):
    if row["pvals_adj"] < 0.05 and row["logfoldchanges"] > 1:
        return "red"
    elif row["pvals_adj"] < 0.05 and row["logfoldchanges"] < -1:
        return "blue"
    else:
        return "gray"

de_results['color'] = de_results.apply(volcano_color, axis=1)

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(
    de_results['logfoldchanges'],
    de_results['neg_log10_pval'],
    c=de_results['color'],
    alpha=0.7,
    edgecolor='k'
)

# Threshold lines
plt.axhline(-np.log10(pval_thresh), color='black', linestyle='--', lw=1)
plt.axvline(logfc_thresh, color='black', linestyle='--', lw=1)
plt.axvline(-logfc_thresh, color='black', linestyle='--', lw=1)

# Labels and title
plt.xlabel('Log2 Fold Change (positive = PCH up, negative = CTRL up)')
plt.ylabel('-log10(p-value)')
plt.title('Volcano Plot of Differentially Expressed Markers')

# Legend
import matplotlib.patches as mpatches
red_patch = mpatches.Patch(color='red', label='Up in PCH')
blue_patch = mpatches.Patch(color='blue', label='Up in CTRL')
gray_patch = mpatches.Patch(color='gray', label='Not significant')
plt.legend(handles=[red_patch, blue_patch, gray_patch], loc='upper right')

# Label top markers
top_markers = de_results.sort_values('neg_log10_pval', ascending=False).head(10)
for _, row in top_markers.iterrows():
    plt.text(row['logfoldchanges'], row['neg_log10_pval'], row['names'],
             fontsize=8, ha='right' if row['logfoldchanges'] < 0 else 'left')

plt.tight_layout()
plt.show()





#Individual volcano plots per cell type:
# Parameters
condition_col = "condition"
celltype_col = "azimuth_broad"

control = "CTRL"
case = "PCH"

logfc_thresh = 1.0
padj_thresh = 0.05


# Make condition categorical
adata.obs[condition_col] = adata.obs[condition_col].astype("category")


# Store results
results_all = {}


# Get cell types
celltypes = sorted(
    adata.obs[celltype_col].dropna().unique()
)


# ============================
# Run DE per cell type
# ============================

for ct in celltypes:

    print(f"\nRunning DE for: {ct}")


    # Subset cell type
    subset = adata[
        adata.obs[celltype_col] == ct
    ].copy()


    # Skip small populations
    if subset.n_obs < 100:
        print("Skipping (<100 cells)")
        continue


    # Check both conditions exist

    if (
        control not in subset.obs[condition_col].values
        or
        case not in subset.obs[condition_col].values
    ):
        print("Skipping (missing condition)")
        continue



    # ============================
    # Wilcoxon DE
    # ============================

    sc.tl.rank_genes_groups(
        subset,
        groupby=condition_col,
        groups=[case],
        reference=control,
        method="wilcoxon"
    )


    # Extract results

    de_results = sc.get.rank_genes_groups_df(
        subset,
        group=case
    )


    # Add marker names

    if "feature_name" in subset.var.columns:

        de_results["marker_name"] = de_results["names"].map(
            lambda x: subset.var.loc[x, "feature_name"]
            if x in subset.var.index else x
        )

    else:

        de_results["marker_name"] = de_results["names"]



    # Calculate significance

    de_results["neg_log10_padj"] = (
        -np.log10(
            de_results["pvals_adj"] + 1e-300
        )
    )


    # ============================
    # Classification
    # ============================

    def volcano_color(row):

        if (
            row["pvals_adj"] < padj_thresh
            and row["logfoldchanges"] > logfc_thresh
        ):
            return "red"

        elif (
            row["pvals_adj"] < padj_thresh
            and row["logfoldchanges"] < -logfc_thresh
        ):
            return "blue"

        else:
            return "gray"


    de_results["color"] = de_results.apply(
        volcano_color,
        axis=1
    )


    de_results["celltype"] = ct


    results_all[ct] = de_results



    # Save individual results

    de_results.to_csv(
        f"Scanpy_DE_{ct.replace(' ','_')}.csv",
        index=False
    )



    # ============================
    # Volcano plot
    # ============================

    plt.figure(figsize=(8,6))


    plt.scatter(
        de_results["logfoldchanges"],
        de_results["neg_log10_padj"],
        c=de_results["color"],
        alpha=0.6,
        edgecolor="k"
    )


    # Thresholds

    plt.axhline(
        -np.log10(padj_thresh),
        color="black",
        linestyle="--",
        lw=1
    )

    plt.axvline(
        logfc_thresh,
        color="black",
        linestyle="--",
        lw=1
    )

    plt.axvline(
        -logfc_thresh,
        color="black",
        linestyle="--",
        lw=1
    )


    # Label top genes

    top_markers = (
        de_results[
            (de_results["pvals_adj"] < padj_thresh)
            &
            (abs(de_results["logfoldchanges"]) > logfc_thresh)
        ]
        .sort_values(
            "pvals_adj"
        )
        .head(10)
    )


    for _, row in top_markers.iterrows():

        plt.text(
            row["logfoldchanges"],
            row["neg_log10_padj"],
            row["marker_name"],
            fontsize=8,
            ha="right" if row["logfoldchanges"] < 0 else "left"
        )


    # Labels

    plt.xlabel(
        "Log2 Fold Change (PCH vs CTRL)"
    )

    plt.ylabel(
        "-log10(adjusted p-value)"
    )


    plt.title(
        f"Volcano Plot: {ct}\nPCH vs CTRL"
    )


    # Legend

    red_patch = mpatches.Patch(
        color="red",
        label="Up in PCH"
    )

    blue_patch = mpatches.Patch(
        color="blue",
        label="Up in CTRL"
    )

    gray_patch = mpatches.Patch(
        color="gray",
        label="Not significant"
    )


    plt.legend(
        handles=[
            red_patch,
            blue_patch,
            gray_patch
        ],
        loc="upper right"
    )


    plt.tight_layout()


    plt.savefig(
        f"Volcano_Scanpy_{ct.replace(' ','_')}.png",
        dpi=300
    )


    plt.show()