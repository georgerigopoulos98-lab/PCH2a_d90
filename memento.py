#Compare Memento to DGE scanpy (maybe also others in R)
# pip install memento-de

import memento
from statsmodels.stats.multitest import multipletests

#Import annotated dataset:
adata = sc.read_h5ad("data/processed_data/PCH2a_d90_annotated.h5ad")
adata

#Check my Anndata:
adata.obs["sample"] # S17818Nr6 (CTRL) vs S17818Nr5 (PCH)
adata.obs["condition"] # CTRL vs PCH
adata.layers["counts"] # raw counts in layers
adata.obs["cell_type"] # Cell annotations

# Restore raw counts
adata.X = adata.layers["counts"].copy()

# Parameters
condition_col = "condition"
celltype_col = "cell_type"

control = "CTRL"
case = "PCH"

# Memento capture/detection parameter q
# Currently set to 0.15 based on the initial analysis setup:
capture_rate = 0.15
num_boot = 5000
num_cpus = 8

# DEG thresholds
lnfc_thresh = np.log(2)   # 2-fold change on natural-log scale
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
    # CTRL = 0
    # PCH = 1
    subset.obs["stim"] = (
        subset.obs[condition_col] == case
    ).astype(int)

    # Need cells from both groups
    if subset.obs["stim"].nunique() != 2:
        print("Skipping (only one condition present)")
        continue

    # Run Memento
    # This tests BOTH:
    # 1. Differential mean expression
    # 2. Differential variability
    res = memento.binary_test_1d(
        adata=subset,
        treatment_col="stim",
        capture_rate=capture_rate,
        num_boot=num_boot,
        num_cpus=num_cpus
    )

    # Convert Memento natural-log coefficients

    # Mean expression:
    # de_coef is on the natural-log fold-change scale
    res["mean_fold_change"] = np.exp(res["de_coef"])

    # Convert to log2 fold change
    res["log2FC"] = (
        res["de_coef"] / np.log(2)
    )

    # Variability:
    # dv_coef is on the natural-log variance fold-change scale
    res["variance_fold_change"] = np.exp(res["dv_coef"])

    # Convert to log2 variance fold change
    res["log2_variance_FC"] = (
        res["dv_coef"] / np.log(2)
    )

    # Multiple testing correction

    # Mean expression FDR
    res["DE_FDR"] = multipletests(
        res["de_pval"],
        method="fdr_bh"
    )[1]

    # Variability FDR
    res["DV_FDR"] = multipletests(
            res["dv_pval"],
            method="fdr_bh"
        )[1]

    res["celltype"] = ct

    # --------------------------------------------------
    # Significant differential mean expression
    # --------------------------------------------------

    res["DE_significant"] = (
        (res["DE_FDR"] < fdr_thresh) &
        (abs(res["de_coef"]) >= lnfc_thresh)
    )

    # --------------------------------------------------
    # Significant differential variability
    # --------------------------------------------------

    res["DV_significant"] = (
        (res["DV_FDR"] < fdr_thresh) &
        (abs(res["dv_coef"]) > lnfc_thresh)
    )

    # --------------------------------------------------
    # Classify genes based on mean and variability
    # --------------------------------------------------

    res["pattern"] = "Neither"

    res.loc[
        (res["DE_significant"]) &
        (~res["DV_significant"]),
        "pattern"
    ] = "Mean only"

    res.loc[
        (~res["DE_significant"]) &
        (res["DV_significant"]),
        "pattern"
    ] = "Variability only"

    res.loc[
        (res["DE_significant"]) &
        (res["DV_significant"]),
        "pattern"
    ] = "Mean + variability"

    res["DE_direction"] = "Not significant"

    res.loc[
        (res["DE_significant"]) &
        (res["de_coef"] > lnfc_thresh),
        "DE_direction"
    ] = "PCH up"

    res.loc[
        (res["DE_significant"]) &
        (res["de_coef"] < -lnfc_thresh),
        "DE_direction"
    ] = "CTRL up"

    res["DV_direction"] = "Not significant"

    res.loc[
        (res["DV_significant"]) &
        (res["dv_coef"] > lnfc_thresh),
        "DV_direction"
    ] = "More variable in PCH"

    res.loc[
        (res["DV_significant"]) &
        (res["dv_coef"] < -lnfc_thresh),
        "DV_direction"
    ] = "More variable in CTRL"


    # Store results
    results_all[ct] = res

    # Print summary
    print(
        f"  Significant mean-expression genes: "
        f"{res['DE_significant'].sum()}"
    )

    print(
        f"  Significant variability genes: "
        f"{res['DV_significant'].sum()}"
    )

    print(
        res["pattern"].value_counts()
    )

print("Done.")

#Combine all results:
combined = pd.concat(
    results_all.values(),
    ignore_index=True
)


#Significant mean-expression genes
sig_DE = combined[
    combined["DE_significant"]
].sort_values(
    ["celltype","DE_FDR"]
)


print(
    f"Total significant differential mean genes: "
    f"{sig_DE.shape[0]}"
)

#Significant variability genes
sig_DV = combined[
    combined["DV_significant"]
].sort_values(
    ["celltype", "DV_FDR"]
)

print(
    f"Total significant differential variability genes: "
    f"{sig_DV.shape[0]}"
)

# Save results

combined.to_csv(
    "Memento_all_results.csv",
    index=False
)

sig_DE.to_csv(
    "Memento_significant_mean_genes.csv",
    index=False
)

sig_DV.to_csv(
    "Memento_significant_variability_genes.csv",
    index=False
)




#Visualize results
# >0 means gene is upregulated in PCH relative to CTRL
# <0 means gene is downregulated in PCH


# ==========================================================
# MEAN EXPRESSION VOLCANO PLOTS
# ==========================================================

for ct, res in results_all.items():

    plt.figure(figsize=(8, 6))

    # Plot all genes
    plt.scatter(
        res["log2FC"],
        -np.log10(res["DE_FDR"] + 1e-300),
        alpha=0.5
    )

    # 2-fold effect-size thresholds
    plt.axvline(
        1,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axvline(
        -1,
        linestyle="--",
        color="black",
        linewidth=1
    )

    # FDR threshold
    plt.axhline(
        -np.log10(fdr_thresh),
        linestyle="--",
        color="black",
        linewidth=1
    )

    # Label significant genes
    top = (
        res[
            res["DE_significant"]
        ]
        .sort_values("DE_FDR")
        .head(15)
    )

    for _, row in top.iterrows():

        plt.text(
            row["log2FC"],
            -np.log10(row["DE_FDR"] + 1e-300),
            row["gene"],
            fontsize=8
        )

    plt.xlabel(
        "log2 fold change (PCH vs CTRL)"
    )

    plt.ylabel(
        "-log10(FDR)"
    )

    plt.title(
        f"Memento differential mean expression: {ct}"
    )

    plt.tight_layout()

    plt.savefig(
        f"Volcano_Memento_Mean_{ct.replace(' ', '_').replace('/', '_')}.png",
        dpi=300
    )

    plt.show()


# ==========================================================
# DIFFERENTIAL VARIABILITY VOLCANO PLOTS
# ==========================================================

for ct, res in results_all.items():

    plt.figure(figsize=(8, 6))

    # Plot all genes
    plt.scatter(
        res["log2_variance_FC"],
        -np.log10(res["DV_FDR"] + 1e-300),
        alpha=0.5
    )

    # 2-fold variance effect-size thresholds
    plt.axvline(
        1,
        linestyle="--",
        color="black",
        linewidth=1
    )

    plt.axvline(
        -1,
        linestyle="--",
        color="black",
        linewidth=1
    )

    # FDR threshold
    plt.axhline(
        -np.log10(fdr_thresh),
        linestyle="--",
        color="black",
        linewidth=1
    )

    # Label significant genes
    top = (
        res[
            res["DV_significant"]
        ]
        .sort_values("DV_FDR")
        .head(15)
    )

    for _, row in top.iterrows():

        plt.text(
            row["log2_variance_FC"],
            -np.log10(row["DV_FDR"] + 1e-300),
            row["gene"],
            fontsize=8
        )

    plt.xlabel(
        "log2 variance fold change (PCH vs CTRL)"
    )

    plt.ylabel(
        "-log10(FDR)"
    )

    plt.title(
        f"Memento differential variability: {ct}"
    )

    plt.tight_layout()

    plt.savefig(
        f"Volcano_Memento_Variability_{ct.replace(' ', '_').replace('/', '_')}.png",
        dpi=300
    )

    plt.show()


# ==========================================================
# TOP DIFFERENTIAL MEAN GENES
# ==========================================================

for ct, res in results_all.items():

    top = (
        res[
            res["DE_significant"]
        ]
        .sort_values("DE_FDR")
        .head(20)
        .sort_values("log2FC")
    )

    if top.empty:
        continue

    plt.figure(figsize=(8, 6))

    plt.barh(
        top["gene"],
        top["log2FC"]
    )

    plt.xlabel(
        "log2 fold change (PCH vs CTRL)"
    )

    plt.title(
        f"Top Memento differential mean genes: {ct}"
    )

    plt.tight_layout()

    plt.show()


# ==========================================================
# TOP DIFFERENTIAL VARIABILITY GENES
# ==========================================================

for ct, res in results_all.items():

    top = (
        res[
            res["DV_significant"]
        ]
        .sort_values("DV_FDR")
        .head(20)
        .sort_values("log2_variance_FC")
    )

    if top.empty:
        continue

    plt.figure(figsize=(8, 6))

    plt.barh(
        top["gene"],
        top["log2_variance_FC"]
    )

    plt.xlabel(
        "log2 variance fold change (PCH vs CTRL)"
    )

    plt.title(
        f"Top Memento differential variability genes: {ct}"
    )

    plt.tight_layout()

    plt.show()










# ==========================================================
# DIFFERENTIAL GENE-GENE CORRELATION / COEXPRESSION
# ==========================================================

from itertools import combinations
from statsmodels.stats.multitest import multipletests

corr_num_boot = 1000

results_corr = {}

for ct in celltypes:

    print("\n===================================================")
    print(f"Running differential correlation for {ct}")

    subset = adata[
        adata.obs[celltype_col] == ct
    ].copy()

    print(f"Cells: {subset.n_obs}")

    # Skip tiny populations
    if subset.n_obs < 100:
        print("Skipping (<100 cells)")
        continue

    # Encode condition
    subset.obs["stim"] = (
        subset.obs[condition_col] == case
    ).astype(int)

    # Need both conditions
    if subset.obs["stim"].nunique() != 2:
        print("Skipping (only one condition present)")
        continue

    # --------------------------------------------
    # Candidate genes FOR THIS CELL TYPE ONLY
    # --------------------------------------------

    de = results_all[ct]

    # Top differential mean-expression genes
    top_de = (
        de[
            de["DE_significant"]
        ]
        .sort_values("DE_FDR")
        .head(50)
    )

    # Top differential variability genes
    top_dv = (
        de[
            de["DV_significant"]
        ]
        .sort_values("DV_FDR")
        .head(50)
    )

    # Union of both sets
    candidate_genes = sorted(
        set(top_de["gene"]).union(
            top_dv["gene"]
        )
    )

    print(
        f"Candidate genes: {len(candidate_genes)}"
    )

    print(
        "Top DE genes:",
        top_de["gene"].tolist()[:10]
    )

    print(
        "Top DV genes:",
        top_dv["gene"].tolist()[:10]
    )

    # Need at least two genes
    if len(candidate_genes) < 2:
        print("Skipping (not enough candidate genes)")
        continue

    # --------------------------------------------
    # Create gene pairs
    # --------------------------------------------

    gene_pairs = list(
        combinations(
            candidate_genes,
            2
        )
    )

    print(
        f"{ct}: {len(candidate_genes)} genes -> {len(gene_pairs)} gene pairs"
    )

    # --------------------------------------------
    # Run Memento
    # --------------------------------------------

    try:

        corr_res = memento.binary_test_2d(
            adata=subset,
            gene_pairs=gene_pairs,
            capture_rate=capture_rate,
            treatment_col="stim",
            num_boot=corr_num_boot,
            num_cpus=num_cpus
        )

    except Exception as e:

        print(f"Failed: {e}")
        continue

    print(
        f"Memento returned {len(corr_res)} tests"
    )

    if corr_res.empty:

        print("No valid gene pairs returned.")
        continue

    # --------------------------------------------
    # Multiple testing correction
    # --------------------------------------------

    corr_res["FDR"] = multipletests(
        corr_res["corr_pval"],
        method="fdr_bh"
    )[1]

    corr_res["celltype"] = ct

    # Standardized effect
    corr_res["z"] = np.where(
        corr_res["corr_se"] > 0,
        corr_res["corr_coef"] / corr_res["corr_se"],
        np.nan
    )

    # Significant differential correlation
    corr_res["significant"] = (
        (corr_res["FDR"] < fdr_thresh) &
        (corr_res["corr_se"] < 0.20)
    )

    results_corr[ct] = corr_res

    print(
        f"Significant differential correlations: "
        f"{corr_res['significant'].sum()}"
    )

print("\nDifferential correlation analysis finished.")

# ==========================================================
# Combine results
# ==========================================================

if len(results_corr) == 0:

    print("No correlation results were generated.")

else:

    combined_corr = pd.concat(
        results_corr.values(),
        ignore_index=True
    )

    sig_corr = combined_corr[
        combined_corr["significant"]
    ].sort_values(
        ["celltype", "FDR"]
    )

    combined_corr.to_csv(
        "Memento_all_differential_correlation_results.csv",
        index=False
    )

    sig_corr.to_csv(
        "Memento_significant_differential_correlations.csv",
        index=False
    )

    print()
    print(
        f"Total tests: {combined_corr.shape[0]}"
    )

    print(
        f"Significant correlations: {sig_corr.shape[0]}"
    )

    print(sig_corr.head(20))


#Visualize correlation:

#Volcano plots:
for ct, res in results_corr.items():

    plt.figure(figsize=(8,6))

    plt.scatter(
        res["corr_coef"],
        -np.log10(res["FDR"] + 1e-300),
        alpha=0.5
    )

    plt.axhline(
        -np.log10(fdr_thresh),
        color="black",
        linestyle="--"
    )

    top = (
        res[res["significant"]]
        .sort_values("FDR")
        .head(15)
    )

    for _, row in top.iterrows():

        plt.text(
            row["corr_coef"],
            -np.log10(row["FDR"] + 1e-300),
            f"{row['gene_1']}\n{row['gene_2']}",
            fontsize=7
        )

    plt.xlabel("Differential correlation coefficient")
    plt.ylabel("-log10(FDR)")
    plt.title(f"Differential correlation: {ct}")

    plt.tight_layout()
    plt.show()





#Heatmap:

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================================
# DIFFERENTIAL CORRELATION HEATMAPS
# ==========================================================

for ct in sorted(combined_corr["celltype"].unique()):

    print(f"Plotting {ct}")

    df = combined_corr[
        combined_corr["celltype"] == ct
    ].copy()

    # ----------------------------------------
    # Use significant pairs if available
    # ----------------------------------------

    plot_df = df[df["significant"]].copy()

    if plot_df.empty:

        print(
            "  No significant pairs; plotting top 50 by FDR."
        )

        plot_df = (
            df.sort_values("FDR")
              .head(50)
        )

    # ----------------------------------------
    # Get genes
    # ----------------------------------------

    genes = sorted(
        set(plot_df["gene_1"]).union(
            plot_df["gene_2"]
        )
    )

    n = len(genes)

    if n < 2:
        print("  Not enough genes.")
        continue

    gene_to_idx = {
        g: i for i, g in enumerate(genes)
    }

    # ----------------------------------------
    # Build matrix
    # ----------------------------------------

    mat = np.full((n, n), np.nan)

    np.fill_diagonal(mat, 0)

    for _, row in plot_df.iterrows():

        i = gene_to_idx[row["gene_1"]]
        j = gene_to_idx[row["gene_2"]]

        mat[i, j] = row["corr_coef"]
        mat[j, i] = row["corr_coef"]

    # ----------------------------------------
    # Plot
    # ----------------------------------------

    plt.figure(figsize=(10, 8))

    im = plt.imshow(
        mat,
        cmap="coolwarm",
    )

    plt.xticks(
        range(n),
        genes,
        rotation=90,
        fontsize=8
    )

    plt.yticks(
        range(n),
        genes,
        fontsize=8
    )

    plt.title(
        f"Differential gene-gene correlation\n{ct}"
    )

    plt.colorbar(
        im,
        label="Correlation coefficient change"
    )

    # Annotate values if matrix is small
    if n <= 20:

        for i in range(n):
            for j in range(n):

                if np.isnan(mat[i, j]):
                    continue

                plt.text(
                    j,
                    i,
                    f"{mat[i, j]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=6
                )

    plt.tight_layout()

    filename = (
        "Correlation_Heatmap_"
        + ct.replace("/", "_")
             .replace(" ", "_")
        + ".png"
    )

    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )

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
celltype_col = "cell_type"

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


    plt.show()

