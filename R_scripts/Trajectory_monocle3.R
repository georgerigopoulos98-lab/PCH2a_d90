if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

install.packages("devtools")
remotes::install_github("bnprks/BPCells/r")

devtools::install_github('cole-trapnell-lab/monocle3')

remotes::install_github("satijalab/seurat-wrappers")


# Monocle3 Pseudotime / Trajectory Analysis

# Load packages

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(SeuratWrappers)
library(tidyr)
library(pheatmap)


# Create output directories

dir.create(
  "Monocle3_trajectory",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Monocle3_trajectory/Plots",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Monocle3_trajectory/Results",
  recursive = TRUE,
  showWarnings = FALSE
)

# Check metadata

cat("Available cell types:\n")
print(unique(seurat_obj$cell_type))

cat("\nCells by cell type:\n")
print(table(seurat_obj$cell_type))

cat("\nCells by condition:\n")
print(table(seurat_obj$condition))


# Convert FULL Seurat object to Monocle 3
# IMPORTANT:
# Do NOT subset by cell type.

cds <- SeuratWrappers::as.cell_data_set(
  seurat_obj
)

# Add metadata explicitly

colData(cds)$cell_type <- seurat_obj$cell_type

colData(cds)$condition <- seurat_obj$condition

colData(cds)$sample_id <- seurat_obj$sample_id

# Transfer EXISTING Seurat UMAP

seurat_umap <- Seurat::Embeddings(
  seurat_obj,
  reduction = "X_umap_scVI"
)

# Make sure cell order matches
seurat_umap <- seurat_umap[
  colnames(cds),
  ,
  drop = FALSE
]

# Transfer UMAP to Monocle 3
reducedDims(cds)$UMAP <- seurat_umap

# Check that UMAP was transferred correctly

cat(
  "\nNumber of cells in Seurat UMAP:",
  nrow(seurat_umap),
  "\n"
)

cat(
  "Number of cells in Monocle 3:",
  ncol(cds),
  "\n"
)

stopifnot(
  identical(
    rownames(seurat_umap),
    colnames(cds)
  )
)

cat("\nExisting Seurat UMAP successfully transferred.\n")

# Cluster cells using existing UMAP

cds <- monocle3::cluster_cells(
  cds,
  reduction_method = "UMAP"
)

# Learn ONE trajectory graph across ALL cell types

cds <- monocle3::learn_graph(
  cds,
  use_partition = FALSE
)

# Plot trajectory colored by cell type

p_celltype <- monocle3::plot_cells(
  cds,
  color_cells_by = "cell_type",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE
) +
  ggplot2::ggtitle(
    "Monocle3 trajectory across all cell types"
  )

ggsave(
  filename = file.path(
    "Monocle3_trajectory",
    "Plots",
    "Trajectory_All_CellTypes.pdf"
  ),
  plot = p_celltype,
  width = 10,
  height = 8
)

# Plot trajectory colored by condition

p_condition <- monocle3::plot_cells(
  cds,
  color_cells_by = "condition",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
) +
  ggplot2::ggtitle(
    "Monocle3 trajectory by condition"
  )

ggsave(
  filename = file.path(
    "Monocle3_trajectory",
    "Plots",
    "Trajectory_All_CellTypes_Condition.pdf"
  ),
  plot = p_condition,
  width = 10,
  height = 8
)

# Choose root cells for pseudotime

root_cells <- colnames(cds)[
  colData(cds)$condition == "CTRL" &
  colData(cds)$cell_type == "Progenitor/ Astrocytes"
]

# Check root cell number

cat(
  "\nNumber of CTRL cells used as root candidates:",
  length(root_cells),
  "\n"
)

if (
  length(root_cells) == 0
) {
  
  stop(
    "No CTRL cells were found for root selection."
  )
  
}

# Order cells in pseudotime

cds <- monocle3::order_cells(
  cds,
  root_cells = root_cells
)

# Extract pseudotime values

pseudotime_values <- monocle3::pseudotime(
  cds
)

# Add pseudotime to metadata

colData(cds)$pseudotime <- pseudotime_values

# Save pseudotime values

pseudotime_results <- data.frame(
  cell = colnames(cds),
  pseudotime = pseudotime_values,
  condition = colData(cds)$condition,
  sample_id = colData(cds)$sample_id
)

write.csv(
  pseudotime_results,
  file = file.path(
    "Monocle3_trajectory",
    "Results",
    paste0(
      "Pseudotime_All_CellTypes.csv"
    )
  ),
  row.names = FALSE
)

# Plot cells colored by pseudotime

p_pseudotime <- monocle3::plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
) +
  ggplot2::ggtitle(
    paste(
      "Pseudotime trajectory"
    )
  )

ggsave(
  filename = file.path(
    "Monocle3_trajectory",
    "Plots",
    paste0(
      "Pseudotime.pdf"
    )
  ),
  plot = p_pseudotime,
  width = 8,
  height = 7
)

# Plot pseudotime by condition

pseudotime_plot_data <- data.frame(
  pseudotime = pseudotime_values,
  condition = colData(cds)$condition
)

p_condition_pseudotime <- ggplot(
  pseudotime_plot_data,
  aes(
    x = condition,
    y = pseudotime
  )
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.2,
    alpha = 0.3
  ) +
  theme_classic() +
  labs(
    title = paste(
      "Pseudotime by condition"
    ),
    x = "Condition",
    y = "Pseudotime"
  )

ggsave(
  filename = file.path(
    "Monocle3_trajectory",
    "Plots",
    paste0(
      "Pseudotime_by_Condition.pdf"
    )
  ),
  plot = p_condition_pseudotime,
  width = 7,
  height = 6
)

# Identify genes associated with pseudotime

gene_fits <- monocle3::graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores = 1
)

# Convert gene results to data frame

gene_fits <- as.data.frame(
  gene_fits
)

# Sort genes by adjusted p-value

gene_fits <- gene_fits[
  order(
    gene_fits$q_value
  ),
  ,
  drop = FALSE
]

# Save all pseudotime-associated genes

write.csv(
  gene_fits,
  file = file.path(
    "Monocle3_trajectory",
    "Results",
    paste0(
      "Pseudotime_Associated_Genes.csv"
    )
  ),
  row.names = TRUE
)

# Identify significant pseudotime-associated genes
# q-value < 0.05

significant_genes <- gene_fits[
  !is.na(
    gene_fits$q_value
  ) &
    gene_fits$q_value < 0.05,
  ,
  drop = FALSE
]

# Save significant pseudotime-associated genes

write.csv(
  significant_genes,
  file = file.path(
    "Monocle3_trajectory",
    "Results",
    paste0(
      "Pseudotime_Associated_Genes_FDR05.csv"
    )
  ),
  row.names = TRUE
)

# Print summary

cat(
  "\nNumber of genes tested for pseudotime association:",
  nrow(gene_fits),
  "\n"
)

cat(
  "Number of significant pseudotime-associated genes:",
  nrow(significant_genes),
  "\n"
)

# Plot top pseudotime-associated genes

rowData(cds)$gene_short_name <- rownames(cds)

if (
  nrow(significant_genes) > 0
) {
  
  top_genes <- rownames(
    head(
      significant_genes,
      10
    )
  )
  
  p_genes <- monocle3::plot_genes_in_pseudotime(
    cds[
      top_genes,
    ],
    min_expr = 0.5
  )
  
  ggsave(
    filename = file.path(
      "Monocle3_trajectory",
      "Plots",
      paste0(
        "Top_Pseudotime_Genes.pdf"
      )
    ),
    plot = p_genes,
    width = 10,
    height = 8
  )
  
}

# Save Monocle3 object

monocle3::save_monocle_objects(
  cds,
  file.path(
    "Monocle3_trajectory",
    "Results",
    paste0(
      "Monocle3"
    )
  )
)

