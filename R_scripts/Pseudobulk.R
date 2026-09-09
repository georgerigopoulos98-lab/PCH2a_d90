#Packages
install.packages("Seurat")
install.packages("remotes")
install.packages("ggrepel")
devtools::install_github('immunogenomics/presto')
devtools::install_github('GfellerLab/EPIC')
devtools::install_github('quadbiolab/voxhunt')
remotes::install_github("mojaveazure/seurat-disk")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


#Load libraries
library(Seurat)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(Matrix)
library(ggrepel)
library(ggnewscale)


#Load h5ad file and convert to seurat_obj
sce <- readH5AD("PCH2a_d90_annotated.h5ad")
seurat_obj <- as.Seurat(
  sce,
  counts = "counts",
  data = "log1p_norm"
)

# Check the object
seurat_obj
dim(seurat_obj)
head(rownames(seurat_obj))
head(colnames(seurat_obj))

# View all metadata columns
colnames(seurat_obj@meta.data)

# Look at sample/condition information
head(seurat_obj@meta.data)

# Check amount of cells:
table(seurat_obj$condition)
table(seurat_obj$cell_type)


# Get raw counts from the RNA assay
counts <- LayerData(
  seurat_obj,
  assay = "originalexp",
  layer = "counts"
)

# Check dimensions
dim(counts)
counts[1:5, 1:5]


# Analysis parameters:
# Set biological coefficient of variation
BCV <- 0.30

# Convert BCV to dispersion
dispersion_value <- BCV^2

cat("BCV:", BCV, "\n")
cat("Dispersion:", dispersion_value, "\n")


# Extract raw counts from Seurat object:
counts <- LayerData(
  object = seurat_obj,
  assay = "originalexp",
  layer = "counts"
)

# Check dimensions
cat(
  "Number of genes:",
  nrow(counts),
  "\n"
)

cat(
  "Number of cells:",
  ncol(counts),
  "\n"
)

# Extract metadata
metadata <- seurat_obj@meta.data %>%
  dplyr::select(
    sample_id,
    condition,
    cell_type
  )

# Make sure counts and metadata are in the same order
if (!all(colnames(counts) == rownames(metadata))) {
  
  stop(
    "Cell names/order in counts do not match metadata."
  )
  
}

# Set condition factor

metadata$condition <- factor(
  metadata$condition,
  levels = c(
    "CTRL",
    "PCH"
  )
)

# Check sample/condition structure

cat("\nSample by condition:\n")

print(
  table(
    metadata$sample_id,
    metadata$condition
  )
)

# Check cell numbers per cell type and condition

cat("\nCell counts by cell type and condition:\n")

print(
  table(
    metadata$cell_type,
    metadata$condition
  )
)


# Get list of cell types

cell_types <- sort(
  unique(
    metadata$cell_type
  )
)

cat("\nCell types to analyze:\n")

print(cell_types)

# Create output directories

dir.create(
  "edgeR_DGE_results",
  showWarnings = FALSE
)

dir.create(
  "edgeR_DGE_results/Volcano_plots",
  showWarnings = FALSE
)

# Initialize results list

dge_results <- list()

# Initialize list to store genes retained after filtering for each cell type

tested_gene_universe <- list()

# Run edgeR separately for each cell type

for (ct in cell_types) {
  
  cat("\n")
  cat("============================================================\n")
  cat("Analyzing cell type:", ct, "\n")
  cat("============================================================\n")
  
  # Select cells belonging to this cell type

  cells_ct <- rownames(metadata)[
    metadata$cell_type == ct
  ]

  # Extract counts for this cell type

  counts_ct <- counts[
    ,
    cells_ct,
    drop = FALSE
  ]

  # Extract metadata for this cell type

  meta_ct <- metadata[
    cells_ct,
    ,
    drop = FALSE
  ]

  # Check condition counts

  cat("\nNumber of cells per condition:\n")
  
  print(
    table(
      meta_ct$condition
    )
  )

  # Skip cell type if one condition is missing

  if (length(unique(meta_ct$condition)) < 2) {
    
    cat(
      "Skipping",
      ct,
      "- fewer than two conditions are present.\n"
    )
    
    next
  }

  # Create edgeR DGEList

  y <- edgeR::DGEList(
    counts = counts_ct,
    group = meta_ct$condition
  )

  # Filter lowly expressed genes.
  # Keep genes expressed at CPM > 1 in at least 10 cells.
  
  cpm_values <- edgeR::cpm(
    y
  )
  
  keep <- rowSums(
    cpm_values > 1
  ) >= 10
  
  # Save genes retained after filtering.
  # These genes represent the tested gene universe for downstream GO enrichment analysis.

  tested_gene_universe[[ct]] <- rownames(y)[keep]
  
  cat(
    "Genes retained after filtering:",
    length(tested_gene_universe[[ct]]),
    "\n"
  )

  # Apply filter

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  cat(
    "Genes retained after filtering:",
    nrow(y),
    "\n"
  )

  # TMM normalization

  y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
  )

  # Create design matrix.
  # CTRL is the reference level.

  design <- model.matrix(
    ~ condition,
    data = meta_ct
  )
  
  cat("\nDesign matrix:\n")
  
  print(
    colnames(design)
  )

  # Fit quasi-likelihood negative binomial GLM
  # Fixed BCV = 0.30
  # Fixed dispersion = 0.09
  # We do not estimate dispersion from biological replicates because there is only one sample per condition.

  fit <- edgeR::glmQLFit(
    y,
    design = design,
    dispersion = dispersion_value,
    robust = TRUE
  )

  # Perform quasi-likelihood F-test
  
  qlf <- edgeR::glmQLFTest(
    fit,
    coef = "conditionPCH"
  )

  # Extract results
  # topTags() returns:
  # logFC
  # logCPM
  # F = quasi-likelihood F-statistic
  # PValue
  # FDR

  results <- edgeR::topTags(
    qlf,
    n = Inf,
    sort.by = "none"
  )$table

  # Add gene names

  results$gene <- rownames(
    results
  )
  
  # Move gene column to first position
  
  results <- results[
    ,
    c(
      "gene",
      setdiff(
        colnames(results),
        "gene"
      )
    )
  ]

  # Add cell type

  results$cell_type <- ct

  # Add analysis parameters

  results$BCV <- BCV
  
  results$dispersion <- dispersion_value

  # Define exploratory significance categories.
  # FDR < 0.05 & absolute logFC > 1
  # The test itself is based on the QLF framework.

  results$significance <- "Not_significant"
  
  results$significance[
    results$FDR < 0.05 &
      results$logFC > 1
  ] <- "Up_in_PCH"
  
  results$significance[
    results$FDR < 0.05 &
      results$logFC < -1
  ] <- "Down_in_PCH"

  # Rank genes by QLF F-statistic

  results <- results[
    order(
      results$F,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  # Store results

  dge_results[[ct]] <- results

  # Create safe filename

  safe_ct <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    ct
  )

  # Save results for this cell type

  output_file <- file.path(
    "edgeR_DGE_results",
    paste0(
      "edgeR_",
      safe_ct,
      "_PCH_vs_CTRL.csv"
    )
  )
  
  write.csv(
    results,
    file = output_file,
    row.names = FALSE
  )

  # Print top 20 genes ranked by QLF F-statistic

  cat(
    "\nTop 20 genes ranked by QLF F-statistic:\n"
  )
  
  print(
    head(
      results[
        ,
        c(
          "gene",
          "logFC",
          "logCPM",
          "F",
          "PValue",
          "FDR"
        )
      ],
      20
    )
  )
  
}


# Combine results from all cell types

all_results <- bind_rows(
  dge_results
)

# Save combined results

write.csv(
  all_results,
  file = "edgeR_DGE_results/edgeR_all_cell_types_PCH_vs_CTRL.csv",
  row.names = FALSE
)

# Create summary table

dge_summary <- all_results %>%
  group_by(cell_type) %>%
  summarise(
    n_genes = n(),
    
    n_FDR_05 = sum(
      FDR < 0.05,
      na.rm = TRUE
    ),
    
    n_up = sum(
      FDR < 0.05 &
        logFC > 1,
      na.rm = TRUE
    ),
    
    n_down = sum(
      FDR < 0.05 &
        logFC < -1,
      na.rm = TRUE
    ),
    
    max_F = max(
      F,
      na.rm = TRUE
    ),
    
    BCV = dplyr::first(BCV),
    
    dispersion = dplyr::first(dispersion),
    
    .groups = "drop"
  )

# Print summary

cat("\nSummary of DGE results:\n")

print(
  dge_summary
)

# Save summary

write.csv(
  dge_summary,
  file = "edgeR_DGE_results/edgeR_summary.csv",
  row.names = FALSE
)


# Generate volcano plots

# Genes highlighted as significant:
# FDR < 0.05
# absolute logFC > 1

for (ct in names(dge_results)) {

  # Extract results for this cell type
  
  results <- dge_results[[ct]]

  # Create safe filename

  safe_ct <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    ct
  )
  
  #Identify top 10 DE genes
  
  top10_genes <- results %>%
    dplyr::arrange(
      dplyr::desc(F)
    ) %>%
    dplyr::slice_head(
      n = 10
  )

  # Create QLF volcano plot

  p <- ggplot(
    results,
    aes(
      x = logFC,
      y = -log10(PValue)
    )
  ) +
    
  # Plot all genes
  geom_point(
    aes(
      color = significance
    ),
    alpha = 0.6,
    size = 1.5
  ) +
    
    # Add vertical lines for logFC thresholds
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    
    #Label top 10 DE genes
    ggrepel::geom_text_repel( 
      data = top10_genes,
      aes(
        label = gene
      ), 
      size = 3, 
      box.padding = 0.5, 
      point.padding = 0.3, 
      max.overlaps = Inf, 
      show.legend = FALSE 
    ) +
    
    # Apply a clean theme
    theme_classic() +
    
    # Set colors 
    scale_color_manual( 
      values = c( 
        "Not_significant" = "grey70", 
        "Up_in_PCH" = "red", 
        "Down_in_PCH" = "blue" 
      ) 
    ) +
    
    # Add labels
    labs(
      title = paste(
        "edgeR QLF:",
        ct,
        "PCH vs CTRL"
      ),
      
      subtitle = paste(
        "TOP 10 genes ranked by QLF F-statistic |",
        "FDR < 0.05 and |logFC| > 1"
      ),
      
      x = "log2 Fold Change (PCH vs CTRL)",
      
      y = "-log10 Quasi-likelihood P-value",
      
      color = "Significance"
    )
  
  # Save volcano plot
  
  ggsave(
    filename = file.path(
      "edgeR_DGE_results",
      "Volcano_plots",
      paste0(
        "Volcano_",
        safe_ct,
        ".pdf"
      )
    ),
    
    plot = p,
    
    width = 7,
    
    height = 6
    
  )
  
  # Print plot to R session

  print(p)
  
}


# Save Seurat object and results:

saveRDS(
  dge_results,
  file = "edgeR_DGE_results/dge_results.rds"
)

saveRDS(
  all_results,
  file = "edgeR_DGE_results/all_results.rds"
)

saveRDS(
  dge_summary,
  file = "edgeR_DGE_results/dge_summary.rds"
)

saveRDS(
  seurat_obj,
  file = "edgeR_DGE_results/seurat_obj_after_edgeR.rds"
)

saveRDS(
  tested_gene_universe,
  file = "edgeR_DGE_results/tested_gene_universe.rds"
)




#Check PCH related genes specifically:

# Specify genes of interest

genes_of_interest <- c(
  "TSEN54",
  "TSEN2",
  "TSEN34",
  "TSEN15",
  "SARS2",
  "RARS2",
  "SEPSECS",
  "EXOSC2",
  "EXOSC3",
  "EXOSC8",
  "EXOSC9",
  "TOE1",
  "EXOSC1",
  "PPIL1",
  "CDC40",
  "PCLO",
  "SLC25A46",
  "COASY",
  "TBC1D20",
  "TBC1D23",
  "VRK1",
  "AMPD2",
  "CHMP1A",
  "MINPP1",
  "VPS51",
  "VPS53",
  "ASAH1",
  "TWNK",
  "SPTBN2",
  "PMPCA"
)

# Check that requested genes are present

all_genes <- unique(
  unlist(
    lapply(
      dge_results,
      function(x) x$gene
    )
  )
)

genes_found <- genes_of_interest[
  genes_of_interest %in% all_genes
]

genes_missing <- genes_of_interest[
  !genes_of_interest %in% all_genes
]

cat(
  "\nGenes found:",
  length(genes_found),
  "\n"
)

if (length(genes_missing) > 0) {
  
  cat(
    "\nGenes not found in edgeR results:\n"
  )
  
  print(
    genes_missing
  )
  
}

# Combine edgeR results from all cell types

all_results <- dplyr::bind_rows(
  dge_results
)

# Extract genes of interest

gene_results <- all_results %>%
  dplyr::filter(
    gene %in% genes_found
  ) %>%
  dplyr::select(
    gene,
    cell_type,
    logFC,
    logCPM,
    F,
    PValue,
    FDR,
    significance
  )

# Sort results

gene_results <- gene_results %>%
  dplyr::arrange(
    gene,
    cell_type
  )

# Print results

cat(
  "\n============================================================\n"
)

cat(
  "Gene-level edgeR results:\n"
)

cat(
  "============================================================\n"
)

print(
  gene_results
)

# Save detailed results

write.csv(
  gene_results,
  file = "edgeR_DGE_results/genes_of_interest_by_cell_type.csv",
  row.names = FALSE
)

# Create summary across cell types

gene_summary <- gene_results %>%
  dplyr::group_by(
    gene
  ) %>%
  dplyr::summarise(
    
    n_cell_types = dplyr::n(),
    
    n_up = sum(
      significance == "Up_in_PCH",
      na.rm = TRUE
    ),
    
    n_down = sum(
      significance == "Down_in_PCH",
      na.rm = TRUE
    ),
    
    n_significant = sum(
      significance != "Not_significant",
      na.rm = TRUE
    ),
    
    mean_logFC = mean(
      logFC,
      na.rm = TRUE
    ),
    
    median_logFC = median(
      logFC,
      na.rm = TRUE
    ),
    
    min_FDR = min(
      FDR,
      na.rm = TRUE
    ),
    
    .groups = "drop"
    
  )

# Print summary

cat(
  "\n============================================================\n"
)

cat(
  "Summary across cell types:\n"
)

cat(
  "============================================================\n"
)

print(
  gene_summary
)

# Save summary

write.csv(
  gene_summary,
  file = "edgeR_DGE_results/genes_of_interest_summary.csv",
  row.names = FALSE
)

# Create logFC matrix
# Rows = genes
# Columns = cell types
# Values = logFC (PCH vs CTRL)

logFC_matrix <- gene_results %>%
  dplyr::select(
    gene,
    cell_type,
    logFC
  ) %>%
  tidyr::pivot_wider(
    names_from = cell_type,
    values_from = logFC
  )

# Convert to matrix

logFC_matrix <- as.data.frame(
  logFC_matrix
)

rownames(
  logFC_matrix
) <- logFC_matrix$gene

logFC_matrix$gene <- NULL

logFC_matrix <- as.matrix(
  logFC_matrix
)

# Save logFC matrix

write.csv(
  logFC_matrix,
  file = "edgeR_DGE_results/genes_of_interest_logFC_matrix.csv"
)

# Create logFC heatmap
# Red = Upregulated in PCH
# Blue = Upregulated in CTRL
# Colour intensity reflects magnitude of logFC
# Convert logFC matrix to long format

logFC_long <- as.data.frame(
  logFC_matrix
)

logFC_long$gene <- rownames(
  logFC_long
)

logFC_long <- logFC_long %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "cell_type",
    values_to = "logFC"
  )

# Create direction variable

logFC_long$direction <- dplyr::case_when(
  
  is.na(logFC_long$logFC) ~ "No_data",
  
  logFC_long$logFC > 0 ~ "Up_in_PCH",
  
  logFC_long$logFC < 0 ~ "Up_in_CTRL",
  
  TRUE ~ "No_change"
  
)

# Create heatmap

p_logFC <- ggplot(
  logFC_long,
  aes(
    x = cell_type,
    y = gene,
    fill = logFC
  )
) +
  
  geom_tile(
    color = "white"
  ) +
  
  # Use a diverging colour scale
  # Blue = Up in CTRL
  # White = No change
  # Red = Up in PCH
  # midpoint = 0 ensures that the direction of the fold change is visually clear.
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 8
    )
    
  ) +
  
  labs(
    title = "PCH vs CTRL: Genes of Interest",
    
    subtitle = paste(
      "Red = Up in PCH |",
      "Blue = Up in CTRL"
    ),
    
    x = "Cell Type",
    
    y = "Gene",
    
    fill = "log2 Fold Change"
    
  )

# Save logFC heatmap

ggsave(
  filename = "edgeR_DGE_results/genes_of_interest_logFC_heatmap.pdf",
  plot = p_logFC,
  width = 10,
  height = 6
)

# Print heatmap

print(
  p_logFC
)



# Create FDR heatmap.
# Red = Upregulated in PCH
# Blue = Upregulated in CTRL
# Gray = Not significant
# FDR < 0.05 = Significant
# Colour intensity reflects statistical significance.
# Darker colour = smaller FDR

# Create plotting data directly from gene_results.
# This retains:
# logFC = direction of change
# FDR   = statistical significance

FDR_long <- gene_results %>%
  dplyr::select(
    gene,
    cell_type,
    logFC,
    FDR
  )

# Calculate -log10(FDR)
# Smaller FDR values produce larger values.
# This determines colour intensity.

FDR_long$neglog10FDR <- -log10(
  FDR_long$FDR
)

# Create significance categories

FDR_long$significant <- dplyr::case_when(
  
  is.na(FDR_long$FDR) ~ "No_data",
  
  FDR_long$FDR < 0.05 ~ "Significant",
  
  TRUE ~ "Not_significant"
  
)

# Create FDR heatmap

p_FDR <- ggplot() +

# Background:
# All genes start as gray.
# This represents non-significant genes.

geom_tile(
  data = FDR_long,
  aes(
    x = cell_type,
    y = gene
  ),
  fill = "grey80",
  color = "white"
) +
  
# Significant genes upregulated in PCH
# Red colour gradient.

geom_tile(
  data = FDR_long %>%
    dplyr::filter(
      FDR < 0.05 &
        logFC > 0
    ),
  aes(
    x = cell_type,
    y = gene,
    fill = neglog10FDR
  ),
  color = "white"
) +
  
  scale_fill_gradient(
    low = "lightpink",
    high = "red",
    name = "-log10(FDR)\nUp in PCH"
  ) +
  
# Add a new fill scale
  
newscale::new_scale_fill() +

# Significant genes upregulated in CTRL
# Blue colour gradient.
  
geom_tile(
  data = FDR_long %>%
    dplyr::filter(
      FDR < 0.05 &
        logFC < 0
    ),
  aes(
    x = cell_type,
    y = gene,
    fill = neglog10FDR
  ),
  color = "white"
) +
  
  scale_fill_gradient(
    low = "lightblue",
    high = "blue",
    name = "-log10(FDR)\nUp in CTRL"
  ) +
  
# Theme

theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 8
    )
    
  ) +
  
# Labels

labs(
  title = "Statistical Significance of Genes of Interest",
  
  subtitle = paste(
    "Red = Up in PCH |",
    "Blue = Up in CTRL |",
    "Gray = FDR ≥ 0.05"
  ),
  
  x = "Cell Type",
  
  y = "Gene"
  
)

# Save FDR heatmap

ggsave(
  filename = "edgeR_DGE_results/genes_of_interest_FDR_heatmap.pdf",
  plot = p_FDR,
  width = 10,
  height = 6
)

# Print FDR heatmap

print(
  p_FDR
)

