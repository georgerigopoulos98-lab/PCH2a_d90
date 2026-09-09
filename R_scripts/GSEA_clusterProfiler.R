if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")


# Reload objects from previous edgeR analysis

seurat_obj <- readRDS(
  "edgeR_DGE_results/seurat_obj_after_edgeR.rds"
)

dge_results <- readRDS(
  "edgeR_DGE_results/dge_results.rds"
)

all_results <- readRDS(
  "edgeR_DGE_results/all_results.rds"
)

dge_summary <- readRDS(
  "edgeR_DGE_results/dge_summary.rds"
)

# Check that objects loaded correctly

cat(
  "Number of cell types with DGE results:",
  length(dge_results),
  "\n"
)

cat(
  "Cell types:\n"
)

print(
  names(dge_results)
)


# Load packages

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)


# Specify folder containing GMT files

gmt_folder <- "data/gsea_references"


# Find all GMT files in folder

gmt_files <- list.files(
  path = gmt_folder,
  pattern = "\\.gmt$",
  full.names = TRUE,
  ignore.case = TRUE
)


# Check GMT files

if (length(gmt_files) == 0) {
  
  stop(
    "No .gmt files were found in the specified folder."
  )
  
}


cat(
  "Number of GMT files found:",
  length(gmt_files),
  "\n"
)


cat(
  "\nGMT files:\n"
)


print(
  basename(gmt_files)
)


# Create output directories

dir.create(
  "edgeR_DGE_results/GSEA",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "edgeR_DGE_results/GSEA/Custom_GMT",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "edgeR_DGE_results/GSEA/Plots",
  showWarnings = FALSE,
  recursive = TRUE
)


# Results will be stored as:

gsea_results_all <- list()


# Loop through each GMT file

for (gmt_file in gmt_files) {
  
  # Get GMT filename

  gmt_name <- tools::file_path_sans_ext(
    basename(gmt_file)
  )
  
  # Create safe GMT name
  
  safe_gmt_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    gmt_name
  )
  
  
  cat("\n")
  cat("============================================================\n")
  cat(
    "Processing GMT file:",
    gmt_name,
    "\n"
  )
  cat("============================================================\n")
  
  # Read GMT file

  gene_sets <- clusterProfiler::read.gmt(
    gmt_file
  )
  
  # Check GMT file

  cat(
    "Number of gene sets:",
    length(unique(gene_sets$term)),
    "\n"
  )
  
  
  cat(
    "Number of gene-set/gene associations:",
    nrow(gene_sets),
    "\n"
  )
  
  # Convert GMT gene symbols to Entrez IDs (the GMT files contain HGNC gene symbols).

  cat(
    "\nConverting GMT gene symbols to Entrez IDs...\n"
  )
  
  gene_sets_entrez <- clusterProfiler::bitr(
    unique(gene_sets$gene),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Merge Entrez IDs with GMT gene sets

  gene_sets_entrez <- dplyr::inner_join(
    gene_sets,
    gene_sets_entrez,
    by = c(
      "gene" = "SYMBOL"
    )
  )
  
  # Create TERM2GENE object.
  # TERM2GENE must contain:
  # Column 1 = pathway/gene-set name
  # Column 2 = gene ID
  
  TERM2GENE <- gene_sets_entrez %>%
    dplyr::select(
      term,
      ENTREZID
    ) %>%
    dplyr::distinct()
  
  # Check GMT conversion

  cat(
    "Number of gene sets after ID conversion:",
    length(unique(TERM2GENE$term)),
    "\n"
  )
  
  
  cat(
    "Number of gene-set/gene associations after conversion:",
    nrow(TERM2GENE),
    "\n"
  )
  
  # Create GMT-specific output directories
  
  gmt_output_dir <- file.path(
    "edgeR_DGE_results",
    "GSEA",
    "Custom_GMT",
    safe_gmt_name
  )
  
  
  gmt_plot_dir <- file.path(
    "edgeR_DGE_results",
    "GSEA",
    "Plots",
    safe_gmt_name
  )
  
  
  dir.create(
    gmt_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  dir.create(
    gmt_plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  

  # Initialize results for this GMT

  gsea_results_all[[gmt_name]] <- list()
  
  # Run GSEA separately for each cell type

  for (ct in names(dge_results)) {
    
    cat("\n")
    cat(
      "Running GSEA:",
      gmt_name,
      "| Cell type:",
      ct,
      "\n"
    )

    # Extract edgeR results

    results <- dge_results[[ct]]
    
    # Create signed QLF ranking statistic.
    # Positive: Upregulated in PCH
    # Negative: Upregulated in CTRL

    results$rank_statistic <- sign(
      results$logFC
    ) * results$F
    
    # Remove missing ranking statistics

    results <- results[
      !is.na(results$rank_statistic),
      ,
      drop = FALSE
    ]

    # Convert edgeR gene symbols to Entrez IDs

    gene_conversion <- clusterProfiler::bitr(
      unique(results$gene),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    # Merge Entrez IDs with edgeR results

    results_converted <- dplyr::inner_join(
      results,
      gene_conversion,
      by = c(
        "gene" = "SYMBOL"
      )
    )

    # Remove duplicated Entrez IDs.
    # Keep gene with largest absolute ranking statistic.

    results_converted <- results_converted %>%
      dplyr::arrange(
        dplyr::desc(
          abs(rank_statistic)
        )
      ) %>%
      dplyr::distinct(
        ENTREZID,
        .keep_all = TRUE
      )
    
    # Create ranked gene list.
    # Names = Entrez IDs
    # Values = signed QLF F-statistic

    gene_list <- results_converted$rank_statistic
    
    
    names(gene_list) <- results_converted$ENTREZID
    
    # Sort ranking

    gene_list <- sort(
      gene_list,
      decreasing = TRUE
    )
    
    # Print conversion statistics

    cat(
      "Genes before ID conversion:",
      nrow(results),
      "\n"
    )
    
    
    cat(
      "Genes successfully converted:",
      length(gene_list),
      "\n"
    )

    # Calculate conversion rate

    conversion_rate <- (
      length(gene_list) /
        nrow(results)
    ) * 100
    
    
    cat(
      "Conversion rate:",
      round(
        conversion_rate,
        1
      ),
      "%\n"
    )
    
    # Skip cell type if too few genes

    if (length(gene_list) < 100) {
      
      cat(
        "Skipping",
        ct,
        "- fewer than 100 genes available.\n"
      )
      
      next
      
    }
    
    # Run GSEA using custom GMT
    # pvalueCutoff = 1 retains all results
    # Significant pathways can be filtered later using adjusted p-value.

    gsea_custom <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = TERM2GENE,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    # Store results

    gsea_results_all[[gmt_name]][[ct]] <- gsea_custom

    # Convert results to data frame

    gsea_results <- as.data.frame(
      gsea_custom
    )

    # Create safe cell type name

    safe_ct <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      ct
    )
    
    # Save all GSEA results

    write.csv(
      gsea_results,
      file = file.path(
        gmt_output_dir,
        paste0(
          "GSEA_",
          safe_ct,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    # Identify significant pathways, FDR < 0.05

    significant_results <- gsea_results[
      !is.na(gsea_results$p.adjust) &
        gsea_results$p.adjust < 0.05,
      ,
      drop = FALSE
    ]
    
    # Save significant pathways

    write.csv(
      significant_results,
      file = file.path(
        gmt_output_dir,
        paste0(
          "GSEA_",
          safe_ct,
          "_FDR05.csv"
        )
      ),
      row.names = FALSE
    )

    # Plot GSEA results
    
    gsea_results_FDR05 <- gsea_results[
      !is.na(gsea_results$p.adjust) &
        gsea_results$p.adjust < 0.05,
      ,
      drop = FALSE
    ]
    
    if (nrow(gsea_results_FDR05) > 0) {
      
      gsea_custom_FDR05 <- gsea_custom
      gsea_custom_FDR05@result <- gsea_results_FDR05
      
      p_gsea <- clusterProfiler::dotplot(
        gsea_custom_FDR05,
        x = "NES",
        showCategory = 20
      ) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed"
        ) +
        ggplot2::ggtitle(
          paste(
            gmt_name,
            ":",
            ct,
            "PCH vs CTRL"
          )
        ) +
        ggplot2::xlab(
          "Normalized Enrichment Score (NES)"
        ) +
        ggplot2::scale_y_discrete(
          labels = function(x) {
            stringr::str_wrap(x, width = 45)
          }
        ) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_text(size = 10)
        )
      
      # Display plot in R
      print(p_gsea)
      
    } else {
      
      cat(
        "No pathways with FDR < 0.05 for:",
        gmt_name,
        "|",
        ct,
        "\n"
      )
    }
    
    # Print top pathways

    cat(
      "\nTop pathways:\n"
    )
    
    
    if (nrow(significant_results) > 0) {
      
      top_results <- significant_results[
        order(
          significant_results$p.adjust
        ),
        c(
          "ID",
          "Description",
          "NES",
          "p.adjust"
        ),
        drop = FALSE
      ]
      
      
      print(
        head(
          top_results,
          10
        )
      )
      
    } else {
      
      cat(
        "No pathways detected.\n"
      )
      
    }
    
    # Print number of significant pathways

    cat(
      "\nNumber of pathways with FDR < 0.05:",
      nrow(significant_results),
      "\n"
    )
    
  }
  
}

