
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


# Set significance thresholds

FDR_cutoff <- 0.05

logFC_cutoff <- 1


# Load cell-type-specific tested gene universes.
# These were generated during the edgeR analysis after CPM filtering.
# Each cell type has its own tested gene universe.

tested_gene_universe <- readRDS(
  "edgeR_DGE_results/tested_gene_universe.rds"
)


# Check loaded gene universes

cat(
  "Cell types with saved gene universes:",
  length(tested_gene_universe),
  "\n"
)

cat(
  "\nCell types:\n"
)

print(
  names(tested_gene_universe)
)


# Create output directories

dir.create(
  "edgeR_DGE_results/GO_Enrichment",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "edgeR_DGE_results/GO_Enrichment/Upregulated",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "edgeR_DGE_results/GO_Enrichment/Downregulated",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "edgeR_DGE_results/GO_Enrichment/Plots",
  showWarnings = FALSE,
  recursive = TRUE
)


# Results will be stored as:

go_results_all <- list()

go_results_all$Upregulated <- list()

go_results_all$Downregulated <- list()


# Run GO enrichment separately for each cell type

for (ct in names(dge_results)) {
  
  cat("\n")
  cat("============================================================\n")
  cat(
    "Running GO BP enrichment for cell type:",
    ct,
    "\n"
  )
  cat("============================================================\n")
  
  # Check that a tested gene universe exists for this cell type

  if (
    !ct %in% names(tested_gene_universe)
  ) {
    
    warning(
      paste(
        "No tested gene universe found for cell type:",
        ct,
        "- skipping."
      )
    )
    
    next
    
  }
  
  # Extract edgeR results
  
  results <- dge_results[[ct]]
  
  # Convert tested gene universe to Entrez IDs
  
  tested_genes <- tested_gene_universe[[ct]]
  
  cat(
    "\nGenes retained after edgeR filtering:",
    length(
      tested_genes
    ),
    "\n"
  )
  
  # Convert tested genes from HGNC symbols to Entrez IDs
  
  universe_conversion <- clusterProfiler::bitr(
    unique(tested_genes),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  # Create Entrez ID universe
  # Remove duplicated Entrez IDs

  tested_universe_entrez <- unique(
    universe_conversion$ENTREZID
  )
  
  # Print universe conversion statistics
  
  cat(
    "Tested genes converted to Entrez IDs:",
    length(
      tested_universe_entrez
    ),
    "\n"
  )
  
  universe_conversion_rate <- (
    length(
      tested_universe_entrez
    ) /
      length(
        unique(
          tested_genes
        )
      )
  ) * 100
  
  cat(
    "Universe conversion rate:",
    round(
      universe_conversion_rate,
      1
    ),
    "%\n"
  )
  
  # Check that enough genes remain in the universe
  
  if (
    length(
      tested_universe_entrez
    ) == 0
  ) {
    
    warning(
      paste(
        "No genes in the tested universe could be converted",
        "to Entrez IDs for cell type:",
        ct,
        "- skipping."
      )
    )
    
    next
    
  }
  
  # Identify upregulated genes
  # PCH > CTRL
  # Criteria: FDR < 0.05 & logFC > 1
  
  upregulated <- results %>%
    dplyr::filter(
      !is.na(FDR),
      !is.na(logFC),
      FDR < FDR_cutoff,
      logFC > logFC_cutoff
    )
  
  # Identify downregulated genes
  # PCH < CTRL
  # Criteria: FDR < 0.05 & logFC < -1

  downregulated <- results %>%
    dplyr::filter(
      !is.na(FDR),
      !is.na(logFC),
      FDR < FDR_cutoff,
      logFC < -logFC_cutoff
    )
  
  # Print number of DE genes

  cat(
    "\nNumber of upregulated genes:",
    nrow(upregulated),
    "\n"
  )
  
  cat(
    "Number of downregulated genes:",
    nrow(downregulated),
    "\n"
  )
  
  # Function to perform GO enrichment
  
  run_go_enrichment <- function(
    gene_data,
    direction,
    cell_type,
    tested_universe_entrez
  ) {
    
    # Check whether any genes are available

    if (
      nrow(gene_data) == 0
    ) {
      
      cat(
        "\nNo",
        direction,
        "DE genes for",
        cell_type,
        "\n"
      )
      
      return(NULL)
      
    }
    
    # Convert DE gene symbols to Entrez IDs

    gene_conversion <- clusterProfiler::bitr(
      unique(gene_data$gene),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    # Print conversion statistics

    cat(
      "\n",
      direction,
      "genes before ID conversion:",
      length(
        unique(
          gene_data$gene
        )
      ),
      "\n"
    )
    
    
    cat(
      direction,
      "genes successfully converted:",
      length(
        unique(
          gene_conversion$ENTREZID
        )
      ),
      "\n"
    )
    
    # Check whether any genes were successfully converted

    if (
      nrow(
        gene_conversion
      ) == 0
    ) {
      
      cat(
        "No genes could be converted to Entrez IDs.\n"
      )
      
      return(NULL)
      
    }
    
    # Create vector of Entrez IDs

    gene_list <- unique(
      gene_conversion$ENTREZID
    )
    
    # Make sure DE genes are part of the tested universe.
    # This is a safety check. Genes outside the tested universe should not be included in the enrichment input.

    gene_list <- intersect(
      gene_list,
      tested_universe_entrez
    )
    
    
    cat(
      direction,
      "genes within tested universe:",
      length(
        gene_list
      ),
      "\n"
    )
    
    # Check whether any DE genes remain

    if (
      length(
        gene_list
      ) == 0
    ) {
      
      cat(
        "No DE genes remain after restricting to the tested universe.\n"
      )
      
      return(NULL)
      
    }
    
    # Run GO Biological Process enrichment

    ego <- clusterProfiler::enrichGO(
      
      gene = gene_list,
      
      universe = tested_universe_entrez,
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "BP",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.05,
      
      readable = TRUE
      
    )
    
    # Convert results to data frame
    
    go_results <- as.data.frame(
      ego
    )
    
    # Create safe cell type name

    safe_ct <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      cell_type
    )

    # Determine output directory

    if (
      direction == "Upregulated"
    ) {
      
      output_dir <- file.path(
        "edgeR_DGE_results",
        "GO_Enrichment",
        "Upregulated"
      )
      
    } else {
      
      output_dir <- file.path(
        "edgeR_DGE_results",
        "GO_Enrichment",
        "Downregulated"
      )
      
    }
    
    # Save all GO enrichment results

    write.csv(
      go_results,
      file = file.path(
        output_dir,
        paste0(
          "GO_BP_",
          direction,
          "_",
          safe_ct,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    # Identify significant GO terms, FDR < 0.05
 
    significant_results <- go_results[
      !is.na(
        go_results$p.adjust
      ) &
        go_results$p.adjust < 0.05,
      ,
      drop = FALSE
    ]

    # Save significant GO terms

    write.csv(
      significant_results,
      file = file.path(
        output_dir,
        paste0(
          "GO_BP_",
          direction,
          "_",
          safe_ct,
          "_FDR05.csv"
        )
      ),
      row.names = FALSE
    )

    # Generate dot plot
    
    if (nrow(significant_results) > 0) {
      
      p_go <- clusterProfiler::dotplot(
        ego,
        showCategory = 20
      ) +
        ggplot2::ggtitle(
          paste(
            "GO Biological Process:",
            cell_type,
            "|",
            direction
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
      print(p_go)
    }

    # Print top enriched GO terms

    cat(
      "\nTop GO Biological Process terms:",
      direction,
      "\n"
    )
    
    
    if (
      nrow(
        significant_results
      ) > 0
    ) {
      
      
      top_results <- significant_results[
        order(
          significant_results$p.adjust
        ),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "p.adjust",
          "geneID"
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
        "No significant GO Biological Process terms detected.\n"
      )
      
    }

    # Print number of significant terms

    cat(
      "\nNumber of significant GO BP terms (FDR < 0.05):",
      nrow(
        significant_results
      ),
      "\n"
    )

    # Return enrichGO object

    return(
      ego
    )
    
  }

  # Run GO enrichment for upregulated genes

  go_results_all$Upregulated[[ct]] <- run_go_enrichment(
    gene_data = upregulated,
    direction = "Upregulated",
    cell_type = ct,
    tested_universe_entrez = tested_universe_entrez
  )

  # Run GO enrichment for downregulated genes

  go_results_all$Downregulated[[ct]] <- run_go_enrichment(
    gene_data = downregulated,
    direction = "Downregulated",
    cell_type = ct,
    tested_universe_entrez = tested_universe_entrez
  )
  
}

