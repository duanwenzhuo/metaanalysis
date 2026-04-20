#!/usr/bin/env Rscript
# install_r_packages.R
# Install R packages needed for bioskills differential expression analysis

required_packages <- c(
  "DESeq2",        # Differential expression
  "edgeR",         # Differential expression  
  "limma",         # Linear models for DE
  "sva",           # Surrogate variable analysis
  "clusterProfiler", # Enrichment analysis
  "org.Hs.eg.db",  # Human gene annotation
  "org.Mm.eg.db",  # Mouse gene annotation
  "AnnotationDbi"
)

message("Installing R packages for MetaAnalysisi...")
message("This may take several minutes...\n")

# Detect if packages are already installed
missing <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing) == 0) {
  message("All packages already installed!")
} else {
  message(paste("Installing", length(missing), "missing package(s):", paste(missing, collapse=", ")))
  BiocManager::install(missing, update = FALSE, ask = FALSE)
}

message("\n✅ R packages installation complete!")
