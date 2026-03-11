# AIM ---------------------------------------------------------------------
# Use mitoClone2 to infer the clonal hierarchy and plot the phylogenetic tree from MQuad mitochondrial variant matrices.

# libraries ---------------------------------------------------------------
# If not installed: BiocManager::install("mitoClone2")
library(tidyverse)
library(Matrix)
library(mitoClone2)

# 1. functions & data loading ---------------------------------------------
read_mquad_mtx <- function(path_MM, path_barcode, path_variant) {
  mtx <- readMM(path_MM)
  barcodes <- read_tsv(path_barcode, col_names = FALSE, show_col_types = FALSE)$X1
  variants <- read_tsv(path_variant, col_names = FALSE, show_col_types = FALSE)$X1
  colnames(mtx) <- barcodes
  rownames(mtx) <- variants
  
  # mitoClone2 requires standard dense matrices, not sparse ones
  return(as.matrix(mtx))
}

# Define your MQuad output paths
path_mtx_AD <- "../../../../results/mquad/W8_24h_CSF-controls_plus_untreated_multiplexed/passed_ad.mtx"
path_mtx_DP <- "../../../../results/mquad/W8_24h_CSF-controls_plus_untreated_multiplexed/passed_dp.mtx"
path_barcodes <- "../../../../results/preprocess/W8_24h_CSF-controls_plus_untreated_multiplexed/outs/barcodes.tsv"
path_variants <- "../../../../results/mquad/W8_24h_CSF-controls_plus_untreated_multiplexed/passed_variant_names.txt"

print("Loading matrices...")
mtx_AD <- read_mquad_mtx(path_mtx_AD, path_barcodes, path_variants)
mtx_DP <- read_mquad_mtx(path_mtx_DP, path_barcodes, path_variants)

# 2. Initialize mitoClone2 Object -----------------------------------------
print("Creating mutationCalls object...")
# We map the Alternative Depth (mutants) and Total Depth directly
mut_calls <- mutationCallsFromMatrix(
  M = mtx_AD,
  N = mtx_DP
)

# run SCITE
mut_calls <- varCluster(mut_calls, tempfolder = "../out/object/SCITE",method='SCITE')


mut_calls <- clusterMetaclones(mut_calls, min.lik = 1)
plotClones(mut_calls)

m2c <- mitoClone2:::getMut2Clone(P1)
print(m2c)
