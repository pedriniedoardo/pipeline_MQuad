# AIM ---------------------------------------------------------------------
# Run cardelino de-novo to assign cells to clones based on MQuad output

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Matrix)
library(cardelino)  # BiocManager::install("cardelino")

# functions ---------------------------------------------------------------
read_mquad_mtx <- function(path_MM, path_barcode, path_variant) {
  mtx <- readMM(path_MM)
  barcodes <- read_tsv(path_barcode, col_names = FALSE, show_col_types = FALSE)$X1
  variants <- read_tsv(path_variant, col_names = FALSE, show_col_types = FALSE)$X1
  colnames(mtx) <- barcodes
  rownames(mtx) <- variants
  return(mtx)
}

# 1. read in the data -----------------------------------------------------
# Update these paths to match your directory structure
path_mtx_AD <- "../../../../results/mquad/W8_24h_CSF-controls_plus_untreated_multiplexed/passed_ad.mtx"
path_mtx_DP <- "../../../../results/mquad/W8_24h_CSF-controls_plus_untreated_multiplexed/passed_dp.mtx"
path_barcodes <- "../../../../results/preprocess/W8_24h_CSF-controls_plus_untreated_multiplexed/outs/barcodes.tsv"
path_variants <- "../../../../results/mquad/W8_24h_CSF-controls_plus_untreated_multiplexed/passed_variant_names.txt"

mtx_AD <- read_mquad_mtx(path_mtx_AD, path_barcodes, path_variants)
mtx_DP <- read_mquad_mtx(path_mtx_DP, path_barcodes, path_variants)

# 2. run cardelino inference ----------------------------------------------
print("Starting cardelino Gibbs sampling... This may take a few minutes!")

# Set a seed so your results are perfectly reproducible
set.seed(42)

# Run clone_id in de novo mode (Config = NULL) for 3 clones
# Note: Gibbs sampling is computationally heavy. If you have >5,000 cells,
# keep min_iter and max_iter relatively low for your first test run.
assignments <- clone_id(
  A = mtx_AD, 
  D = mtx_DP, 
  Config = NULL, 
  n_clone = 3, 
  min_iter = 800, 
  max_iter = 1200
)

# 3. extract and wrangle results ------------------------------------------
# cardelino has a built-in function to assign cells based on a probability threshold (e.g., 0.5 or 0.8)
df_clones <- assign_cells_to_clones(assignments$prob, threshold = 0.8)

# Add the cell barcodes back into the dataframe for easy Seurat/Scanpy integration
df_clones <- df_clones %>%
  mutate(sample_id = colnames(mtx_AD)) %>%
  relocate(sample_id)

print(head(df_clones))

# Save the assignments
# write_csv(df_clones, "../out/cardelino_N3_assignments.csv")

# 4. visualization --------------------------------------------------------
# Calculate raw allele frequency matrix (filling NAs with 0)
AF <- as.matrix(mtx_AD / mtx_DP)
AF[is.na(AF)] <- 0

# Use cardelino's built-in variant-clone heatmap plotter
# It automatically clusters the cells and shows the inferred Config (variant signature)
p_cardelino_hm <- vc_heatmap(
  mat = AF, 
  prob = assignments$prob, 
  Config = assignments$Config_est, 
  show_legend = TRUE
)

# You can save this plot just like your previous ones
# pdf("../out/plot/cardelino_heatmap.pdf", width = 10, height = 8)
# print(p_cardelino_hm)
# dev.off()