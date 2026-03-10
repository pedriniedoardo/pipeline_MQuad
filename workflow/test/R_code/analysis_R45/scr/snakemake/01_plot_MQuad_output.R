# AIM ---------------------------------------------------------------------
# The aim of the script is to load the sample output from MQuad/vireo and explore the results via ComplexHeatmap.

# --- Snakemake Integration ---

# Inputs
in_donor <- snakemake@input[["clone_assignments"]]
# in_donor <- "../../../../results/vireo/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/barcodes_donor_ids.csv"
in_variant <- snakemake@input[["variant_donors"]]
# in_variant <- "../../../../results/vireo/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/variant_donor.csv"
in_elbo <- snakemake@input[["elbo_inits"]]
# in_elbo <- "../../../../results/vireo/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/ELBO.csv"
in_mtx_AD <- snakemake@input[["passed_ad"]]
# in_mtx_AD <- "../../../../results/mquad/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/passed_ad.mtx"
in_mtx_DP <- snakemake@input[["passed_dp"]]
# in_mtx_DP <- "../../../../results/mquad/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/passed_dp.mtx"
in_barcodes <- snakemake@input[["barcodes"]]
# in_barcodes <- "../../../../results/preprocess/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/outs/barcodes.tsv"
in_variant_name <- snakemake@input[["variant_names"]]
# in_variant_name <- "../../../../results/mquad/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/passed_variant_names.txt"

# Outputs
out_plot_elbo <- snakemake@output[["plot_elbo"]]
out_plot_ratio <- snakemake@output[["plot_ratio"]]
out_plot_prob <- snakemake@output[["plot_prob"]]
out_plot_af_all <- snakemake@output[["plot_af_all"]]
out_plot_af_conf <- snakemake@output[["plot_af_conf"]]

# Parameters
id_donor <- paste0("N", snakemake@params[["donors"]])
# id_donor <- c("N2","N3","N4","N5")

# renv integration --------------------------------------------------------
# to load the packages
# source(".Rprofile")
renv::load("workflow/test/R_code/analysis_R45")

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)
library(patchwork)
library(magick)
library(Matrix)
library(pals)

# functions ---------------------------------------------------------------
read_mquad_mtx <- function(path_MM, path_barcode, path_variant) {
  mtx <- readMM(path_MM)
  barcodes <- read_tsv(path_barcode, col_names = FALSE, show_col_types = FALSE)$X1
  variants <- read_tsv(path_variant, col_names = FALSE, show_col_types = FALSE)$X1
  colnames(mtx) <- barcodes
  rownames(mtx) <- variants
  return(mtx)
}

# parameters --------------------------------------------------------------
# define the colors palette
col_fun_alRatio <- colorRamp2(
  c(0, 0.2, 0.5), 
  c("#F9F9E0", "#A68B85", "#4A3631")
)

col_fun_donProb <- colorRamp2(c(0, 0.8), c("#F7FBFF", "#42648B"))

# set seed for the clustering
set.seed(2144)

# read in the data --------------------------------------------------------
# read in the barcodes per assigned donor id
sample_donor <- read_csv(in_donor, show_col_types = FALSE)
# read in the AF table 
sample_variant <- read_csv(in_variant, show_col_types = FALSE)
# red in the ELBO plot
ELBO_df <- read_csv(in_elbo, show_col_types = FALSE)

# redad in the matrices for Allele counts
mtx_AD <- read_mquad_mtx(path_MM = in_mtx_AD, path_barcode = in_barcodes, path_variant = in_variant_name)
mtx_DP <- read_mquad_mtx(path_MM = in_mtx_DP, path_barcode = in_barcodes, path_variant = in_variant_name)

# wrangling ---------------------------------------------------------------
# define the donor id per barcode per normber of donor
df_clone_id <- sample_donor %>%
  select(sample_id, contains("clone_id")) %>%
  pivot_longer(names_to = "id", values_to = "clone_id", -sample_id) %>%
  mutate(id = str_remove(id, pattern = "clone_id_"))

# define the confidence of the donor assignamebt
df_confidence <- sample_donor %>%
  select(sample_id, contains("confident")) %>%
  pivot_longer(names_to = "id", values_to = "confident", -sample_id) %>%
  mutate(id = str_remove(id, pattern = "confident_"))

# merge the information
df_clone <- purrr::reduce(list(df_clone_id, df_confidence), .f = left_join, by = c("sample_id", "id")) %>%
  mutate(clone_id_plot = paste0("clone_", clone_id)) %>%
  arrange(id)

# probabiluty of the donor assignament per barcode
df_prob <- sample_donor %>%
  select(sample_id, contains("prob")) %>%
  pivot_longer(names_to = "id", values_to = "prob", -sample_id) %>%
  separate(col = id, into = c("donor_number", "clone_id"), sep = "_clone_") %>%
  mutate(donor_number = str_remove(donor_number, pattern = "prob_")) %>%
  mutate(clone_id_plot = paste0("clone_", clone_id))

# ELBO plot
ELBO_df_long <- ELBO_df %>%
  mutate(id = row_number()) %>%
  pivot_longer(names_to = "donor_n", values_to = "ELBO", -id) %>%
  mutate(donor_n = str_remove(donor_n, pattern = "_ELBO"))

# calculate the matrix od AF
af_matrix <- as.matrix(mtx_AD / mtx_DP)
af_matrix[is.na(af_matrix)] <- 0

# exploration -------------------------------------------------------------
# 1. Plot ELBO
p01 <- ELBO_df_long %>%
  ggplot(aes(x = donor_n, y = ELBO)) + geom_boxplot() + theme_bw()
ggsave(plot = p01, filename = out_plot_elbo, width = 4, height = 3)

# 2. Plot Allelic Ratio
list_mat <- lapply(id_donor, function(id){
  sample_variant %>% select(Variant_Name, contains(id)) %>%
    column_to_rownames("Variant_Name")
}) %>% setNames(id_donor)

list_hm <- pmap(list(names(list_mat), list_mat), function(nm, mat){
  Heatmap(mat, column_title = nm, col = col_fun_alRatio, name = "Mean \nallelic \nratio",
          show_row_dend = FALSE, show_column_dend = FALSE)
})

list_plothm <- lapply(list_hm, function(x){
  grid.grabExpr(draw(x, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 40), "mm")))
})
p02 <- wrap_plots(list_plothm)
ggsave(plot = p02, filename = out_plot_ratio, width = 12, height = 12)

# 3. Plot Donor Probabilities
list_mat2 <- lapply(id_donor, function(id){
  df_prob %>% filter(donor_number %in% id) %>%
    select(sample_id, prob, clone_id_plot) %>%
    pivot_wider(names_from = clone_id_plot, values_from = prob) %>%
    column_to_rownames("sample_id")
}) %>% setNames(id_donor)

list_hm2 <- pmap(list(names(list_mat2), list_mat2), function(nm, mat){
  Heatmap(mat, column_title = nm, row_km = 20, row_km_repeats = 1, row_title = "Clustered Cells",
          show_row_names = FALSE, col = col_fun_donProb, name = "Clone \nprob",
          show_row_dend = FALSE, show_column_dend = FALSE, use_raster = TRUE, raster_quality = 2)
})

list_plothm2 <- lapply(list_hm2, function(x){
  grid.grabExpr(draw(x, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 40), "mm")))
})
p03 <- wrap_plots(list_plothm2)
ggsave(plot = p03, filename = out_plot_prob, width = 12, height = 12)

# 4. Plot AF All
list_mat3_all <- lapply(id_donor, function(id_don){
  df_clone_id <- df_clone %>% filter(id == id_don)
  sample_id_all <- df_clone_id %>% arrange(clone_id_plot)
  # avoid conversion from matrix to vectro in case of a single variant
  af_matrix_orderAll <- af_matrix[, sample_id_all %>% pull(sample_id),drop = FALSE]
  
  cluster_colors <- alphabet(length(unique(df_clone_id$clone_id_plot)))
  names(cluster_colors) <- unique(unique(df_clone_id$clone_id_plot))
  
  column_ha_all <- HeatmapAnnotation(cluster = sample_id_all %>% pull(clone_id_plot), col = list(cluster = cluster_colors))
  return(list(mat = af_matrix_orderAll, meta = sample_id_all, clust = column_ha_all))
}) %>% setNames(id_donor)

list_hm3_all <- pmap(list(names(list_mat3_all), list_mat3_all), function(nm, obj){
  Heatmap(matrix = obj$mat,
          # column_title = nm,
          name = "AF", top_annotation = obj$clust, col = col_fun_alRatio,
          cluster_columns = FALSE, column_split = obj$meta %>% pull(clone_id_plot), show_row_names = TRUE,
          show_column_names = FALSE, column_title_rot = 45, column_title_side = "top", use_raster = TRUE, raster_quality = 2)
})

list_plothm3_all <- lapply(list_hm3_all, function(x){
  grid.grabExpr(draw(x, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 40), "mm")))
})
p04 <- wrap_plots(list_plothm3_all)
ggsave(plot = p04, filename = out_plot_af_all, width = 20, height = 12)

# 5. Plot AF Confident
list_mat3_confident <- lapply(id_donor, function(id_don){
  df_clone_id <- df_clone %>% filter(id == id_don)
  sample_id_confident <- df_clone_id %>% filter(confident == TRUE) %>% arrange(clone_id_plot)
  # avoid conversion from matrix to vectro in case of a single variant
  af_matrix_orderConfident <- af_matrix[, sample_id_confident %>% pull(sample_id),drop = FALSE]
  
  cluster_colors <- alphabet(length(unique(df_clone_id$clone_id_plot)))
  names(cluster_colors) <- unique(unique(df_clone_id$clone_id_plot))
  
  column_ha_confident <- HeatmapAnnotation(cluster = sample_id_confident %>% pull(clone_id_plot), col = list(cluster = cluster_colors))
  return(list(mat = af_matrix_orderConfident, meta = sample_id_confident, clust = column_ha_confident))
}) %>% setNames(id_donor)

list_hm3_confident <- pmap(list(names(list_mat3_confident), list_mat3_confident), function(nm, obj){
  Heatmap(matrix = obj$mat,
          # column_title = nm,
          name = "AF", top_annotation = obj$clust, col = col_fun_alRatio,
          cluster_columns = FALSE, column_split = obj$meta %>% pull(clone_id_plot), show_row_names = TRUE,
          show_column_names = FALSE, column_title_rot = 45, column_title_side = "top", use_raster = TRUE, raster_quality = 2)
})

list_plothm3_confident <- lapply(list_hm3_confident, function(x){
  grid.grabExpr(draw(x, heatmap_legend_side = "left", padding = unit(c(2, 2, 2, 40), "mm")))
})
p05 <- wrap_plots(list_plothm3_confident)
ggsave(plot = p05, filename = out_plot_af_conf, width = 20, height = 12)