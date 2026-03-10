# AIM ---------------------------------------------------------------------
# the aim of the script is to load the sample output from the MQuand to explore the results

# renv integration --------------------------------------------------------
# to load the packages
# source(".Rprofile")

# libraries ---------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)
library(patchwork)
library(magick)
library(Matrix)
library(pals)

# functions ---------------------------------------------------------------
# funcitno to load the matrix file from MQuad
read_mquad_mtx <- function(path_MM,path_barcode,path_variant) {
  # 1. Load the matrix
  mtx <- readMM(path_MM)
  
  # 2. Load Cell Barcodes (Columns)
  barcodes <- read_tsv(path_barcode, col_names = FALSE)$X1
  
  # 3. Load Variants (Rows)
  variants <- read_tsv(path_variant, col_names = FALSE)$X1
  
  # 4. Assign names
  colnames(mtx) <- barcodes
  rownames(mtx) <- variants
  
  return(mtx)
}

# parameters --------------------------------------------------------------
# fucntion for the palette of colors
col_fun_alRatio <- colorRamp2(
  c(0, 0.2, 0.5), 
  c("#F9F9E0", "#A68B85", "#4A3631")
)

col_fun_donProb <- colorRamp2(c(0, 0.8), c("#F7FBFF", "#42648B"))

# set seed for the clustering
set.seed(2144)

# read in the data --------------------------------------------------------
# read in the sample output from a single sample
# folder <- "../../../../results/standalone/vireo/TNBC1/"
folder <- "../../../../results/vireo/s59/"

sample_donor <- read_csv(paste0(folder,"barcodes_donor_ids.csv"))
sample_variant <- read_csv(paste0(folder,"variant_donor.csv"))
ELBO_df <- read_csv(paste0(folder,"ELBO.csv"))

# matrices for variant per cell
path_mtx_AD <- "../../../../results/mquad/s59/passed_ad.mtx"
path_mtx_DP <- "../../../../results/mquad/s59/passed_dp.mtx"

# cell barcodes
path_barcodes <- "../../../../results/preprocess/s59/outs/barcodes.tsv"
path_variants <- "../../../../results/mquad/s59/passed_variant_names.txt"

# read in the matrices
mtx_AD <- read_mquad_mtx(path_MM = path_mtx_AD,
                         path_barcode = path_barcodes,
                         path_variant = path_variants)

mtx_DP <- read_mquad_mtx(path_MM = path_mtx_DP,
                         path_barcode = path_barcodes,
                         path_variant = path_variants)

# define the range of donors number used in the analysis
id_donor <- c("N2","N3","N4","N5")

# wrangling ---------------------------------------------------------------
# shape the table from wide format to long format
df_clone_id <- sample_donor %>%
  select(sample_id,contains("clone_id")) %>%
  pivot_longer(names_to = "id",values_to = "clone_id",-sample_id) %>%
  mutate(id=str_remove(id,pattern = "clone_id_"))

df_confidence <- sample_donor %>%
  select(sample_id,contains("confident")) %>%
  pivot_longer(names_to = "id",values_to = "confident",-sample_id) %>%
  mutate(id=str_remove(id,pattern = "confident_"))

# join the tables for the clone id with its confidence ad different donors resolution
df_clone <- purrr::reduce(list(df_clone_id,df_confidence),.f = left_join,by=c("sample_id","id")) %>%
  mutate(clone_id_plot = paste0("clone_",clone_id)) %>%
  arrange(id)

# pull the probability per clone per number of donors
df_prob <- sample_donor %>%
  select(sample_id,contains("prob")) %>%
  pivot_longer(names_to = "id",values_to = "prob",-sample_id) %>%
  separate(col = id,into = c("donor_number","clone_id"),sep = "_clone_") %>%
  mutate(donor_number = str_remove(donor_number,pattern = "prob_")) %>%
  mutate(clone_id_plot = paste0("clone_",clone_id))

# adjust the shape of the ELBO dataset
ELBO_df_long <- ELBO_df %>%
  mutate(id = row_number()) %>%
  pivot_longer(names_to = "donor_n",values_to = "ELBO",-id) %>%
  mutate(donor_n= str_remove(donor_n,pattern = "_ELBO"))

# caluclate the allele frequency per barcode
# Avoid division by zero by setting 0/0 to 0
af_matrix <- as.matrix(mtx_AD/mtx_DP)
af_matrix[is.na(af_matrix)] <- 0

# exploration -------------------------------------------------------------
# plot the ELBO score 
p01 <- ELBO_df_long %>%
  ggplot(aes(x=donor_n,y=ELBO)) + geom_boxplot() + theme_bw()
ggsave(plot = p01,filename = "../out/plot/00_ELBO.pdf",width = 4,height = 3)

# plot the mean allelic ration per variant
list_mat <- lapply(id_donor, function(id){
  sample_variant %>% select(Variant_Name,contains(id)) %>%
    column_to_rownames("Variant_Name")
}) %>%
  setNames(id_donor)

# make the heatmap plot in a list
list_hm <- pmap(list(names(list_mat),list_mat), function(nm,mat){
  hm <- Heatmap(mat,
                column_title =  nm,
                col = col_fun_alRatio,
                name = "Mean \nalleleic \nratio",
                show_row_dend = F,
                show_column_dend = F)
  
  return(hm)
})

list_plothm <- lapply(list_hm, function(x){
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",
                           padding = unit(c(2, 2, 2, 40), "mm")))
  return(hm)
})

p02 <- wrap_plots(list_plothm)
ggsave(plot = p02,filename = "../out/plot/00_hm_allelic_ratio.pdf",width = 12,height = 12)

# plot the probability of the assignament per cell to a clone id
list_mat2 <- lapply(id_donor, function(id){
  df_prob %>%
    filter(donor_number %in% id) %>%
    select(sample_id,prob,clone_id_plot) %>%
    pivot_wider(names_from = clone_id_plot,values_from = prob) %>%
    column_to_rownames("sample_id")
}) %>%
  setNames(id_donor)

list_hm2 <- pmap(list(names(list_mat2),list_mat2), function(nm,mat){
  hm <- Heatmap(mat,
                column_title =  nm,
                # --- ROW AGGREGATION SETTINGS ---
                row_km = 20,              # Aggregates 28k rows into 20 clusters
                row_km_repeats = 1,        # Faster clustering
                row_title = "Clustered Cells", 
                # --------------------------------
                show_row_names = F,
                col = col_fun_donProb,
                name = "Clone \nprob ",
                show_row_dend = F,
                show_column_dend = F,
                use_raster = TRUE,
                raster_quality = 2)
  
  return(hm)
})

list_plothm2 <- lapply(list_hm2, function(x){
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",
                           padding = unit(c(2, 2, 2, 40), "mm")))
  return(hm)
})

p03 <- wrap_plots(list_plothm2)
ggsave(plot = p03,filename = "../out/plot/00_hm_donor_prob.pdf",width = 12,height = 12)

# plot the frequency of the allele by clone
# all the barcodes
# id_don <- "N3"
list_mat3_all <- lapply(id_donor, function(id_don){
  # order the heatmap by clone id
  df_clone_id <- df_clone %>%
    filter(id == id_don)
  
  # save the meta with all the barcodes and clone id
  sample_id_all <- df_clone_id %>%
    arrange(clone_id_plot)
  
  # reorder the matrix
  af_matrix_orderAll <- af_matrix[,sample_id_all %>% pull(sample_id),drop = FALSE]
  
  # define the column annotation
  cluster_colors <- alphabet(length(unique(df_clone_id$clone_id_plot)))
  names(cluster_colors) <- unique(unique(df_clone_id$clone_id_plot))
  
  column_ha_all <- HeatmapAnnotation(
    cluster = sample_id_all %>% pull(clone_id_plot),     # Rotates to 45 degrees
    col = list(cluster = cluster_colors)
    )
  
  return(list(mat = af_matrix_orderAll,
              meta = sample_id_all,
              clust = column_ha_all))
}) %>%
  setNames(id_donor)

list_hm3_all <- pmap(list(names(list_mat3_all),list_mat3_all), function(nm,obj){
  # plot the heatamp
  hm <- Heatmap(
    matrix = obj$mat, 
    # column_title =  nm,
    name = "AF", 
    top_annotation = obj$clust,
    col = col_fun_alRatio,
    
    # Clustering & Ordering
    cluster_columns = F,
    column_split = obj$meta %>% pull(clone_id_plot),
    
    # Visuals
    show_row_names = T,
    show_column_names = F,
    
    # rotate the lales
    column_title_rot = 45,
    column_title_side = "top", 
    
    # Performance for 28k cells
    use_raster = T,
    raster_quality = 2
  )
  return(hm)
})

list_plothm3_all <- lapply(list_hm3_all, function(x){
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",
                           padding = unit(c(2, 2, 2, 40), "mm")))
  return(hm)
})

p04 <- wrap_plots(list_plothm3_all)
ggsave(plot = p04,filename = "../out/plot/00_hm_donor_AF_all.pdf",width = 20,height = 12)

# only the confidet ones
list_mat3_confident <- lapply(id_donor, function(id_don){
  # order the heatmap by clone id
  df_clone_id <- df_clone %>%
    filter(id == id_don)
  
  # save the meta with anly confidend assignaments and clone id
  sample_id_confident <- df_clone_id %>%
    filter(confident == T) %>%
    arrange(clone_id_plot)
  
  # order the matrix
  af_matrix_orderConfident <- af_matrix[,sample_id_confident %>% pull(sample_id),drop = FALSE]
  
  # Define the Column Annotation (The colored bars at the top)
  cluster_colors <- alphabet(length(unique(df_clone_id$clone_id_plot)))
  names(cluster_colors) <- unique(unique(df_clone_id$clone_id_plot))
  
  column_ha_confident <- HeatmapAnnotation(
    cluster = sample_id_confident %>% pull(clone_id_plot),
    col = list(cluster = cluster_colors)
  )
  
  return(list(mat = af_matrix_orderConfident,
              meta = sample_id_confident,
              clust = column_ha_confident))
}) %>%
  setNames(id_donor)

list_hm3_confident <- pmap(list(names(list_mat3_confident),list_mat3_confident), function(nm,obj){
  # plot the heatamp
  hm <- Heatmap(
    matrix = obj$mat, 
    # column_title =  nm,
    name = "AF", 
    top_annotation = obj$clust,
    col = col_fun_alRatio,
    
    # Clustering & Ordering
    cluster_columns = FALSE, # Manually ordered by your data prep
    column_split = obj$meta %>% pull(clone_id_plot), # This creates the gaps between clones
    
    # Visuals
    show_row_names = TRUE,
    show_column_names = FALSE,
    
    # rotate the lales
    column_title_rot = 45,
    column_title_side = "top", 
    
    # Performance for 28k cells
    use_raster = TRUE,
    raster_quality = 2
  )
  return(hm)
})

list_plothm3_confident <- lapply(list_hm3_confident, function(x){
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",
                           padding = unit(c(2, 2, 2, 40), "mm")))
  return(hm)
})

p05 <- wrap_plots(list_plothm3_confident)
ggsave(plot = p05,filename = "../out/plot/00_hm_donor_AF_confident.pdf",width = 20,height = 12)
