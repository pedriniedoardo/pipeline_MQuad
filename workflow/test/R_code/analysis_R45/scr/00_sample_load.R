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

# parameters --------------------------------------------------------------
# fucntion for the palette of colors
col_fun <- colorRamp2(
  c(0, 0.2, 0.5), 
  c("#F9F9E0", "#A68B85", "#4A3631")
)

# read in the data --------------------------------------------------------
# read in the sample output from a single sample
# folder <- "../../../../results/standalone/vireo/TNBC1/"
folder <- "../../../../results/vireo/W8_24h_CSF-controls_plus_untreated_multiplexed/"

sample_donor <- read_csv(paste0(folder,"barcodes_donor_ids.csv"))
sample_variant <- read_csv(paste0(folder,"variant_donor.csv"))
ELBO_df <- read_csv(paste0(folder,"ELBO.csv"))

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
  arrange(id)

# pull the probability per clone per number of donors
df_prob <- sample_donor %>%
  select(sample_id,contains("prob")) %>%
  pivot_longer(names_to = "id",values_to = "prob",-sample_id) %>%
  separate(col = id,into = c("donor_number","clone_id"),sep = "_clone_") %>%
  mutate(donor_number = str_remove(donor_number,pattern = "prob_"))

# adjust the shape of the ELBO dataset
ELBO_df_long <- ELBO_df %>%
  mutate(id = row_number()) %>%
  pivot_longer(names_to = "donor_n",values_to = "ELBO",-id) %>%
  mutate(donor_n= str_remove(donor_n,pattern = "_ELBO"))

# exploration -------------------------------------------------------------
# plot the ELBO score 
p01 <- ELBO_df_long %>%
  ggplot(aes(x=donor_n,y=ELBO)) + geom_boxplot() + theme_bw()
ggsave("../out/")

# plot the mean allelc ration per variant
id_donor <- c("N2","N3","N4","N5")
list_mat <- lapply(id_donor, function(id){
  sample_variant %>% select(Variant_Name,contains(id)) %>%
    column_to_rownames("Variant_Name")
}) %>%
  setNames(id_donor)

# make the heatmap plot in a list
list_hm <- pmap(list(names(list_mat),list_mat), function(nm,mat){
  hm <- Heatmap(mat,
                column_title =  nm,
                col = col_fun,
                name = "Mean \nalleleic \nratio",
                show_row_dend = F,
                show_column_dend = F)
  
  return(hm)
})

list_hm2 <- lapply(list_hm, function(x){
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",
                           padding = unit(c(2, 2, 2, 40), "mm")))
  return(hm)
})

wrap_plots(list_hm2)
