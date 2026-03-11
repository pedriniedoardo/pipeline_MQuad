# libraries ---------------------------------------------------------------
library(tidyverse)
# install.packages("ape")
library(ape)        
library(ggtree)     


# read in the data --------------------------------------------------------
# load the output of vireo for the variant frequency per donor
path_variant <- "../../../../results/vireo/W8_24h_CSF-controls_plus_untreated_multiplexed/variant_donor.csv"
df_variants <- read_csv(path_variant)

# load the table with clone assignament
path_clone <- "../../../../results/vireo/SC5v2_Melanoma_5Kcells_Connect_single_channel_gex/barcodes_donor_ids.csv"
df_clone <- read_csv(path_clone)

# Select the specific resolution you want to plot (e.g., N=4 clones)
target_N <- "N4"

# wrangling ---------------------------------------------------------------
# Extract just those columns and transpose the matrix
clone_matrix <- df_variants %>% 
  select(Variant_Name, contains(target_N)) %>% 
  column_to_rownames("Variant_Name")

# define the number of cell per cluster
column_id <- paste0("clone_id_",target_N)
newname <- "label"

# build the summary matrix
df_summary_clones <- df_clone %>%
  select(sample_id, contains(column_id)) %>%
  dplyr::rename_with(~newname, .cols = column_id) %>%
  group_by(label) %>%
  summarise(n = n()) %>%
  mutate(label = paste0("Clone_",label)) %>%
  # add the simulated ancestor reference
  bind_rows(data.frame(label = "Ancestor (sim)", n = 1))

# Transpose so Clones are ROWS and Variants are COLUMNS
mat_clones <- t(clone_matrix)

# Clean up the names (change "N4_clone_0" to just "Clone_0")
rownames(mat_clones) <- str_replace_all(rownames(mat_clones), pattern = paste0(target_N, "_|_VAF"), replacement = "")

# 3. Add the Ancestral Root 
# We create a fake reference cell that has 0% VAF for every mutation. This forces the tree to have a clear starting point (the root).
root_node <- rep(0, ncol(mat_clones))
mat_clones_rooted <- rbind(mat_clones, 'Ancestor (sim)' = root_node)

# 4. Calculate Distance and Build the Tree
# Calculate how different each clone is from the others based on their mutation frequencies
dist_mat <- dist(mat_clones_rooted, method = "euclidean")

# Use Neighbor-Joining (NJ) to build the evolutionary tree
clone_tree <- nj(dist_mat)

# Anchor the tree at root "Ancestor"
clone_tree <- root(clone_tree, outgroup = "Ancestor (sim)", resolve.root = TRUE)

# -------------------------------------------------------------------------
# force the anchestro into the root
# 1. Find the internal ID of the Ancestor tip
anc_index <- which(clone_tree$tip.label == "Ancestor (sim)")
# 2. Find the branch (edge) leading to that tip
anc_edge <- which(clone_tree$edge[, 2] == anc_index)
# 3. Force the length of that branch to exactly 0
clone_tree$edge.length[anc_edge] <- 0
# -------------------------------------------------------------------------

# This creates a 'treedata' object that securely holds the branches AND your cell counts.
tree_with_data <- full_join(clone_tree, df_summary_clones, by = "label")

# plotting ----------------------------------------------------------------
# Plot the Clone Tree
p_tree <- ggtree(tree_with_data, aes(color = label),size = 1) +
  geom_tippoint(aes(color = label, size = n)) +
  geom_tiplab(size = 5, fontface = "bold", offset = 0.1) +
  theme_tree2() +
  ggtitle("Evolutionary Lineage of N4 Clones") +
  coord_cartesian(clip = "off") + 
  scale_size_continuous(range = c(1, 8), name = "Number of Cells",breaks = c(10,100,1000,2000)) +
  theme(
    legend.position = "right",
    plot.margin = margin(t = 20, r = 150, b = 20, l = 20, unit = "pt")
  ) +
  guides(color = "none")

p_tree
# ggsave("../out/plot/00_Clone_Lineage_Tree.pdf", plot = p_tree, width = 6, height = 5)