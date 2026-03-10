message("--- Starting renv test ---")

# 1. Activate the renv environment using your absolute path
message("Activating renv...")
# renv integration --------------------------------------------------------
# to load the packages
# source(".Rprofile")
renv::load("workflow/test/R_code/analysis_R45")

# libraries ---------------------------------------------------------------
library(tidyverse)

# -------------------------------------------------------------------------
mtcars %>%
  filter(carb==4)

# 3. Create the dummy output file so Snakemake marks the rule as complete
file.create(snakemake@output[[1]])

message("--- Test finished successfully ---")