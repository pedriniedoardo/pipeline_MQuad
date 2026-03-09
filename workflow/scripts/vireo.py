# --- Load Libraries ---
import vireoSNP
from vireoSNP import BinomMixtureVB
import numpy as np
from scipy.io import mmread
import pandas as pd

print("vireoSNP version:", vireoSNP.__version__)

# --- Snakemake Integration ---
# Inputs
in_passed_ad     = snakemake.input.passed_ad
in_passed_dp     = snakemake.input.passed_dp
in_barcodes      = snakemake.input.barcodes
in_variant_names = snakemake.input.variant_names

# Outputs
out_clones_df  = snakemake.output.clones_df
out_variant_df = snakemake.output.variant_df
out_ELBO_df = snakemake.output.ELBO_df

# Parameters
n_clones_to_test = snakemake.params.n_donors

# --- 1. Load Data ---
AD = mmread(in_passed_ad).tocsc()
DP = mmread(in_passed_dp).tocsc()

with open(in_barcodes, 'r') as f:
    sample_id = f.read().splitlines()

with open(in_variant_names, 'r') as f:
    variant_names = f.read().splitlines()


# --- 2. Initialize Master DataFrames ---
clones_df = pd.DataFrame({'sample_id': sample_id})
variant_df = pd.DataFrame(index=variant_names)
ELBO_df = pd.DataFrame()
variant_df.index.name = "Variant_Name"


# --- 3. Iterate over the range of clones ---
for n_clones in n_clones_to_test:
    print(f"Running vireoSNP for N={n_clones} clones...")
    
    # Initialize and fit the model
    _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=n_clones)
    _model.fit(AD, DP, min_iter=30, n_init=300)
    
    # --- Process Cells ---
    prob_matrix = _model.ID_prob
    clone_id = np.argmax(prob_matrix, axis=1)
    conf = np.max(prob_matrix, axis=1) >= 0.8
    
    clones_df[f'clone_id_N{n_clones}'] = clone_id
    clones_df[f'confident_N{n_clones}'] = conf
    
    for i in range(n_clones):
        clones_df[f'prob_N{n_clones}_clone_{i}'] = prob_matrix[:, i]
        
    # --- Process Variants ---
    beta_mu_matrix = _model.beta_mu
    for i in range(n_clones):
        variant_df[f'N{n_clones}_Clone_{i}_VAF'] = beta_mu_matrix[:, i]
    
    # Extract the 300 initialization scores and save them as a new column
    ELBO_df[f'N{n_clones}_ELBO'] = _model.ELBO_inits


# --- 4. Save Final Master Tables ---
clones_df.to_csv(out_clones_df, index=False)
variant_df.to_csv(out_variant_df)
ELBO_df.to_csv(out_ELBO_df, index=False)

print(f"\nSuccessfully saved {out_clones_df}, {out_variant_df}, and {out_ELBO_df}")