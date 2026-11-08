# Advanced data analysis project

Stats and Data Science PhD requirement. 

## Data source

Data cleaning and analysis of Lee et al., 2024, available on [Zenodo](https://zenodo.org/records/10888270).

## Notebooks

Brief summaries of the contents of `notebooks/`, reverse chronologically

- `20261028`: Inter-region site sparsity experiments
  - Sparsity calculated on log ratios
- `20261027`: Single site texture discrimination
  - Cohen's d averaged over pairwise families
- `20261007`: Inter-region site sparsity
  - Initial definitions of sparsity (ratio of moments, Gini)
  - Baseline correction with truncation or absolute values
  - Sparsity and split-half correlation (with or without Spearman-Brown correction)
- `20260914`: verification that performance is stable in the most reliable subsets of sites
  - Since IT needs to be omitted, these only verify the relative performance of V1 vs. V2
- `20260909`: demonstrating that pairwise Mahalanobis distance is comparable to results from randomized dichotomies
  - Visual comparison and correlation of family-wise means
  - This will help iterate through sample-wise results faster
- `20260826`: returning to reliable/responsive site characterization
  - New onset time calculations
  - Matched reliability with matching interval of Pearson's correlation
- `20260820`: Normalization of neural responses, Mahalanobis robustness
- `20260813`: Detection of sparsity in neural responses
- `20260723`: Randomized dichotomy within-family d' calculations 
- `20260715`: Pairwise between-family d' calculations
- `20260625`: Onset time and responsive site calculation, regression models
  - Useful reference for data loading
- `20260611`: Split-half PCA reliability
- `20260604`: PSTHs and responsive sites on original 4 repetition data
- `20260521`: redo of MATLAB loading
- `20260303`, `20260401`: discarded data loading procedures in Python
- `20260212`: PSTH sanity check for noisy/naturalistic images
- `20260206`: ADA proposed abstract
- `20260201`: extracting some info for ADA 1.0
- `20260122`: adding some start time extraction
- `20260107`: initial check of eigenvalues from PCA (not SVD)
- `20251120`: initial loading of MATLAB files
