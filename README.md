## This repo contains R scripts and `.ipynb` for spatial_analysis of MBEN data
`./scripts` main folder for scripts and notebooks, details can be found within each R script or in `.ipynb` documentaitons
  - `spatial_binning.R` # contains code to do spatial binning given custom bin size (ie side length of a square shape)
  - `spatial_compare_v2.ipynb` # cell-level analysis
  - `spatial_compare_mols_v2.ipynb` # molecule-level (or bin-level) analysis, spatial binning of transcript coordiantes
  - batch correction/integration using [Harmony](https://github.com/immunogenomics/harmony) on samples for each technology (MC, Vizgen, Xenium)
    - `spatial_integration_Vz_v2.R` # Vizgen
    - `spatial_integration_Re_v1.R` # MC
    - `spatial_integration_Xen_v1` # Xenium
