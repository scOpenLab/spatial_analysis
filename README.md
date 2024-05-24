## This repo contains scripts and notebooks for spatial_analysis of MBEN data
`./scripts` main folder for R scripts and `.ipynb`, details can be found within each R script or notebook documentation
  - `spatial_binning.R`: contains code to do spatial binning given custom bin size (ie side length of a square shape)
  - `spatial_compare_v2.ipynb`: cell-level analysis
  - `spatial_compare_mols_v2.ipynb`: molecule-level (or bin-level) analysis, spatial binning of transcript coordiantes
  - batch correction/integration using [Harmony](https://github.com/immunogenomics/harmony) on samples for each technology (MC, Vizgen, Xenium)
    - `spatial_integration_Vz_v2.R`: Vizgen
    - `spatial_integration_Re_v1.R`: MC
    - `spatial_integration_Xen_v1`: Xenium
  - `compute_nnn_expression.R`: script that extracts expression level per section and median next nearest neighbour distance among the same transcripts
  - clustering and cluster analysis:
    - `cluster_maker.R`: Performs leiden clustering on integrated samples of each technology.
    - `cluster_analysis.ipynb`: Performs plotting of UMAPs, barplots, dotplots, heatmaps and overlays for the rulting clusters.
    - `neighborhood_analysis.R`: Performs the permutation tests for the neighborhood analysis and plots heatmaps with p-values and z-scores. 
