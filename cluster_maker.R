rm(list = ls())

library(Seurat)
library(leiden)
library(igraph)

# https://satijalab.org/seurat/articles/future_vignette.html
options(future.globals.maxSize = 120000 * 1024^2) # ~125G RAM

# Set working drirectory
wd <- ""
setwd(wd)

# Read objects:
message("Reading Objects...")
mc_matched_hm <- readRDS("./objects/segmentation_based/mc_matched_hm.rds") # only Molecular Cartography
vz_matched_hm <- readRDS("./objects/segmentation_based/vz_matched_hm.rds") # only Merscope
xen_matched_hm <- readRDS("./objects/segmentation_based/xen_matched_hm.rds") # only Xenium

# Leiden CLustering function
cluster_maker <- function(obj, resolutions){
    # Remove any previous clustering results
	cols_to_remove <- c("SCT_snn_res.0.1.leiden", "SCT_snn_res.0.2.leiden",
						"SCT_snn_res.0.1.leiden.hm", "SCT_snn_res.0.2.leiden.hm")
    obj@meta.data <- obj@meta.data[ , !(names(obj@meta.data) %in% cols_to_remove)]
    
    # Removed any previous clustering results
	obj <- FindClusters(obj, resolution = resolutions,
			method = "igraph", 
			algorithm = 4) # leiden clustering
	return(obj)
}

# Setup Parallelization
message("Setting up parallelization...")
plan("multicore", workers = 12)
plan()

# Make the clusters, vary the resolution parameter from 0 to 1 with 0.05 increments
message("Making the Clusters...")
resolutions <- seq(from = 0.1, to = 1, by = 0.05)
message("MC...")
mc_matched_hm <- cluster_maker(mc_matched_hm, resolutions)
message("Merscope...")
vz_matched_hm <- cluster_maker(vz_matched_hm, resolutions)
message("Xenium...")
xen_matched_hm <- cluster_maker(xen_matched_hm, resolutions)

# End Parallelization
message("Ending Parallelization...")
plan("sequential")

# Save the objects
message("Saving the Objects...")
message("MC...")
saveRDS(object = mc_matched_hm, file = "./objects/mc_matched_hm_clustered_2.RDS")
message("Merscope...")
saveRDS(object = vz_matched_hm, file = "./objects/vz_matched_hm_clustered_2.RDS")
message("Xenium...")
saveRDS(object = xen_matched_hm, file = "./objects/xen_matched_hm_clustered_2.RDS")
