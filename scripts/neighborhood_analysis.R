# Cleanup
rm(list = ls())

library(Seurat)
library(ggplot2)
library(paletteer)
library(data.table)
library(dbscan)
library(pheatmap)
library(tidyr)

# Set working directory
wd <- ""
setwd(wd)
rm(wd)

# Make folders for output
plot_dir <- "neighbourhood_plots"
dir.create(path = plot_dir, showWarnings = F)

data_dir <- "neighbourhood_data"
dir.create(path = data_dir, showWarnings = F)


# Load the Seurat objects with the clustering
message("Loading data")
mc_matched_hm <- readRDS("objects/clustered/re_matched_hm_clustered_2.RDS") # Molecular Cartography
vz_matched_hm <- readRDS("objects/clustered/vz_matched_hm_clustered_2.RDS") # Merscope
xen_matched_hm <- readRDS("objects/clustered/xen_matched_hm_clustered_2.RDS") # Xenium

# Convert coluns with clusters to numeric and set the resolution
mc_matched_hm@meta.data[grep("SCT", colnames(mc_matched_hm@meta.data))] <- lapply(mc_matched_hm@meta.data[grep("SCT", colnames(mc_matched_hm@meta.data))], as.numeric)
vz_matched_hm@meta.data[grep("SCT", colnames(vz_matched_hm@meta.data))] <- lapply(vz_matched_hm@meta.data[grep("SCT", colnames(vz_matched_hm@meta.data))], as.numeric)
xen_matched_hm@meta.data[grep("SCT", colnames(xen_matched_hm@meta.data))] <- lapply(xen_matched_hm@meta.data[grep("SCT", colnames(xen_matched_hm@meta.data))], as.numeric)

mc_resolution <- "SCT_snn_res.0.55"
vz_resolution <- "SCT_snn_res.0.55"
xen_resolution <- "SCT_snn_res.0.55"

# Manually annotated cells types:
# See: https://www.nature.com/articles/s41467-023-44117-x/figures/7
message("Annotating clusters")


# Molecular Cartography 
mc_cell_types <- c(
    "Early CGNP-like", "Differentiated neuronal-like", "Early CGNP-like, proliferating",
     "Stromal / Meningeal",  "Late CGNP-like", "Stromal / Meningeal (2)",
    "Monocytes",  "Migrating CGNP-like" , "Vascular / Endothelial",
    "Migrating CGNP-like (2)", "Other Immune", "TULP1+ cells",
    "Stromal / Meningeal (3)", "Astrocyte / Astrocytic-like")
mc_matched_hm@meta.data$cell_types <- mc_cell_types[mc_matched_hm@meta.data[[mc_resolution]]]
mc_resolution <- "cell_types"

# Merscope
vz_cell_types <- c(
    "Differentiated neuronal-like", "Early CGNP-like", "Early CGNP-like, proliferating",
    "Stromal / Meningeal", "Migrating CGNP-like", "Astrocyte / Astrocytic-like",
    "Differentiated neuronal-like (2)",  "Monocytes", "Differentiated neuronal-like (3)",
    "Vascular / Endothelial", "Differentiated neuronal-like (4)", "Other Immune",
    "Other Immune (2)", "Other Immune (3)", "Early CGNP-like (2)")

vz_matched_hm@meta.data$cell_types <- vz_cell_types[vz_matched_hm@meta.data[[vz_resolution]]]
vz_resolution <- "cell_types"

# Xenium
xen_cell_types <- c(
    "Differentiated neuronal-like", "Differentiated neuronal-like (2)",  "Other Immune",
    "Other Immune (2)", "B Cells", "Oligodendrocyte-like",
    "Differentiated neuronal-like (3)", "Late CGNP-like", "Early CGNP-like",
    "Vascular / Endothelial", "Stromal / Meningeal", "Monocytes",
     "Astrocyte / Astrocytic-like",  "Differentiated neuronal-like (4)", 
    "Early CGNP-like, proliferating")

xen_matched_hm@meta.data$cell_types <- xen_cell_types[xen_matched_hm@meta.data[[xen_resolution]]]
xen_resolution <- "cell_types"

# Simplify the heatmaps by merging / excluding certain cell populations
# Heatmaps with all cell types are also generated
slimmed_cell_types = c("Astrocyte / Astrocytic-like" = "Astrocyte / Astrocytic-like",
    "Differentiated neuronal-like" = "Differentiated neuronal-like",
    "Differentiated neuronal-like (2)" = "Differentiated neuronal-like",
    "Differentiated neuronal-like (3)" = "Differentiated neuronal-like",
    "Differentiated neuronal-like (4)" = "Differentiated neuronal-like",
    "Oligodendrocyte-like" = "NA",
    "Early CGNP-like" = "CGNP-like",
    "Early CGNP-like (2)" = "CGNP-like",
    "Early CGNP-like, proliferating" = "Early CGNP-like, proliferating",
    "Late CGNP-like" = "CGNP-like",
    "Migrating CGNP-like" = "CGNP-like",
    "Migrating CGNP-like (2)" = "CGNP-like",
    "TULP1+ cells" = "NA",
    "Stromal / Meningeal" = "Stromal / Meningeal",
    "Stromal / Meningeal (2)" = "Stromal / Meningeal",
    "Stromal / Meningeal (3)" = "Stromal / Meningeal",
    "Vascular / Endothelial" = "Vascular / Endothelial",
    "Monocytes" = "Monocytes",
    "Other Immune" = "Other Immune",
    "Other Immune (2)" = "Other Immune",
    "Other Immune (3)" = "Other Immune",
    "Other Immune (4)" = "Other Immune",
    "B Cells" = "Other Immune")

mc_matched_hm@meta.data$cell_types_slim <- slimmed_cell_types[mc_matched_hm@meta.data$cell_types]
vz_matched_hm@meta.data$cell_types_slim <- slimmed_cell_types[vz_matched_hm@meta.data$cell_types]
xen_matched_hm@meta.data$cell_types_slim <- slimmed_cell_types[xen_matched_hm@meta.data$cell_types]

# Alphabetically sorted list for plotting
genes <- sort(rownames(xen_matched_hm@assays$MultiTech@meta.features))

# Functions for neighbourhood and permutation tests
##############################
get_neighbours_all <- function(obj, sample_names, fovs, K = 10, n_perm = 5000, cores = 32)
{
   
    cell_list <-lapply(1:length(fovs), function(n){
        # Get the K nearest neighbours
        knn <- kNN(obj@images[[fovs[[n]]]]$centroids@coords, k = K, sort = FALSE)
    
        # Assign cell types
        obj_cell_types <- obj@meta.data[obj@meta.data$orig.ident == sample_names[[n]],]$cell_types
        us_obj_cell_types <- sort(unique(obj_cell_types))
        nbs <- vapply(X = knn$id, function(x){obj_cell_types[[x]]}, character(1))
        dim(nbs) <- dim(knn$id)
        cbind(data.table(cell_type = obj_cell_types), as.data.table(nbs))
    })    
    cells <- rbindlist(cell_list, use.names = TRUE, fill = FALSE, idcol = FALSE)
    
    # Make the interactions simmetrical and count
    cells <- melt(data = cells, id.vars = "cell_type", value.name = "cell_type2", value.factor = FALSE)
    cells[, variable := NULL]
    counts <- copy(cells)
    counts <- as.data.table(complete(data = counts, cell_type, cell_type2))
    counts[, interaction := fifelse(cell_type < cell_type2, paste(cell_type, cell_type2, sep = ";"), paste(cell_type2, cell_type, sep = ";"))]
    counts <- unique(counts[, .N, by = interaction])
    counts[, c("cell_type1", "cell_type2") := tstrsplit(x = interaction, split = ";")]
    counts[, interaction := NULL]
    
    # Make permutations by shuffling the interacting cell types
    permutations <- parallel::mclapply(1:n_perm, function(n){
        pcells <- copy(cells)
        pcells$cell_type <- sample(x = pcells$cell_type, size = length(pcells$cell_type), replace = FALSE)
        pcells$cell_type2 <- sample(x = pcells$cell_type2, size = length(pcells$cell_type2), replace = FALSE)
        pcells <- as.data.table(complete(data = pcells, cell_type, cell_type2))
        pcells[, interaction := fifelse(cell_type < cell_type2, paste(cell_type, cell_type2, sep = ";"), paste(cell_type2, cell_type, sep = ";"))]
        pcells <- unique(pcells[, .N, by = interaction])
        pcells[, c("cell_type1", "cell_type2") := tstrsplit(x = interaction, split = ";")]
        pcells[, interaction := NULL]
        return(pcells)
    }, mc.cores = cores, mc.allow.recursive = TRUE)
    permutations <- rbindlist(permutations, idcol = "permutation")
    
    return(list(cells = cells, counts = counts, permutations = permutations))
    
}

get_pvals_all <- function(permutations, counts, n_perm) {    
    pvals <- copy(permutations)
    pvals[, difference := counts$N - N, by = permutation]
    pvals[, zeta :=  mean(difference) / sd(difference), by = c("cell_type1", "cell_type2")] # µ(a - X) = a - µ(X) and sd(a - X) = sd(X)
    pvals[, greater := (sum(difference < 0) + 1) / (n_perm + 1), by = c("cell_type1", "cell_type2")] # https://pubmed.ncbi.nlm.nih.gov/21044043/
    pvals[, lower := (sum(difference > 0) + 1) / (n_perm + 1), by = c("cell_type1", "cell_type2")]
    pvals[, greater_bh := p.adjust(greater, method = "BH")]
    pvals[, lower_bh := p.adjust(lower, method = "BH")]
    pvals[, pvalue := fifelse(lower < greater, lower, greater)]
    pvals[, pvalue_bh := fifelse(lower_bh < greater_bh, lower_bh, greater_bh)]
    pvals[, score := fifelse(lower < greater, -1 + lower, 1 - greater)]
    pvals[, score := fifelse(pvalue > 0.05, NA, score)]
    pvals[, score_bh := fifelse(lower_bh < greater_bh, -1 + lower_bh, 1 - greater_bh)]
    pvals[, score_bh := fifelse(pvalue_bh > 0.05, NA, score_bh)]
    pvals[, zeta := fifelse(pvalue_bh > 0.05, NA, zeta)]
    pvals <- unique(pvals[, .(cell_type1, cell_type2, zeta, greater, lower, pvalue, score, greater_bh, lower_bh, pvalue_bh, score_bh)])   
    return(pvals)
}

# Heatmaps drawing functions with p-value / z-score
##################################################
nn_heatmapper <- function(pm, title = NULL){
    ggplot(data = pm[order(cell_type1, cell_type2)], aes(x=cell_type1, y=cell_type2, fill=score_bh)) + 
      geom_tile() +  theme_minimal() + coord_fixed() + labs(x = NULL, y = NULL, title = title) +
      scale_fill_gradient2(name = NULL, low = "blue", high = "red", mid = "white",  
       midpoint = 0.0, limit = c(-1,1), space = "Lab") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust = 1))
}

nn_heatmapper_z <- function(pm, title = NULL){
    ggplot(data = pm[order(cell_type1, cell_type2)], aes(x=cell_type1, y=cell_type2, fill=zeta)) + 
      geom_tile() +  theme_minimal() + coord_fixed() + labs(x = NULL, y = NULL, title = title) +
      scale_fill_gradient2(name = NULL, low = "blue", high = "red", mid = "white",  
       midpoint = 0.0, space = "Lab") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust = 1))
}
####################################################

###################################################################
n_perm = 10000 # Permutations
K = 10         # Number of nearest neighbours  

message("Subsetting objects")
xen_nn_obj <- subset(xen_matched_hm, subset = cell_types_slim != "NA")
mc_nn_obj <- subset(mc_matched_hm, subset = cell_types_slim != "NA")
vz_nn_obj <- subset(vz_matched_hm, subset = cell_types_slim != "NA")

xen_nn_obj@meta.data$cell_types <- slimmed_cell_types[xen_nn_obj@meta.data$cell_types]
mc_nn_obj@meta.data$cell_types <- slimmed_cell_types[mc_nn_obj@meta.data$cell_types]
vz_nn_obj@meta.data$cell_types <- slimmed_cell_types[vz_nn_obj@meta.data$cell_types]

message("Neighbours Xenium all clusters")
xen_neighbours_all <- get_neighbours_all(xen_matched_hm, c("mb266.Xenium", "mb295.Xenium"), c("mb266.Xenium", "mb295.Xenium"), K = K, n_perm = n_perm)
message("Pvalues Xenium all clusters")
xen_pvals_all <- get_pvals_all(xen_neighbours_all$permutations, xen_neighbours_all$counts, n_perm = n_perm)
message("Saving data")
fwrite(x = xen_pvals_all, file = file.path(data_dir, "/xen_nn_10-10000_all_pvals.csv"))
message("Making heatmaps")
xen_nn_heatmap_all_p <- nn_heatmapper(xen_pvals_all, "Xenium")
xen_nn_heatmap_all_z <- nn_heatmapper_z(xen_pvals_all, "Xenium")
message("Saving heatmaps")
ggsave(plot = xen_nn_heatmap_all_p, filename = file.path(plot_dir, "/xen_nn_10-10000_all_p_bh_heatmap.pdf"), width = 20, height = 20)
ggsave(plot = xen_nn_heatmap_all_z, filename = file.path(plot_dir, "/xen_nn_10-10000_all_z_bh_heatmap.pdf"), width = 20, height = 20)
rm(xen_neighbours_all)
gc()

message("Neighbours Xenium selected clusters")
xen_neighbours_slim <- get_neighbours_all(xen_nn_obj, c("mb266.Xenium", "mb295.Xenium"), c("mb266.Xenium", "mb295.Xenium"), K = K, n_perm = n_perm)
message("Pvalues Xenium selected clusters")
xen_pvals_slim <- get_pvals_all(xen_neighbours_slim$permutations, xen_neighbours_slim$counts, n_perm = n_perm)
message("Saving data")
fwrite(x = xen_pvals_slim, file = file.path(data_dir, "/xen_nn_10-10000_slim_pvals.csv"))
message("Making heatmaps")
xen_nn_heatmap_slim_p <- nn_heatmapper(xen_pvals_slim, "Xenium")
xen_nn_heatmap_slim_z <- nn_heatmapper_z(xen_pvals_slim, "Xenium")
message("Saving heatmaps")
ggsave(plot = xen_nn_heatmap_slim_p, filename = file.path(plot_dir, "/xen_nn_10-10000_slim_p_bh_heatmap.pdf"), width = 20, height = 20)
ggsave(plot = xen_nn_heatmap_slim_z, filename = file.path(plot_dir, "/xen_nn_10-10000_slim_z_bh_heatmap.pdf"), width = 20, height = 20)
rm(xen_neighbours_slim)
gc()

message("Neighbours MC all clusters")
mc_neighbours_all <- get_neighbours_all(mc_matched_hm, c("mb266.w5a1.MC", "mb295.w6a1.MC", "mb295.w6a2.MC"), c("mb266.w5a1.Resolve", "mb295.w6a1.Resolve", "mb295.w6a2.Resolve"), K = K, n_perm = n_perm)
message("Pvalues MC all clusters")
mc_pvals_all <- get_pvals_all(mc_neighbours_all$permutations, mc_neighbours_all$counts, n_perm = n_perm)
message("Saving data")
fwrite(x = mc_pvals_all, file = file.path(data_dir, "/mc_nn_10-10000_all_pvals.csv"))
message("Making heatmaps")
mc_nn_heatmap_all_p <- nn_heatmapper(mc_pvals_all, "Molecular Cartography")
mc_nn_heatmap_all_z <- nn_heatmapper_z(mc_pvals_all, "Molecular Cartography")
message("Saving heatmaps")
ggsave(plot = mc_nn_heatmap_all_p, filename = file.path(plot_dir, "/mc_nn_10-10000_all_p_bh_heatmap.pdf"), width = 20, height = 20)
ggsave(plot = mc_nn_heatmap_all_z, filename = file.path(plot_dir, "/mc_nn_10-10000_all_z_bh_heatmap.pdf"), width = 20, height = 20)
rm(mc_neighbours_all)
gc()

message("Neighbours MC selected clusters")
mc_neighbours_slim <- get_neighbours_all(mc_nn_obj, c("mb266.w5a1.MC", "mb295.w6a1.MC", "mb295.w6a2.MC"), c("mb266.w5a1.Resolve", "mb295.w6a1.Resolve", "mb295.w6a2.Resolve"), K =K, n_perm = n_perm)
message("Pvalues MC selected clusters")
mc_pvals_slim <- get_pvals_all(mc_neighbours_slim$permutations, mc_neighbours_slim$counts, n_perm = n_perm)
message("Saving data")
fwrite(x = mc_pvals_slim, file = file.path(data_dir, "/mc_nn_10-10000_slim_pvals.csv"))
message("Making heatmaps")
mc_nn_heatmap_slim_p <- nn_heatmapper(mc_pvals_slim, "Molecular Cartography")
mc_nn_heatmap_slim_z <- nn_heatmapper_z(mc_pvals_slim, "Molecular Cartography")
message("Saving heatmaps")
ggsave(plot = mc_nn_heatmap_slim_p, filename = file.path(plot_dir, "/mc_nn_10-10000_slim_p_bh_heatmap.pdf"), width = 20, height = 20)
ggsave(plot = mc_nn_heatmap_slim_z, filename = file.path(plot_dir, "/mc_nn_10-10000_slim_z_bh_heatmap.pdf"), width = 20, height = 20)
rm(mc_neighbours_slim)
gc()

message("Neighbours Merscope all clusters")
vz_neighbours_all <- get_neighbours_all(vz_matched_hm, c("mb266.region1.MERFISH", "mb295.region0.MERFISH"), c("mb266.region1.Vizgen", "mb295.region0.Vizgen"), K = K, n_perm = n_perm)
message("Pvalues Merscope all clusters")
vz_pvals_all <- get_pvals_all(vz_neighbours_all$permutations, vz_neighbours_all$counts, n_perm = n_perm)
message("Saving data")
fwrite(x = vz_pvals_all, file = file.path(data_dir, "/vz_nn_100-10000_all_pvals.csv"))
message("Making heatmaps")
vz_nn_heatmap_all_p <- nn_heatmapper(vz_pvals_all, "Merscope")
vz_nn_heatmap_all_z <- nn_heatmapper_z(vz_pvals_all, "Merscope")
message("Saving heatmaps")
ggsave(plot = vz_nn_heatmap_all_p, filename = file.path(plot_dir, "/vz_nn_10-10000_all_p_bh_heatmap.pdf"), width = 20, height = 20)
ggsave(plot = vz_nn_heatmap_all_z, filename = file.path(plot_dir, "/vz_nn_10-10000_all_z_bh_heatmap.pdf"), width = 20, height = 20)
rm(vz_neighbours_all)
gc()

message("Neighbours Merscope selected clusters")
vz_neighbours_slim <- get_neighbours_all(vz_nn_obj, c("mb266.region1.MERFISH", "mb295.region0.MERFISH"), c("mb266.region1.Vizgen", "mb295.region0.Vizgen"), K = K, n_perm = n_perm)
message("Pvalues Merscope selected clusters")
vz_pvals_slim <- get_pvals_all(vz_neighbours_slim$permutations, vz_neighbours_slim$counts, n_perm = n_perm)
message("Saving data")
fwrite(x = vz_pvals_slim, file = file.path(data_dir, "/vz_nn_10-10000_slim_pvals.csv"))
message("Making heatmaps")
vz_nn_heatmap_slim_p <- nn_heatmapper(vz_pvals_slim, "Merscope")
vz_nn_heatmap_slim_z <- nn_heatmapper_z(vz_pvals_slim, "Merscope")
message("Saving heatmaps")
ggsave(plot = vz_nn_heatmap_slim_p, filename = file.path(plot_dir, "/vz_nn_10-10000_slim_p_bh_heatmap.pdf"), width = 20, height = 20)
ggsave(plot = vz_nn_heatmap_slim_z, filename = file.path(plot_dir, "/vz_nn_10-10000_slim_z_bh_heatmap.pdf"), width = 20, height = 20)
rm(vz_neighbours_slim)
gc()

