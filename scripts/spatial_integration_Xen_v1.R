
## Xenium
# Analysis, Batch-correction & Clustering ----

### make sure all libs are installed a priori!
## load libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(scater)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(png)
  library(cowplot)
  library(parallel)
  library(harmony)
  library(gridExtra)
  library(jpeg)
  library(scales) 
  library(fields)
  library(spatstat)
  library(patchwork)
  library(Matrix)
  library(BiocParallel)    
})

# set R options, paths, object params, names, functions, etc..
message("Setting options, paths and object params..")
Sys.setenv(LANG = "en")
# ..to visualize huge array, standard notation vs scientific one, digits after comma.
options(max.print = 6e+5, scipen = 500, digits = 6)
# set maximum for object size
# https://satijalab.org/seurat/articles/future_vignette.html
options(future.globals.maxSize = 120000 * 1024^2) # ~125G RAM

# function to return metadata of Seurat obj
callmeta <- function(object = NULL) {
  return(object@meta.data) 
}

# Load object ----
dir_path <- "./segmentation_based/objects/"
obj <- readRDS(file.path(dir_path, "vz.re.sc.xen_matched.rds"))
obj

# Update metadata ----
if (grep("SCT_snn", 
         obj %>% 
         callmeta %>% 
         names) %>% any) {
  # remove clusters
  obj@meta.data %<>%
    select(!contains(c("SCT_snn")))
}
obj %>% callmeta %>% str

# Clean object ----
DefaultAssay(obj) <- "MultiTech"
obj %<>% DietSeurat(assays = "MultiTech")

# Subset to keep only Xenium data ----
cells.use <- 
  obj %>%
  callmeta %>%  
  # 1 tech only
  slice(pull(., spatial_tech) %>% grep("Xenium", .)) %>%
  # 2 samples only
  slice(pull(., samples) %>% grep("295|266", .)) %>% rownames
obj$spatial_tech %>% table

# subset obj and FOVs ----
# load modified subset function
source("./scripts/subset_obj_seurat_v2.R")
obj %<>% subset_opt(cells = cells.use)

if (!any(Assays(obj) == "SCT")) {
  start.time <- Sys.time()
  # fast analysis
  message(">>> running `SCTransform`")
  obj %<>% 
    SCTransform(assay = Assays(obj), 
                #clip.range = c(-10, 10), 
                verbose = FALSE)
  end.time <- Sys.time()
  end.time - start.time
} else { 
  message(">>> SCT assay is already present, setting default Assay to SCT")
  DefaultAssay(obj) <- "SCT"
}
obj

# Fast analysis ----
gc() %>% invisible()
start.time <- Sys.time()
obj %<>% 
  RunPCA(npcs = 30, reduction.name = "pca", 
         features = rownames(obj), 
         verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = "pca", reduction.name = "umap.uncorrected", 
          verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
end.time <- Sys.time()
end.time - start.time
obj

if (TRUE) { # set to TRUE to run this chunk
  # Find clusters - Leiden ----
  # run using parallelization
  #plan("sequential") # setback the plan
  plan("multicore", workers = 20)
  plan()
  start.time <- Sys.time()
  obj %<>%
    FindClusters(resolution = c(0.1, 0.2),
                 method = "igraph", 
                 algorithm = 4 # leiden clustering
    ) %>% suppressWarnings()
  end.time <- Sys.time()
  end.time - start.time
  obj
  
  # change cluster vars name
  obj@meta.data %<>%
    select(!seurat_clusters) %>%
    rename(SCT_snn_res.0.1.leiden = SCT_snn_res.0.1,
           SCT_snn_res.0.2.leiden = SCT_snn_res.0.2)
  obj %>% callmeta() %>% str
} 

# Batch-correction using Harmony ----
use.harmony <- TRUE # set to TRUE to run harmony integration
if (use.harmony) {
  DefaultAssay(obj) <- "SCT"
  group.by.vars <- c("samples")
  gc()
  obj <-
    obj %>% harmony::RunHarmony(
      group.by.vars = group.by.vars, 
      assay.use = "SCT",
      theta = rep(2, length(group.by.vars)), 
      lambda = rep(1, length(group.by.vars)),
      reduction = "pca", 
      reduction.save = "harmony", 
      dims.use = 1:30,
      max.iter.harmony = 30, 
      plot_convergence = F)
  obj
  
  # run UMAP on Harmony-corrected PCA ----
  gc()
  obj <- 
    RunUMAP(obj, 
            reduction.name = "umap", 
            #n.epochs = 400, spread = 1, min.dist = 0.3,
            #n.neighbors = 30L,
            umap.method = "uwot", 
            dims = 1:30,
            reduction = "harmony",
            n.components = 2L, 
            verbose = F)
  
  if (TRUE) { # set to TRUE to run this chunk
    # get SNN graph on corrected-PCA
    obj %<>% 
      FindNeighbors(reduction = "harmony", dims = 1:30)
    
    # Find clusters on Harmony-corrected PCA - Leiden ----
    # NOTE: might take 2-3 hours!
    # run using parallelization
    #plan("sequential") # setback the plan
    plan("multicore", workers = 20)
    plan()
    start.time <- Sys.time()
    obj %<>% 
      FindClusters(resolution = c(0.1, 0.2),
                   method = "igraph", 
                   algorithm = 4 # leiden clustering
      ) %>% suppressWarnings()
    end.time <- Sys.time()
    end.time - start.time
    obj
    
    # change cluster vars name
    obj@meta.data %<>%
      select(!seurat_clusters) %>%
      rename(SCT_snn_res.0.1.leiden.hm = SCT_snn_res.0.1,
             SCT_snn_res.0.2.leiden.hm = SCT_snn_res.0.2)
    obj %>% callmeta() %>% str
  }
}

# update..
obj %<>% UpdateSeuratObject()

# Update metadata ----
obj@meta.data %<>%
  mutate(samples = stringi::stri_replace_all_regex(samples, 
                                                   pattern = "region[0-1].|w(.)[^:][0-9].",
                                                   replacement = ""))
# checks on metadata vars
obj %>% 
  callmeta() %>% 
  select(orig.ident, 
         samples, 
         spatial_tech) %>%
  apply(., 2, unique)

# save obj
saveRDS(obj, file = paste0(dir_path, sep = "/", "xen_matched_hm.rds"))
message("Saving - done!")
list.files(dir_path)
