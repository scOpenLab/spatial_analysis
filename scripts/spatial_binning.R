#' @import dplyr
#' @import tibble
#' @import purrr
#' @import furrr
#' @import tidyr
#' @import magrittr
#' @importFrom Matrix import t
NULL

#' Read and Load..
#' NOTE: input data must have these variables -> \code{x, y, gene}


#' Function to do spatial binning with a custom bin.size
#' params/args:
#' @param file.path Provide path \code{data.frame} object containg \code{$x, $y} coordinates and 
#'   \code{$molecule}. Eg, extraction from Seurat object for 1st FOV 
#'   \code{obj[[Images(obj)[1]]][["molecules"]] %>% GetTissueCoordinates }
#' @param input.df Alternatively, if \code{file.path} is empty, provide \code{data.frame} object
#' @param bin.size Size (in µm or px) of the square bin, eg \code{500} means -> \code{500 x 500} bins size.
#' @param coord.space Define in which space coordinates are located, Default to \code{"microns"}
#' @param spatial.technology Provide spatinal technology name from where the molecule coordinates orignate; 
#' Default is to \code{NaN}, this will use simple \code{data.frame} with \code{x, y, gene} coords.
#' choose only one of:
#' \itemize{
#'  \item \dQuote{Vizgen}: ..from Vizgen MERSCOPE
#'  \item \dQuote{Xenium}: ..from 10X Xenium
#'  \item \dQuote{Resolve}: ..from Resolve Bioscience
#'  \item \dQuote{RNAScope}: ..from RNAScope 
#' }
#' @param MulticoreParam.cores Number of cores to use for \code{BiocParallel::bplapply()}, default is to use \code{50%} of total available cores.
#' @param ..

## ==============================================================
## Reading molecule coordinates data
## ==============================================================

BinMols <- function(
    file.path = NULL,
    input.df = NULL, #supports df input
    bin.size = 500, 
    coord.space = "micron",
    to.microns.RNAScope = FALSE, # use only if RNAScope data is in pixel space, and want to convert to micron space.
    to.microns.Resolve = FALSE, # use only if Resolve data is in pixel space, and want to convert to micron space.
    spatial.technology = NaN,
    default.data = FALSE, # set to TRUE only if input data contains only x,y coords and at least 1 character variable
    checks.requireNamespace = TRUE,
    show.progress = TRUE, # If to show progress bar
    use.BiocParallel = TRUE, # if to use parallelization for fast computation, else using apply family
    MulticoreParam.cores = quantile(BiocParallel::multicoreWorkers() %>% seq())[3] %>% ceiling()
)
{
  # packages that needs to be installed a priori
  if (checks.requireNamespace) {
    pkgs <- c("data.table", "arrow", "sfarrow", 
              "stringr", "plyr", "pbapply", "future",
              "tidyverse", "dtplyr", "furrr", "BiocParallel")
    lapply(pkgs %>% length %>% seq, function(i)
    { !requireNamespace(pkgs[i], quietly = TRUE) } ) %>% 
      unlist %>% 
      { if (c(which(.) > 0) %>% any()) 
        { c("Please install ->", "\n",
            paste0("'", pkgs[which(.)], "'", collapse = ", "), " for this function") %>% 
          stop(., call. = FALSE) } }
    }
  
  # Read molecule coordinated (coords) ----	
  # check if file path or dataframe object is provides	
  if (is.null(input.df) && is.null(file.path)) { 
    stop("Provide (molecule) coordinates in `file.path` or in `input.df`")
  } else if (!is.null(input.df) && !is.null(file.path)) {
    stop("Provide (molecule) coordinates either in `file.path` or in `input.df`")
  }
  
  ## Sanity checks on data.frame or file.path provided ----  	
  if (!is.null(input.df) && is.null(file.path) &&
      any(class(input.df) == c("data.table", "data.frame"))) {
    # use dataframe as input
    mols.coord <- input.df %>% 
      # convert for faster dplyr operation
      dtplyr::lazy_dt(.)
  } else if (spatial.technology == "Xenium" && !length(grep("parquet", file.path)) == 0) {
    # read coords data - `.parquet`
    mols.coord <- 
      arrow::read_parquet(file.path, as_data_frame = FALSE) %>%
        mutate(feature_name = cast(feature_name, arrow::string())) %>%
        as_tibble() %>%
        # convert for faster dplyr operation
        dtplyr::lazy_dt(.)
  } else if (spatial.technology == "Resolve" && 
             !length(grep("Panorama|run", file.path)) == 0) {
    # read Resolve-specific coords data
    mols.coord <- 
      data.table::fread(input = file.path, blank.lines.skip = TRUE,
                        #sep = ",",
                        data.table = TRUE) %>%
      # add col names
      # NOTE: x,y coords needs to be flipped to matche the original image
      data.table::setnames(., c("y", "x", "z", "gene")) %>%
      # convert for faster dplyr operation
      dtplyr::lazy_dt(.)
  } else if (spatial.technology == "RNAScope") {
    # read RNAScope-specific coords data
    mols.coord <- 
      data.table::fread(input = file.path, blank.lines.skip = TRUE,
                        #sep = ",",
                        data.table = TRUE) %>%
      # get only first 4 cols
      select(1:4) %>%
      # add col names
      # NOTE: x,y coords needs to be flipped to matche the original image
      data.table::setnames(., c("x", "y", "z", "gene")) %>%
      # convert for faster dplyr operation
      dtplyr::lazy_dt(.)
  } else {
    # read Resolve-specific coords data
    mols.coord <- 
      data.table::fread(input = file.path, blank.lines.skip = TRUE,
                        #sep = ",",
                        data.table = TRUE) %>%
      # convert for faster dplyr operation
      dtplyr::lazy_dt(.)
  }
  
  ## Sanity checks on coordinates data ----
  # guess what type of spatial technology data it is
  #..based on variable names (currently only Vizgen and Xenium)
  if(is.nan(spatial.technology)) {
    if (grep("global_", 
             mols.coord %>% names) %>% length() == 3) {
      # set technology name
      message("Spatial Technology is likely to be Vizgen/Merfish")
      spatial.technology <- "Vizgen"
      } else if (grep("cell_id|overlaps_nucleus|feature_name|_location|molecule",
                    mols.coord %>% names) %>% length == 6) {
      # set technology name
      message("Spatial Technology is likely to be 10X Xenium")
      spatial.technology <- "Xenium"
      }
    }
  
  if (!default.data && is.nan(spatial.technology)) {
    mols.coord %>%
      select(., contains(c("x", "y", "gene"))) %>% 
      # sanity on if 2 vars are numeric/integers and 1 var is char
      select_if(., ~ is.character(.) | is.numeric(.) | is.integer(.)) %>% 
      { if (!ncol(.) == 3)
        { stop("Provided coordinate data does NOT contain `x, y, gene` variable (cols)") }} 
    } else if (default.data && is.nan(spatial.technology)) {
      mols.coord %>%
        # bind together x,y coords and a character variable
        bind_cols(.,
                  select(., contains(c("x", "y"))) %>%
                  select_if(., ~ is.character(.))
                  ) %>%
        { if (!ncol(.) == 3) { 
          stop("Provided coordinate data must contain `x, y` & character variable (cols)")}}
      }
  
  ## Filter - for technology-specific or default data ----
  # keep all z-stacks, select only x,y coords and gene col
  mols.coord %<>%
    { 
      # Vizgen - in micron space
      if (spatial.technology == "Vizgen" && coord.space == "micron")
        { if (length(grep("detected_transcripts.csv", file.path)) == 0) {
          warning("Vizgen:", "\n", 
                  ">>> Standard pattern 'detected_transcripts.csv' is NOT detected in `file.path`", 
                  immediate. = TRUE)
          } else {
            message("Reading Data:", "\n", ">>> Using Vizgen data (in ", coord.space, " space)",
                    if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                      } else { paste("") })
            select(., starts_with(c("global", "gene"))) %>%
              # rename all coords to "x", "y"
              rename(x = global_x, y = global_y) }
        
      # Vizgen - in pixel space
      } else if (spatial.technology == "Vizgen" && coord.space == "pixel")
        { if (length(grep("detected_transcripts.csv", file.path)) == 0) {
          warning("Vizgen:", "\n", 
                  ">>> Standard pattern 'detected_transcripts.csv' is NOT detected in `file.path`", 
                  immediate. = TRUE)
          } else {
            message("Reading Data:", "\n", ">>> Using Vizgen data (in ", coord.space, " space)", 
                    if (!is.null(file.path)) { paste(" from: ", "\n", file.path)
                      } else { paste("") })
            select(., starts_with(c("x", "y", "gene"))) }
        
      # Resolve - in micron space
      } else if (spatial.technology == "Resolve" && coord.space == "micron")
        { if (length(grep("Panorama|run", file.path)) == 0) {
          warning("Resolve:", "\n", 
                  ">>> Standard pattern 'Panorama|run' is NOT detected in `file.path`", 
                  immediate. = TRUE)
        } else {
          message("Reading Data:", "\n", ">>> Using Resolve data (in ", coord.space, " space)", 
                  if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                    } else { paste("") })
            select(., starts_with(c("x", "y", "gene"))) }
        
      # Resolve - convert from pixel to micron space 
      } else if (spatial.technology == "Resolve" && coord.space == "pixel" && 
                 to.microns.Resolve)
      { if (length(grep("Panorama|run", file.path)) == 0) {
        warning("Resolve:", "\n", 
                ">>> Standard pattern 'Panorama|run' is NOT detected in `file.path`", 
                immediate. = TRUE)
      } else {
        message("Reading Data:", "\n", ">>> Using Resolve data (in ", coord.space, " space)", 
                if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                  } else { paste("") }, "\n", 
                ">>> Converting Resolve data from pixel to micron space")
        select(., starts_with(c("x", "y", "gene"))) %>%
          # convert coords to microns -> pixel size is 0.138 µm
          mutate(x = x * 0.138, y = y * 0.138) }
      
      # Resolve - in pixel space
      } else if (spatial.technology == "Resolve" && coord.space == "pixel")
      { if (length(grep("Panorama|run", file.path)) == 0) {
        warning("Resolve:", "\n", 
                ">>> Standard pattern 'Panorama|run' is NOT detected in `file.path`", 
                immediate. = TRUE)
      } else {
        message("Reading Data:", "\n", ">>> Using Resolve data (in ", coord.space, " space)", 
                if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                  } else { paste("") })
        select(., starts_with(c("x", "y", "gene"))) }  
        
      # Xenium - in micron space only
      } else if (spatial.technology == "Xenium" && coord.space == "micron")
      { if (length(grep("transcripts.", file.path)) == 0 && !is.null(input.df)) {
        warning("Xenium:", "\n",
                ">>> Standard pattern 'transcripts.' is NOT detected in `file.path`", 
                immediate. = TRUE)
        { if (tibble::as_tibble(.) %>% colnames() %>% 
              grep("molecule", x = .) %>% any) {
          select(., contains(c("molecule", "x", "y"))) %>%
            # rename all coords to "x", "y"
            rename(., gene = molecule)}
          }
        } else {
        message("Reading Data:", "\n", ">>> Using Xenium data (in micron space)", 
                if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                  } else { paste("") })
        select(., contains(c("feature", "location"))) %>%
          # rename all coords to "x", "y"
          rename(., x = x_location, y = y_location, gene = feature_name) }
        
      # RNAScope - convert from pixel to micron space 
      } else if (spatial.technology == "RNAScope" && 
                 coord.space == "pixel" &&
                 to.microns.RNAScope) {
      { message("Reading Data:", "\n", ">>> Using RNAScope data (in ", coord.space, " space)", 
                if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                } else { paste("") }, "\n", 
                ">>> Converting RNAScope data from pixel to micron space")
        select(., starts_with(c("x", "y", "gene"))) %>%
          # convert coords to microns -> pixel size is 0.3014844 µm
          mutate(x = x * 0.3014844, y = y * 0.3014844) }
      
      # RNAScope - in pixel space
      } else if (spatial.technology == "RNAScope" && coord.space == "pixel") {
      { message("Reading Data:", "\n", ">>> Using RNAScope data (in ", coord.space, " space)", 
                if (!is.null(file.path)) { paste(" from: ", "\n", file.path) 
                } else { paste("") })
        select(., starts_with(c("x", "y", "gene"))) } 
        
      # Default data - from some spatial technology, eg Visium etc..
      } else if (is.nan(spatial.technology) && !default.data) 
        { message("Reading Data:", "\n", ">>> Using Default data (in ", coord.space, " space)", 
                  if (!is.null(file.path)) { paste(" from: ", "\n", file.path)
                    } else { paste("") })
        select(., contains(c("x", "y", "gene"))) %>% 
        # sanity: if 2 vars are numeric/integers and 1 var is char
        select_if(., ~ is.character(.) | is.numeric(.) | is.integer(.)) %>% 
          { if (!ncol(.) == 3) { 
          stop("Provided coordinate data must contain `x, y, gene` variables!")
          } else { (.) } }
        
      # Default data - any technology  
      } else if (is.nan(spatial.technology) && default.data) 
        { message("Reading Data:", "\n", ">>> Using Default data (in ", coord.space, " space)",
                  if (!is.null(file.path)) { paste(" from: ", "\n", file.path)
                    } else { paste("") })
        # bind together x,y coords and a character variable
        bind_cols(.,
                  select(., contains(c("x", "y"))) %>%
                  select_if(., ~ is.character(.))
                  ) %>%
          { if (!ncol(.) == 3) { 
          stop("Provided coordinate data must contain `x, y` & character variable (cols)") 
            } else { (.) } }
        }
    } %>%
      select(!contains("z")) %>% # keep only xy coords
      # add unique ID for each mol coord..
      #..to trace later for duplicates
    
      # use fast adding of bin ids ----
      mutate(coord_ID = stringi::stri_c(nrow(.) %>% seq,
                                        "_coord")) %>%
      # convert to data.table format
      data.table::as.data.table()
      message(">>> Added unique ID to coordinates")
      
  ## Custom helper funtion for parallel or standard processing ----
  custom.lapply <- function(...,
                         show.progressbar = show.progress, 
                         use.parallelization = use.BiocParallel) {
    if (use.parallelization) { 
      # using parallelization
      BiocParallel::bplapply(...,
                             BPPARAM = BiocParallel::MulticoreParam(MulticoreParam.cores,
                                                                    tasks = 20L,  
                                                                    force.GC = FALSE, 
                                                                    progressbar = show.progressbar)) 
    } else {
      if (show.progressbar) {
        #..will show progress bar 
        pbapply::pboptions(type = "timer", style = 1, char = "=")
        pbapply::pblapply(...) # use apply family..
      } else {
        lapply(...) # use apply family..
      }
    }
    }
  
  ## Sanity check on molecule coords - check for duplicates ----
  message("Sanity check:", "\n",
          ">>> Checking data coordinates for duplicates - per gene or (char) variable")		
  # get gene or char var list	
  markers <- 
    mols.coord %>%
    { if (default.data && is.nan(spatial.technology)) { 
      select_if(., is.character) %>% 
        # extract char variable
        pull(1) 
      } else { pull(., gene) }} %>%
    unique()
  
  # check if any mols coords are duplicated - per gene
  test.dupl <-
    custom.lapply(markers %>% length %>% seq, function(i) {
      mols.coord %>%
        dtplyr::lazy_dt(.) %>%
        filter(gene == markers[i]) %>%
        select(!contains("coord_")) %>%
        as_tibble() %>%
        duplicated.data.frame() %>% any()
    })
  #test.dupl %>% unlist %>% any # should return FALSE if no duplicated
  if(!any(test.dupl %>% unlist)) {
    message(">>> Coordinate data is clean - no duplicates!")
  } else {	
    warning("Some coordinates are duplicated", "\n", 
            ">>> Consider checking for duplicated coords - per gene or (char) variable", 
            immediate. = TRUE)
  }
  
  ## Generate Spatial bins of a custom bin size ----
  ### outline of the workflow:
  # generate min-max bins ranges, eg 500x500 microns
  # ie area of 250000 µm^2 or area of 0.5x0.5 = 0.25 mm^2
  # adding/subtracting c(bin.size/2) to/from max/min
  
  message(">>> Generating Spatial Bins according to provided coords..")	
  # generate min-max bins ranges
  # function to make bins
  make.bins <- function(descartes, # a coordiate
                        bin_size = NULL # spatial bin size
                        ) 
    { c(min(descartes) - c(bin_size / 2), # subtract half of bin_size from min of coord
        seq(min(descartes) - c(bin_size / 2) + bin_size,
            max(descartes) + c(bin_size / 2) - bin_size,
            by = bin_size),
        max(descartes) + c(bin_size / 2) # add half of bin_size to max of coord
    )
    }
  coords_bins <-
    lapply(grep("x|y", names(mols.coord)), function(i) 
      { mols.coord[[i]] %>% make.bins(., bin_size = bin.size) }
      )
  
  # find which vector of coords has the smallest length
  min.coord.index <- 
    lapply(coords_bins %>% seq, function(i) 
    { coords_bins[[i]] %>% length }) %>% 
    unlist %>% which.min
  
  # calculate difference in lenght of coords
  coord.length.diff <-
    lapply(coords_bins %>% seq, function(i) 
    { coords_bins[[i]] %>% length }) %>% 
    unlist %>% sort
  coord.length.diff <- coord.length.diff[2] - coord.length.diff[1]
  
  # combine: added bin.size to each last value
  coords_bins[[min.coord.index]] <-
    c(coords_bins[[min.coord.index]],
      seq(tail(coords_bins[[min.coord.index]], n = 1),
          c(tail(coords_bins[[min.coord.index]], n = 1) + bin.size * coord.length.diff),
          length.out = coord.length.diff + 1)) %>% unique
  names(coords_bins) <- c("x", "y")
  # convert to matrix or dataframe
  coords_bins %<>% as.data.frame.list	
  
  # initial 1st spatial bin
  i <- 0
  init_bin <-
    bind_cols(rep(coords_bins[1:2 + i, 1], 2) %>% sort,
              c(coords_bins[1:2 + i, 2], coords_bins[1:2 + i, 2] %>% rev)) %>% 
    #as.data.frame.matrix %>%
    mutate(id = i) %>% 
    data.table::setnames(new = c("x", "y", "id")) %>%
    suppressMessages()
  
  # TODO: (optionally) optimize bin construction with either of: ----
    # - data.table::CJ
    # - tidyr::expand_grid
  
  # generate - spatial bins/areas along y, ie starting from lower left corner
  coords_bins_along.y <- 
    lapply(c(0, seq(nrow(coords_bins) - 1)),
           function(i) { 
             bind_cols(init_bin[,1], # keep x as it is
                       init_bin[,2] + bin.size * i, # add bin.size µm to y coord
                       init_bin[,3] + i + 1 # add id
             ) %>% 
               data.table::setnames(new = c("x", "y", "id")) %>%
               suppressMessages()
           })
  
  # ..move 1 bin along x, then along y again, adding each time a given bin.size
  bin.size.add <- seq(bin.size, 
                      length.out = length(coords_bins_along.y), 
                      by = bin.size)
  coords_bins_along.x <- 
    custom.lapply(coords_bins_along.y %>% seq(),
               function(i) { 
                 bind_cols(coords_bins_along.y[[1]][,1] + bin.size.add[i], # add bin.size µm to x coord
                           coords_bins_along.y[[1]][,2], # keep y as it is
                           length(coords_bins_along.y) + i # add id
                 ) %>%
                   data.table::setnames(new = c("x", "y", "id")) %>%
                   suppressMessages()
               }, show.progressbar = show.progress
               )
  
  # generate - full spatial bins/areas
  # make steps for unique id, ie each id corresponds to single bin
  gc() %>% invisible()
  steps_id <- 
    seq(unique(coords_bins_along.x[[length(coords_bins_along.x)]]$id) + 1, 
        length.out = length(coords_bins_along.x),
        by = length(coords_bins_along.x))
  # move each of coords_bins_along.x along y
  coords_bins_along.x_moved <-
    custom.lapply(length(coords_bins_along.x) %>% seq(), function(j) {
      # substract 1 to avoid duplicated bins
      custom.lapply(length(coords_bins_along.x) %>% seq() %>% .[-1],
                 function(i) {
                   bind_cols(coords_bins_along.x[[j]] %>% dtplyr::lazy_dt(.) %>% pull(x), # keep x the same
                             # add bin.size to y coord, except 1st one
                             coords_bins_along.x[[j]] %>% dtplyr::lazy_dt(.) %>% pull(y) + bin.size.add[i] - bin.size, 
                             seq(steps_id[j], 
                                 length.out = length(coords_bins_along.x))[i]
                   ) %>% 
                     data.table::setnames(., new = c("x", "y", "id")) %>%
                     suppressMessages()
                 }, 
                 show.progressbar = FALSE,
                 use.parallelization = use.BiocParallel
              )
      },
      show.progressbar = show.progress,
      use.parallelization = FALSE
      )
  
  # combine all binned coords/dfs
  coords_binned <- 
    lapply(coords_bins_along.x_moved %>% seq(), 
           function(i) { coords_bins_along.x_moved[[i]] %>% data.table::rbindlist(.) }) %>% 
    c(coords_bins_along.y, coords_bins_along.x, .) %>% data.table::rbindlist(.)
  
  # rename bins id
  coords_binned %<>%
    # remove NAs if present
    tidyr::drop_na(.) %>%
    mutate(id = lapply(id %>% unique() %>% seq(), 
                       function(i) rep_len(i, 4)) %>% unlist()) %>%
    # convert to tibble format
    as_tibble()

  ## Filter molecule coords using Spatial bins, ie ..remove empty bins ----
  # sanity check if the mols coords are of class "data.table"
  if (!any(class(mols.coord) == c("dtplyr_step_first", "dtplyr_step"))) {
    message(">>> Converting `data.table` or `data.frame` to `dtplyr`-like class")
    mols.coord %<>% dtplyr::lazy_dt(.)
  }
  
  message(">>> Spatial Binning of coordinates..")
  mols.coord.df <- 
    custom.lapply(coords_binned %>% pull(id) %>% unique(), function(i) {
      # get bin range
      x_range <- filter(coords_binned %>% dtplyr::lazy_dt(.), id == i) %>% 
        pull(x) %>% unique()
      y_range <- filter(coords_binned %>% dtplyr::lazy_dt(.), id == i) %>% 
        pull(y) %>% unique()
      # df of mols within the bin range
      df_temp <- 
        filter(mols.coord,
               between(x, x_range[1], x_range[2]) & 
                 between(y, y_range[1], y_range[2])
               )
      
      # check if df is emply
      if (!plyr::empty(df_temp %>% as_tibble())) {
        # pass df further only if not empty
        df_temp %>%
          # add area/spot/bin name, to be treated as a "cell"
          mutate(spatial_bin = paste0("spatial_bin_", i)) %>%
          data.table::as.data.table(.)
      }
      }, 
      show.progressbar = show.progress,
      use.parallelization = use.BiocParallel
    )
  
  # clean up list - remove empty dfs from the list
  #mols.coord.df %<>% purrr::discard(., ~ is.null(.))

  # convert spatially binned molecule coords to single dataframe ----
  mols.coord.df %<>%
    data.table::rbindlist(.) %>% 
    # remove duplicates
    distinct(coord_ID, .keep_all = TRUE)
  
  # Sanity check on binned molecule coordinates for duplicates - per gene ----
  message("Sanity check:", "\n",
          ">>> Checking Binned coordinates for duplicates - per gene or (char) variable")		
  # get gene or char var list
  markers <- 
    mols.coord %>%
    { if (default.data && is.nan(spatial.technology)) { 
      select_if(., is.character) %>% 
        # extract char variable
        pull(1) 
    } else { pull(., gene) }} %>%
    unique()
  
  # check if any mols coords are duplicated - per gene
  test.dupl <- 
    custom.lapply(markers %>% length %>% seq, function(i) {
      mols.coord.df %>% 
        dtplyr::lazy_dt(.) %>%
        filter(gene == markers[i]) %>%
        select(!contains(c("spatial_", "coord_"))) %>%
        as_tibble() %>%
        duplicated.data.frame() %>% any()
        })
  #test.dupl %>% unlist %>% any # should return FALSE if no duplicated
  if(!any(test.dupl %>% unlist)) {
    message(">>> Binned coordinates data is clean - no duplicates!")
  } else {	
    warning("Some Binned coordinates are duplicated", "\n", 
            ">>> Consider checking for duplicated coords - per gene or (char) variable", 
            immediate. = TRUE)
  }
  
  outs <- list(bin.size = bin.size, # size of the bin
               bins.min.max = coords_bins, # min & max of bin coords
               bins.boxes = coords_binned, # box of bin coords
               mols.coord.df = mols.coord.df, # spatially binned (molecule) coords as single dataframe
               custom.lapply = custom.lapply # store a custom function as well
  )
  return(outs)
}


## ==============================================================
## Loading molecule coordinates data - output is:
### - merged Seurat object with spatial bins (as cells or spots)
## ==============================================================

BinMols.Seurat <- function(
    file.path = NULL,
    fov = NULL,
    assay = NULL,
    spatial.technology = NaN,
    use.furrr = FALSE,
    update.object = TRUE,
    MulticoreParam.cores = 12,
    ...)
  
{
  data <- BinMols(file.path = file.path, 
                  spatial.technology = spatial.technology,
                  ...)
  
  # set spatially binned (molecule) coords df
  mols.coord.df <- data[["mols.coord.df"]]
  
  # extract funtion 
  custom.lapply <- data$custom.lapply
  
  # set fov and assay names if not provided
  if (is.null(fov) && is.null(assay)) {	
    if (spatial.technology == "Vizgen") {
      fov <- "vz"; assay <- "Vizgen"  
    } else if (spatial.technology == "Resolve") { fov <- "re"; assay <- "Resolve" 
    
    # for RNAScope
    } else if (spatial.technology == "RNAScope") { fov <- "rs"; assay <- "RNAScope"   
    } else if (is.null(spatial.technology)) { fov <- "st"; assay <- "Spatial" }	
  }
  
  # get original planned session if exists
  if (is(future::plan(), "multisession")) {
    orig.plan <- future::plan()
  }
  
  message(">>> Generating count matrix from Binned coords..")
  # make a single df for count matrix

  mols.counts <- 
    mols.coord.df %>%
    # fast dplyr operations
    dtplyr::lazy_dt(.) %>% 
    group_by(gene, spatial_bin) %>% 
    summarise(nCounts = n()) %>% # make counts per gene & bin
    suppressMessages() %>%
    ungroup() %>% 
    # convert back to data.table
    data.table::as.data.table(.) %>%
    group_by(spatial_bin) %>%
    group_split() %>% # split df into list of dfs per bin
    ## apply fast `purrr` functions
    { 
      if (use.furrr) {
        # set temporary workers
        f.plan <- future::plan("multisession", workers = 5L, gc = TRUE)
        #on.exit(f.plan %>% future::plan()) # exiting doesn't to work with furrr
        # using furrr (faster purrr with future)
        # add spatial bin name to counts col
        furrr::future_imap(., function(x, y) x %>% 
                             rename_with(stringr::str_replace, 
                                         pattern = "nCounts",
                                         replacement = pull(., spatial_bin) %>% 
                                           unique)) %>%
          # remove spatial bin col
          furrr::future_map(., ~ select(., !spatial_bin))
        } else {
          purrr::imap(., function(x, y) x %>% 
                        rename_with(stringr::str_replace,
                                    pattern = "nCounts",
                                    replacement = pull(., spatial_bin) %>% 
                                      unique)) %>%
            # remove spatial bin col
            purrr::map(., ~ select(., !spatial_bin))
          }
      } %>% 
    suppressPackageStartupMessages() %>%
    # merge all dfs into single df
    purrr::reduce(., full_join, by = "gene") %>%
    # add rownames
    tibble::column_to_rownames(., var = "gene") %>%
    # replace NA values with 0
    replace(is.na(.), 0)
    
    # order spatial bins
    index.order <- 
      mols.counts %>% 
      names() %>% 
      gsub("spatial_bin_", "", .) %>% 
      as.integer() %>% 
      order()
    mols.counts %<>% 
      select(index.order %>% all_of())
    
    # set originally plannned `future` session back
    if (exists("orig.plan")) { 
      future::plan(orig.plan)
    } else { 
      # or plan sequential
      future::plan("sequential")
    }
    
  # create single merged `Seurat` obj with Spatial bins
  message(">>> Generating Seurat object..")
  obj <- 
    mols.counts %>% 
    as.sparse() %>%
    CreateSeuratObject(counts = ., assay = assay) %>%
    # add spatial bin name to obj
    AddMetaData(metadata = Cells(.) %>% 
                  gsub("spatial_bin_", "", x = .) %>% as.integer,
                col.name = mols.coord.df %>% 
                  names %>% grep("spatial_", ., value = TRUE))
 
  ## Add FOVs to `Seurat` obj ----
  message(">>> Creating FOVs:", 
          paste0("\n", sep = "\t",
                 c("..centroids","..bounding boxes","..molecules")))
  
  # NOTE: for Resolve flip x coord (ie mutate(x = -x)) in all FOVs
  
  # FOV - boxes ----
  # keep only non-empty bins
  data[["bins.boxes"]] %<>% 
    filter(id %in% obj$spatial_bin) %>%
    rename(cell = id)
  # get cell ids
  cell.ids <- lapply(Cells(obj) %>% length %>% seq(), 
					 function(i) { rep_len(Cells(obj)[i], 4) })
  # add cell ids to boxes
  data[["bins.boxes"]] %<>% mutate(cell = cell.ids %>% unlist)
  bound.boxes <- 
    data[["bins.boxes"]] %>%
    
    { if (spatial.technology == "Resolve") 
        { mutate(., x = -x) } else { (.) }
      } %>%
    CreateSegmentation()
  
  # FOV - box centroids ----
  gc() %>% invisible() # free up memory	
  cents <- 
    custom.lapply(Cells(obj) %>% length %>% seq, function(i) {
      data[["bins.boxes"]] %>% 
        filter(cell == Cells(obj)[i]) %>% 
        select(x, y) %>% 
        colMeans %>%
        rbind %>% 
        as.data.frame %>% 
        mutate(cell = Cells(obj)[i]) 
    }) %>% data.table::rbindlist(.) %>% as.data.frame()
  cents %<>%
    { if (spatial.technology == "Resolve") 
    { mutate(., x = -x) } else { (.) }} %>%
    CreateCentroids()
  
  # Create FOVs (including molecules) ----
  # use binned molecules version
  bound.boxes.data <- list(centroids = cents, boxes = bound.boxes)
  coords <- CreateFOV(coords = bound.boxes.data, 
                      type = c("boxes", "centroids"), 
                      molecules = mols.coord.df %>%
                        as.data.frame %>% 
                        select(contains(c("x", "y", "gene"))) %>%
                        { if (spatial.technology == "Resolve")
                          { mutate(., x = -x) } else { (.) }},
                      assay = assay)
  # only consider the bins that have counts and bounding boxes
  coords <- subset(x = coords, 
                   cells = intersect(x = Cells(x = coords[["boxes"]]),
                                     y = Cells(x = obj)))	
  
  # sanity on fov name
  fov %<>% gsub("_|-", ".", .)
  
  message(">>> Adding FOV..")
  obj[[fov]] <- coords
  
  # add bin.size info to obj
  obj$bin.size <- data[["bin.size"]]
  
  if (update.object) { 
    message("\n", ">>> Updating object..")
    obj %<>% UpdateSeuratObject() } 
  
  message("Object is ready!")    
  return(obj)
  }
