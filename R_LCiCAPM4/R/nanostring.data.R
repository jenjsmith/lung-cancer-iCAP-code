#' @omit_threshold omit samples based on median of housekeeping genes.
load_ns_plexset <- function(
  xls, 
  metadata, 
  rcc_update=list(),
  omit_regex = NULL,
  omit_precyte_uids=c(),
  omit_threshold=200
){
  
  require(SummarizedExperiment)
  
  # load nCounter normalized counts
  ns_data <- read_excel_clean(xls)
  
  # recode column names if rcc_update values are provided
  if (length(rcc_update) > 0)
    colnames(ns_data) <- colnames(ns_data) |> recode_character_vector(rcc_update)
  
  ## omit partial plates
  if ( ! is.null(omit_regex) ){
    
    ns_data <- ns_data[,-grep(omit_regex, colnames(ns_data))]
  
  }else{

    # omit based on housekeeping signal threshold
    # detect partial/empty rows and omit 
    rcc_idx_omit <- ns_data |> 
      dplyr::filter(class_name=="Housekeeping") |>
      dplyr::select(starts_with('set')) |>
      mutate(across(where(is.character), as.numeric)) |>
      summarise(across(where(is.numeric), mean)) |>
      tidyr::gather('rcc_uid', 'mean') |>
      mutate(rcc_idx = row_number()) |>
      dplyr::filter(mean < omit_threshold) |>
      pull(rcc_idx)
      
    # split, filter, combine
    if (length(rcc_idx_omit) > 0){
      
      # ns_meta <- ns_data |> 
      #   dplyr::select(!starts_with('set'))
      # ns_counts <- ns_data |> 
      #   dplyr::select(starts_with('set'))    
      # ns_data <- bind_cols(ns_meta, ns_counts[,-rcc_idx_omit])
      
      ns_data <- bind_cols(
        ns_data[,-grep('^set', names(ns_data))],
        ns_data[, grep('^set', names(ns_data))][,-rcc_idx_omit]
      )
    }    

  }
    
  # update column names of samples
  #   expected format: set_a_a02_20211013_pre_cyte_20211012_01_2_02_rcc
  #   
  names(ns_data)[grep('^set', names(ns_data))] <- unlist(
    lapply(names(ns_data)[grep('^set', names(ns_data))], 
           FUN=function(x) toupper(paste(unlist(strsplit(x, '_'))[c(8,3)], collapse='_'))
    )
  )
  
  # create nanoString plate names
  ns_id <-
    paste(metadata$ns_plate_uid |> stringr::str_pad(pad="0", width=2),
          gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", metadata$ns_well_position, perl = TRUE),
          sep='_'
    )
  
  metadata <- metadata |> tibble::add_column(ns_id, .after='barcode')
  
  # subset metadata for cases where nCounter flagged outlier samples were removed
  metadata <- metadata[metadata$ns_id %in% names(ns_data),]
  
  # add rcc_filename 
  metadata <- metadata |> 
    tibble::add_column(
      rcc_filename=as.character(ns_data[1, metadata$ns_id]), 
      .after = 'ns_id'
    )
  
  # create SummarizedExperiment object
  #   ns_data is read as a string, and first requires conversion to double
  #   omit rows 1,2 which contain rcc file information
  x <- ns_data[-1:-2, metadata$ns_id]
  m <- matrix(as.double(unlist(x)), ncol=ncol(x), dimnames = list(ns_data[-1:-2,]$probe_name, metadata$precyte_uid))
  ex <- SummarizedExperiment::SummarizedExperiment(
    assays=list(normalized_counts=m),
    colData = metadata,
    rowData = ns_data[-1:-2,1:19] |> 
      dplyr::select(-comments, -class_source, -annotation) |>
      dplyr::rename(gene_name=probe_name)
  )
  # assert colData and rowData correspond to existing names
  assertthat::assert_that(all(rownames(colData(ex)) == colData(ex)$precyte_uid))
  assertthat::assert_that(all(rownames(rowData(ex)) == rowData(ex)$gene_name))
  
  # add log2 normalized counts
  assay(ex, 'log2_normalized_counts') <- log2(assay(ex, 'normalized_counts'))
  
  # apply filter omit_precyte_uids
  if (length(omit_precyte_uids) > 0){
    ex <- ex[,!ex$precyte_uid %in% omit_precyte_uids]
  }
    
  return (ex)
}

#' @description substitute values in a character vector
#' @vec character vector
#' @values list(c('pattern', 'replacement'), ...)
#' @returns character vector with replacements
recode_character_vector <- function(vec, values=list()){
  for (val in values){
    vec <- recode_character_vector(sub(val[1], val[2], vec), tail(values, -1))
  }
  return(vec)
}



