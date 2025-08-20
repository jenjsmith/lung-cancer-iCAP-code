
# update column names using column definitions
#   df: <tibble>
#   column_names: <character>, current, new, current, new, ...
rename_columns <- function(df, column_names){
  
  # format input list
  m <- matrix(column_names, nrow=length(column_names)/2, ncol=2, byrow=T)
  colnames(m) <- c('from', 'to')
  m <- m |> tibble::as_tibble()
  
  # only include if present
  m <- m[m$from %in% names(df),]
  
  # rename columns
  df <- df |> dplyr::rename_at(dplyr::vars(m$from), function(x) m$to)
  
  return(df)
}

# rename columns to standardize names
#   - rename ID to _id
#   - remove camelCase
#   - remove trailing periods
#   - replace periods with underscores
#   - all lowercase
#   - remove '___' for repeated column names
#   - replace '/' with underscore
read_excel_clean <- function(fp, range=NULL, include_rows=NULL, sheet=1, skip=0){
  
  column_names.rename.pre <-c(
    'PreCytePlate', 'icap_plate',
    'PreCyte group', 'project_group',
    'iCAP plate', 'icap_plate',
    'iCAP batch', 'icap_batch',
    'iCAP well', 'icap_well'
  )
  
  column_names.rename.post <-c(
    'bar_code', 'barcode',
    'expbatch', 'icap_batch',
    'exp_batch', 'icap_batch',
    'i_cap_group', 'icap_group',
    'i_cap_batch', 'icap_batch',
    'i_cap_plate', 'icap_plate',
    'i_cap_well', 'icap_well',
    'smoking', 'smoking_status'
    
  )
  
  # read excel with minimal name repair
  df <- readxl::read_excel(fp, .name_repair = 'minimal', range=range, sheet = sheet, skip=skip)
  
  # drop duplicate columns
  df <- df[!duplicated(as.list(df))]
  
  # include rows
  if (!is.null(include_rows)){
    df <- df[include_rows,]
  }
  
  # drop empty columns
  df <- df[,colnames(df) != ""]
  
  # rename prior to substitutions
  df <- df |> rename_columns(column_names.rename.pre)
  
  # clean up column names
  df <- df |> 
    dplyr::rename_all(~gsub("ID$", "_id", .)) |>
    dplyr::rename_all(~gsub("([a-z])([A-Z])", "\\1.\\L\\2", ., perl = TRUE)) |>  # rename camel case
    dplyr::rename_all(~gsub('\\.$', "", .)) |>
    dplyr::rename_all(~gsub('\\.| |\\/|\\(|\\)|\\-', "_", .)) |>
    dplyr::rename_all(~tolower(.)) |>
    dplyr::rename_all(~gsub('_{1,}', '_', .)) |>
    dplyr::rename_all(~gsub('_$', '', .)) |>
    dplyr::rename_all(~gsub('\\%', 'percent', .)) |>
    dplyr::rename_all(~gsub('\\#', 'number', .))
  
  # rename using standard definition
  df <- df |> rename_columns(column_names.rename.post)
  
  return (df)
}


# fixes rstudio/knitr image rendering
#   usage: knitr::include_graphics(ggimage(gg, 'plot1'))
ggimage <- function(gg, fn, width=14, height=8, scale=2, base_path='img'){
  fp <- file.path(base_path, paste(fn, 'png', sep = '.'))
  ggsave(fp, plot = gg, width = width, height=height, scale=scale)
  knitr::include_graphics(fp)
}

