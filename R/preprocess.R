#' Import Data from File
#' @title Import Data
#' @description
#' Imports data from various file formats (.csv, .xlsx, .xls).
#' @param filepath Path to the file to import
#' @return A dataframe containing the imported data
#' @export
#' @examples
#' # data <- importData("data.csv")
importData <- function(filepath) {
  if (grepl(".csv", filepath)) {
    data <- readr::read_csv(filepath)
  } else if (grepl(".xlsx", filepath)) {
    data <- readxl::read_excel(filepath)
  } else if (grepl(".xls", filepath)) {
    data <- readxl::read_excel(filepath)
  } else {
    stop("File type not supported")
  }
    return(data)
}

#' Filter Rows with Missing Values
#' @title Filter Missing Values
#' @description
#' Removes rows with more than a specified threshold of missing values.
#' @param df A dataframe containing the data
#' @param threshold The minimum fraction of non-missing values required (default 0.7)
#' @return A filtered dataframe
#' @export
#' @examples
#' # df <- data.frame(Gene=c("A", "B"), s1=c(1, NA), s2=c(NA, NA), s3=c(2, 3))
#' # filtered <- filterMissingValues(df, threshold=0.5)
filterMissingValues <- function(df, threshold = 0.7) {
  df_numeric <- df %>% select(where(is.numeric))
  df_character <- df %>% select(where(is.character))
  df_numeric <- df_numeric %>% mutate(coverage = rowSums(!is.na(.)) / ncol(df_numeric))
  df_merged <- cbind(df_character, df_numeric) %>% filter(coverage > threshold) %>% 
    select(-coverage)
  return(df_merged)
}

#' Extract Gene Name and Protein Name from Description column
#' @title Extract Gene Name and Protein Name
#' @description
#' Given a dataframe of PD Result with a 'Description' column, this function extracts the
#' Gene Name and Protein Name
#' @param df Data frame containing a 'Description' column
#' @return Data frame with added 'Gene' and 'Protein' columns
#' @export 
#' @examples
#' df <- data.frame(Description = c("sp|P12345|PROT1_HUMAN Protein Name 1 GN=GENE1",
#'                                  "sp|P67890|PROT2_HUMAN Protein Name 2 GN=GENE2"))
#' df_extracted <- extract_gene_protein(df)
extractGeneName <- function(df) {
  if (!"Description" %in% colnames(df)) {
    stop("Data frame must contain a 'Description' column")
  }
  
  # Extract Gene Name
  df$Gene <- sub(".* GN=([^ ]+).*", "\\1", df$Description)
  
  # Extract Protein Name by splitting on OS= and taking the first part
  df$Protein <- sub(" OS=.*", "", df$Description)
  
  return(df)
}


#' Normalize Proteomics Data
#' @title Normalize Data
#' @description
#' Normalizes data using various methods: log2, quantile, or combined approaches.
#' @param df A dataframe containing the data to normalize
#' @param method Normalization method: "log2quantile", "log2", "quantile", or "relative"
#' @return A normalized dataframe
#' @export
#' @examples
#' # df <- data.frame(Gene=c("A", "B"), s1=c(100, 200), s2=c(150, 250))
#' # normalized <- normalizeData(df, method="log2quantile")
normalizeData <- function(df, method = "log2quantile") {
  original_colnames <- colnames(df)
  df_chr <- df %>% select(where(is.character))
  df_num <- df %>% select(where(is.numeric))
  num_colnames <- colnames(df_num)
  if (method == "log2quantile") {
    df_quantile_matrix <- normalize.quantiles(as.matrix(df_num))
    df_quantile <- as.data.frame(df_quantile_matrix)
    colnames(df_quantile) <- num_colnames 
    df_log <- log2(df_quantile)
    result <- cbind(df_chr, df_log)
  } else if (method == "log2") {
    df_log <- log2(df_num + 1)
    result <- cbind(df_chr, df_log) 
  } else if (method == "quantile") {
    df_quantile_matrix <- normalize.quantiles(as.matrix(df_num))
    df_quantile <- as.data.frame(df_quantile_matrix)
    colnames(df_quantile) <- num_colnames 
    result <- cbind(df_chr, df_quantile)
  } else if (method == "relative") {
    df_relative <- as.data.frame(
      apply(df_num, 2, function(x) x / sum(x, na.rm = TRUE) * 100)
    )
    colnames(df_relative) <- num_colnames 
    result <- cbind(df_chr, df_relative)
  } else {
    stop("Method not supported")
  }
  result <- result[, original_colnames] 
  return(result)
}

