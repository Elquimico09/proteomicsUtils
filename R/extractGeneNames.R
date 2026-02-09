#' Extract Gene Names from Protein Descriptions
#'
#' Extracts gene names from the Description column of a dataframe (typically from proteomics data)
#' and cleans up the description by removing organism information (OS=). Reorders columns to place
#' Gene, Description, and Accession at the beginning.
#'
#' @param df Dataframe containing at least a "Description" column with protein descriptions.
#'   The description is expected to contain gene name information and organism data (OS=).
#'
#' @return A dataframe with an added "Gene" column, cleaned "Description" column, and reordered
#'   columns with Gene, Description, and Accession first.
#'
#' @export
#'
#' @examples
#' # Extract gene names from proteomics data
#' cleaned_data <- extractGeneNames(protein_df)
extractGeneNames <- function(df) {
  df$Gene <- sapply(df$Description, extractGeneName)

  # split description by OS= and take first part
  df$Description <- sapply(df$Description, function(desc) {
    strsplit(desc, " OS=")[[1]][1]
  })

  df <- df  %>%
    dplyr::select(Gene, Description, Accession, dplyr::everything())
  return(df)
}