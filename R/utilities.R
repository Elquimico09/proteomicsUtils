

# Function: performGO
# Takes a list of dataframes and performs gene ontology on each dataframe and returns
# a list of gene ontology results
#
# Args:
#   - data: list of dataframes containing each dataframe to perform GO Analysis on
#
# Returns:
#   - goresult: list containing dataframes with the go_results
#' Perform Gene Ontology on Multiple Dataframes
#' @title Perform GO Analysis on List
#' @description
#' Takes a list of dataframes and performs gene ontology analysis on each,
#' across all three ontologies (BP, MF, CC).
#' @param data A list of dataframes, each containing a 'Gene' column
#' @return A list of lists containing GO results for each dataframe and ontology
#' @export
#' @examples
#' # df1 <- data.frame(Gene = c("ACTB", "GAPDH"))
#' # df2 <- data.frame(Gene = c("TP53", "MYC"))
#' # results <- performGO(list(df1, df2))
performGO <- function(data) {
  ontologies <- c("BP", "MF", "CC")
  goresult <- list() # Empty list to store the final results
  for (i in seq_along(data)) {
    gene_ontology_data <- data[[i]]
    gene_ont <- list() # Empty list to store result for each go before merging with goresult
    n <- 1
    for (onto in ontologies) {
      gene_ont[[n]] <- geneOntology(gene_ontology_data, ontology = onto, database = org.Rn.eg.db)
      n <- n + 1    # move on to the next ontology
    }
    goresult[[i]] <- gene_ont
    }
  return(goresult)
}





















  

  












