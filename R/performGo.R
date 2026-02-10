# Function: performGO
# Takes a dataframe with a gene column and performs gene ontology on the
# specified ontology
#
# Args:
#   - df: dataframe containing a column labelled gene with the hgnc gene symbol
#   - ontology = "BP", other options are "CC" for cellular component and
#                "MF" for molecular function
# Returns:
#   - dataframe containing GO results
#' Perform Gene Ontology Analysis
#' @title Gene Ontology Analysis
#' @description
#' Performs gene ontology enrichment analysis on a set of genes using clusterProfiler.
#' @param df A dataframe containing a 'Gene' column with gene symbols
#' @param ontology The ontology to use: "BP" (Biological Process), "CC" (Cellular Component), or "MF" (Molecular Function)
#' @param database The organism database to use (default: org.Rn.eg.db for rat)
#' @param keytype The type of gene identifier (default: "SYMBOL")
#' @return A dataframe containing GO enrichment results
#' @export
#' @examples
#' # df <- data.frame(Gene = c("ACTB", "GAPDH", "TP53"))
#' # go_results <- geneOntology(df, ontology = "BP")
geneOntology <- function(df, ontology = "BP", database = org.Rn.eg.db, keytype = "SYMBOL") {
  genes <- df$Gene
  go_analysis <- as.data.frame(clusterProfiler::enrichGO(gene = genes, OrgDb = database,
                                        keyType = keytype, ont = ontology))
  return(go_analysis)
}