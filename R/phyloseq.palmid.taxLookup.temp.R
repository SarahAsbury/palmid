#' @title taxLookup
#' @description Return a list of all runs processed by Serratus that match
#' tax. Used internally by getVirome.
#  NOT EXPORTED
#' @param tax A taxon defined in NCBI taxonomy. Must be type char.
#' @param con A connection to the Serratus database
#' @return A character vector of SRA accessions.
#' @import dplyr

taxLookup <- function(tax, con) {
  # Get ranking of taxonomic term
  class <- taxize::classification(tax, db = 'ncbi')[[1]]
  # Check if a species was provided
  rank <- class[class$name == tax, 'rank']
  if (is.null(rank)) {
    stop("Error: could not find taxonomic term in NCBI taxonomy database")
  }
  else if (rank == 'species') {
    searchTerms <- class[class$name == tax, 'id']
  }
  else {
    # Collect all child taxa
    taxid <- as.character(class[class$name == tax, 'id'])
    searchTerms <- taxize::downstream(taxid, db ='ncbi', downto = 'species')
    searchTerms <- searchTerms[[1]][,'childtaxa_id']
  }
  # Get all SRA accessions
  query <- tbl(con, "srarun") %>%
    dplyr::filter(tax_id %in% searchTerms) %>%
    dplyr::select(run) %>%
    dplyr::collect()
  return(query$run)
}


# === temporary file 
# this is AH's code, added here for compatability with phyloseq.palmid until taxLookup fun is added to palmid