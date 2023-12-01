#' @title phyloseq.cladeogram
#' @description
#' Add a cladeogram to a phyloseq object phylo-tree slot object for all Serratus IDs found in input 
#' @param dat phyloseq object derived from phyloseq.palmid()
#' @param con Serratus connection
#' @return phyloseq object with phyloseq-tree slot filled with taxa 
#' @seealso [phyloseq.palmid, SerratusConnect()]
#' @import tidyverse
#' @import phyloseq


#TODO: not all taxa show-up in palm_graph - why? Halted function development.
phyloseq.cladeogram <- function(dat, 
                                con = SerratusConnect()){
  
  # get Serratus palm_id pairwise sequence identity 
  query_taxa <- dat %>% phyloseq::tax_table() %>% rownames # taxa to found in phyloseq object  
  
  taxa_pident <- tbl(con, "palm_graph") %>%       # query Serratus
    dplyr::select(palm_id1, palm_id2, pident) %>% 
    dplyr::filter(palm_id1 %in% query_taxa & palm_id2 %in% query_taxa)
  
  
  # convert to distance matrix 
  
  return(out)
}