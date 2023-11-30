#' @title phyloseq.palmid
#' @description
#' Create a phyloseq object for all viruses (palm_id) detected in given SRR runs.
#' @param srr character vector of SRA run IDs to generate phyloseq object. 
#' @param con Serratus SQL database connection 
#' @param metadata list of metadata columns to extract from Serratus table sraruns. Reducing number of columns can decrease query time.
#' @param sample sample column to use from palm_sra2 table as character string. Recommend to use either run_id (each run as sample) or Biosample (each biological sample, including replicates, as sample) **WARNING** currently only tested with run_id.
#' @param taxa one of c("sotu", "palm_id"). **WARNING** currently only implemented/tested with palm_id
#' @return Phyloseq object containig all samples that can be found in the palm_sra2 Serratus table
#' @import tidyverse
#' @import Matrix
#' @import phyloseq
#' @import Biostrings
#' @export


#TODO: test with other samples/taxa combinations 

phyloseq.palmid <- function(srr,
                            con = SerratusConnect(),
                            metadata = c("run", "bio_project", "bio_sample", "spots", "bases", "tax_id", "scientific_name"),
                            sample = "run_id",
                            taxa = "sotu", 
                            qc_filter = T,
                            threads = 1
                            )
  {
  # add mandatory metadata columns
  metadata_input <- metadata
  metadata <- c(metadata_input, "run", "bio_sample", "spots", "tax_id") %>% unique
  
  # query Serratus
  data_import <- phyloseq.palmid.help.query(srr = srr,
                                            con = con,
                                            metadata = metadata,
                                            qc_filter = qc_filter)
  # derive count matrix
  count_mat <- phyloseq.palmid.help.count(data_import$palm,
                                          sample = sample,
                                          taxa   = taxa,
                                          threads = threads
  )
  
  # format taxonomic metadata 
  # note: if using sotu, taxonomic ranks and palm_id are derived from the tax_id and palm sequence associated with that OTUs palm_id 
  #    to be added to the refseq slot? or returned separately? 
  # TODO: add collapsed palm_id 
  tax_metadata <- left_join(data_import$palm %>% select(taxa) %>% distinct,                              # sotu/palm_id found in SRR runs
                          data_import$tax %>% select(palm_id, tax_id, gb_acc, percent_identity, pp_cov), # add taxonomic metadata 
                          by = setNames("palm_id",  # to be matched to the tax_id assigned to the palm_id (i.e. even if sotu, take each otu's palm_id tax assignment)
                                        taxa)) %>%  # for sotu or palm_id, depending on user input
    distinct

  
  # format taxonomic ranks
  tax_mat <- tax_metadata %>% 
    # filter to distinct IDs 
    select(taxa, 
           tax_id) %>% 
    distinct %>%
    # add taxonomy 
    left_join(data_import$tax_ranks %>% select(-row_id),
              by = "tax_id") %>% 
    # format for phyloseq
    select(-c(tax_id)) %>% 
    column_to_rownames(var = taxa)
    
  
  # add palm_id sequence
  tax_refseq_pre <- tax_mat %>% 
    rownames_to_column(var = taxa) %>%
    select(taxa) %>%
    # add sequence
    left_join(data_import$palm_seq %>% select(palm_id, palmprint),
              by = setNames("palm_id", # always use palm_id for match with palm_seq (i.e. even if sotu, take each otu's palm_id palm sequence assignment)
                            taxa)      # use sotu or palm_id from palm_id reads (palmsra2), depending on user input
              ) %>% 
    # format for phyloseq 
    column_to_rownames(var = taxa)
  
  tax_refseq <- tax_refseq_pre$palmprint %>%  # extract palmprint sequence 
    setNames(rownames(tax_refseq_pre)) %>%    # use taxa name (sotu or palm_id) as names
    Biostrings::AAStringSet(use.names = T)    # store as Biostrings object
  
    
  
  
    
  # format metadata 
  sample_mat <- data_import$sra %>% 
    mutate(run_id = run) %>%
    column_to_rownames(var = sample)
  
  

  # create phyloseq object 
  OTU = phyloseq::otu_table(as.matrix(count_mat), taxa_are_rows = TRUE)
  TAX = phyloseq::tax_table(as.matrix(tax_mat))
  samples = phyloseq::sample_data(sample_mat)
  dat_pre <- phyloseq::phyloseq(OTU, TAX, samples)
  
  # add refseq
  dat <- dat_pre %>% phyloseq::merge_phyloseq(tax_refseq)
  
  
  return(dat)
}
