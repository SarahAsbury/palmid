#' @title phyloseq.palmid.help.query
#' @description
#' Takes an input of SRR run ids. Returns the 4 data tables from Serratus required to generate phyloseq object. Helper function for phyloseq.palmid. 
#' @param srr character vector of SRA run IDs to generate phyloseq object. 
#' @param con a connection to the Serratus database. use SerratusConnect by default 
#' @param metadata character vector of sraruns table columns that should be extracted for each SRR sample 
#' @return 4 filtered Serratus tables: 
#' palm = palm_id counts (coverage) for each SRR sample. Derived from Serratus table palm_sra2.
#' sra = SRR sample metadata. Derived from Serratus table srarun. 
#' tax = contains tax_id for each palm_id. Derived from Serratus table palm_gb.
#' tax_ranks = maps tax_id to taxonomic rank. Derived from Serratus table tax_lineage.
#' @seealso [SerratusConnect()]
#' @import tidyverse
#' @import progress

phyloseq.palmid.help.query <- function(srr, 
                                       con = SerratusConnect(),
                                       metadata = c("run", "bio_project", "bio_sample", "spots", "bases", "tax_id", "scientific_name"),
                                       qc_filter = T
)
{
  

  ### sanity checks - input 
  srr_regex <- "[E|S|D]RR[:digit:]{6,7}"
  if(!(srr %>% str_detect(srr_regex) %>% all))
  {srr %>% str_subset(srr_regex, negate = T) %>% print
    warning("SRR runs are in unexpected format. Example of expected format: 'SRR003659'")
  }
  
  
  ### processing
  # save results to list 
  out <- list()
  
  # progress bar
  assign("pb", 
         progress::progress_bar$new(
           format = "querying serratus [:bar] :percent eta: :eta",
           total = 5,
           show_after = 0,
           clear = T),
         envir = .GlobalEnv)
  invisible(pb$tick(0))
  
  
  # - palmid 
  # query palm_id
  out$palm <- tbl(con, "palm_sra2") %>% 
    dplyr::filter(run_id %in% srr) %>%              # run id in user input 
    dplyr::filter(coverage != 0) 
  
  if(qc_filter == T){
    out$palm <- out$palm %>% filter(qc_pass == TRUE)
    
  }
  out$palm <- out$palm %>% dplyr::collect()
  
  
  
  # - metadata
  # get SRR with palm_id
  srr_palm <- out$palm$run_id # speed-up by only querying sraruns that have palm_id data
  
  # query metadata
  out$sra <- tbl(con, "srarun") %>% 
    dplyr::filter(run %in% srr_palm) %>% 
    dplyr::select(metadata) %>% 
    dplyr::collect()
  assign("pb", pb$tick(), envir = "Global.Env")
  # pb$tick()
  Sys.sleep(1/100)
  


  # - taxa
  # get palm_ids
  query_taxa <- out$palm %>% 
    dplyr::select(palm_id, sotu) %>% # get all taxonomic data for all match of either: palm_id or sotu
    unlist %>%
    unique
  
  pb$tick()
  Sys.sleep(1/100)
  
  
  
  # query palmid_id tax_id
  out$tax <- tbl(con, "palm_gb") %>% 
    dplyr::filter(palm_id %in% query_taxa) %>% 
    dplyr::collect()
  pb$tick()
  Sys.sleep(1/100)
  
  # get viral taxonomy id
  query_tax_id <- out$tax %>% 
    dplyr::pull(tax_id)

  #query tax_id taxonomic ranks 
  out$tax_ranks <- tbl(con, "tax_lineage") %>% 
    dplyr::filter(tax_id %in% query_tax_id) %>% 
    dplyr::collect()
  pb$tick()
  Sys.sleep(1/100)
  
  # query for palm_id sequence 
  out$palm_seq <- tbl(con, "palmdb2") %>% 
    dplyr::filter(palm_id %in% query_taxa) %>% 
    dplyr::collect()
  pb$tick
  Sys.sleep(1/100)

  ### ctrl
  message(sprintf("\n\n%s/%s SRR ids were found in palm_id database (QC filter = %s).", length(srr_palm %>% unique), length(srr %>% unique), qc_filter))
  
  
  
  
return(out)
  
  
}

