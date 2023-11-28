#' @title phyloseq.palmid.help.query
#' @description
#' Takes a dataframe queried from Serratus table palm_sra2 and derives a taxa x sample count matrix. Counts are aggregated from coverage per taxa per sample. 
#' @param palm character vector of SRA run IDs to generate phyloseq object. 
#' @param sample user-defined sample column. Recommend to use either run_id (each run as sample) or Biosample (each biological sample, including replicates, as sample)
#' @param taxa one of c("sotu", "palm_id").
#' @return Taxa x sample count matrix. Taxa are rows. Columns are samples.
#' @import tidyverse
#' @import Matrix


phyloseq.palmid.help.count <- function(palm,
                                       sample     = "run_id",
                                       taxa       = "sotu"
                                       )
{
  
  #TODO: add input sanity checks
  # sample column found in palm 
  # sample column is one of recommended columns (else return warning)
  # taxa column is one of sotu or palmid (else return error) 
  
  # calculate total number of reads per OTU
  reads <- suppressMessages(
    palm %>% 
      dplyr::select(all_of(c(sample, taxa, "coverage"))) %>% 
      
      dplyr::group_by(!!!syms(sample), !!!syms(taxa)) %>% 
      dplyr::summarize(total_reads = sum(coverage)) %>% 
      ungroup
  )
  
  # count matrix
  count_mat <- reads %>%
    # partition to save memory
    dplyr::arrange(!!!syms(sample)) %>% 
    tibble::rownames_to_column(var = "row_id") %>% 
    dplyr::mutate(partition = (as.numeric(row_id)/1000) %>% ceiling) %>% 
    dplyr::group_by(!!!syms(sample)) %>%  # take the minimum partition assignment for each sample so that all rows for a sample are processed together 
    dplyr::mutate(partition = min(partition)) %>% 
    ungroup %>% 
    split(.$partition) %>% 
    
    # count matrix
    furrr::future_map(~.x %>% 
                        dplyr::select(-all_of(c("partition", "row_id"))) %>% 
                        tidyr::pivot_wider(names_from = taxa, 
                                           values_from = total_reads) %>% 
                        tibble::column_to_rownames(sample) %>%
                        replace(is.na(.), 0) %>%                                           # replace NA with 0
                        as.matrix %>% 
                        as("sparseMatrix")
                        
    ) 
  
  # return(count_mat) # for testing
  
  # bind partitioned matrices
  # chunk size = 10
  count_final <- count_mat %>% reduce(rbind.sparse)


  # phyloseq formatting
  out <- count_final %>%
    Matrix::t() 


  return(out)
}