#' @title phyloseq.summary
#' @description
#' Get common summary statistics per OTU and per sample from phyloseq object. Expects metadata derived from physeq.palmid().
#' @param dat phyloseq object derived from physeq.palmid(). taxa are rows, samples are columns 
#' @return Two lists:  
#' viral_summary  
#' - taxa_told_detected = total number of counts detected for virus across all samples (same units as otu_table(dat)) 
#' - taxa_n_samples = count of number of samples a taxa is detected within
#' - taxa_prevalence = count of number of samples taxa is detected withim, divided by all samples in the study (taxa_n_samples/# sample in phyloseq)
#' 
#' sample_summary  
#' - total_detected = total viral counts detected within sample (same units as otu_sable(dat))
#' - viral_community_n = total number of unique taxa detected per sample
#' 
#' @seealso [phyloseq.normalize(), phyloseq.palmid()]
#' @import tidyverse
#' @import phyloseq
#' @export 


phyloseq.summary <- function(dat){

  viral_summary <- data.frame(
    taxa_total_detected = dat %>% otu_table() %>% colSums(),
    taxa_n_samples      = (dat %>% otu_table() %>% replace(.!=0, 1) %>% colSums)         # number of samples taxa detected within
  ) %>% 
    mutate(taxa_prevalence     = taxa_n_samples/(nrow(dat %>% otu_table))                # number of samples taxa detected within/total number of samples
    ) %>% 
    cbind(dat %>% tax_table()) %>% 
    rownames_to_column(var = "sotu")
  
  
  sample_summary <- data.frame(
    total_detected      = dat %>% otu_table() %>% rowSums(),
    viral_community_n   = (dat %>% otu_table() %>% replace(.!=0, 1) %>% rowSums)                            # total number of unique taxa detected per sample
    
  ) %>% 
    cbind(dat %>% sample_data())
  
  
  return(list(viral_summary = viral_summary,
              sample_summary = sample_summary))
  
}
