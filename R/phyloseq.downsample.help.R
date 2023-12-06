#' @title phyloseq.downsample.help
#' @description
#' Helper function for phyloseq.downsample. Performs the downsample for a single dataset
#' 
#' @param dat phyloseq object derived from phyloseq.palmid
#' @param disaggregated_otu OTU table converted as list with taxa repeated n times. Produced by phyloseq.downsample.
#' @param downsample Numeric. Downsample to library depth x. 
#' @return downsampled otu table (matrix)
#' @seealso [phyloseq.palmid, phyloseq.downsample]
#' @import tidyverse
#' @import phyloseq
#' @import magrittr
#' @export




phyloseq.downsample.help <- function(dat, disaggregated_otu, downsample){
  
  
  print("Start sample_n")
  Sys.time() %>% print
  # sample viral count n
  # sample n number of viral counts depending on ratio between: total viruses per sample:sample library depth
  sample_n <- dat %>% sample_data %>% data.frame %>% 
    mutate(downsample = downsample,
           total_viral_counts = dat %>% otu_table %>% colSums()
    ) %>%
    dplyr::rowwise() %>% 
    mutate(sample_n = total_viral_counts/spots * min(downsample, spots)) %>% # downsize to downsamples depth or library size, whichever is smaller
    pull(sample_n, name = ) %>% 
    round(digits = 0) %>% # round to the nearest whole number
    setNames(dat %>% sample_data %>% data.frame %>% rownames())
  
  # sanity check
  if((names(sample_n) != names(disaggregated_otu)) %>% any) stop("Unexpected error. Names do not match between derived data.")
  
  
  
  print("Start Downsample")
  Sys.time() %>% print
  
  # sample n number of viral counts from each sample
  downsampled_list <- map2(disaggregated_otu,
                           sample_n,
                           
                           ~.x %>% sample(size = .y, replace = F)
                           
  ) 
  
  print("Start matrix formatting/rbind")
  Sys.time() %>% print
  
  
  # format downsampled counts as matrix
  downsampled_table <- map2(downsampled_list,
                            names(downsampled_list),
                            
                            ~.x %>% table %>%
                              as.matrix %>% 
                              magrittr::set_colnames(.y)) %>% 
    keep(\(x) nrow(x) > 0) %>% # only keep matrices with at least 1 taxa sampled
    map(~.x %>%
          t() %>%
          as("sparseMatrix") 
    ) %>%
    
    # bind rows
    purrr::reduce(rbind.sparse) %>%

    # final format
    Matrix::t() %>% 
    as.matrix %>%
    replace(is.na(.), 0) # replace NA with 0 viral counts

  
    
  print("Fini")
  Sys.time() %>% print
  
  
  
  return(downsampled_table)
}





