#' @title phyloseq.downsample
#' @description
#' Downsample phyloseq object per decile. Phyloseq object should count raw coverage count data (non-normalized)
#' 
#' n viral counts sampled is dependent on the sample library size. For example, if the sample library size = 100 and total viral counts = 10, then when downsampling to 10 counts, only 1 (10/100 * 10 = 1) viral count will be sampled. This is to account for the total library size, and not just the number of viruses identified. 
#' Proportions are rounded to the nearest whole number (e.g 22 viral counts/100 library size * 10 downsample = 2.2  = 2 counts sampled. 
#' Currently, any samples that have 0 counts when downsamples are removed.
#' @param dat non-normalized phyloseq object derived from phyloseq.palmid()
#' @param decile.range downsample decile range. numeric vector of length 2. default is 0 (min) to 90% 
#' @return list:
#' physeq = list of 11 phyloseq objects downsampled to min/max (n = 2) and each decile (n = 9)
#' downsample_depths = read counts used for each downsample
#' @seealso [phyloseq.palmid]
#' @import tidyverse
#' @import phyloseq
#' @export

#TODO: provide support to return 100% of samples, even if there are 0 viral counts after downsampling
phyloseq.downsample <- function(dat,
                                decile.range = c(0, 0.9)
                                ){
  
  ### Sanity check
  # confirm raw counts input: 
  if(!(((dat %>% otu_table)%%1 %>% colSums) == 0) %>% all){ # check if numbers in otu_table are whole numbers
    stop("Input is not counts.")
  } 
  
  
  ### Run
  # sequencing depth (spots) minimum, maximum, and each decile 
  downsample_depths <- quantile(dat %>% sample_data %>% data.frame %>% pull(spots),
                                probs = seq(decile.range[1], decile.range[2], by = .1))
  
  
  # disaggregate the otu table 
  # repeat each OTU by the number of times is appears 
  disaggregated_otu <- map(1:(otu_table(dat) %>% ncol), # for each col (sample) in OTU table
                           
                           ~otu_table(dat)[,.x] %>% .[.!=0] %>% rep(rownames(.), .)
  ) %>% setNames(dat %>% otu_table %>% colnames)
  
  
  
  # progress bar
  pb <- progress::progress_bar$new(
    format = "downsampling [:bar] :percent eta: :eta",
    total = length(downsample_depths),
    show_after = 0,
    clear = T)
  invisible(pb$tick(0))
  
  
  # downsample
  otu_table_downsample <- map(downsample_depths,
                              
                              function(downsample){
                                
                                out <- phyloseq.downsample.help(dat = dat, disaggregated_otu = disaggregated_otu, downsample = downsample)
                                
                                # progress bar
                                pb$tick()
                                Sys.sleep(1/100)
                                
                                return(out)
                                }
                              
                              
                              
                        
                        
  )
  
  
  
  dat_downsample <- map(otu_table_downsample,
                        
                        function(downsampled_otu_table){
                          downsampled_otu_table %<>% otu_table(taxa_are_rows = TRUE)
                          otu_table(dat) <- downsampled_otu_table
                          
                          return(dat)
                        }
  )
  
  return(list(physeq = dat_downsample,
              downsample_depths = downsample_depths))
  

  
  
  
  
  
}

