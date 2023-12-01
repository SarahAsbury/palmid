#' @title phyloseq.rarefy
#' @description
#' Rarefy phyloseq object per decile
#' @param dat phyloseq object derived from phyloseq.palmid()
#' @return list of 11 phyloseq objects rarefied to min/max and each decile 
#' @seealso [phyloseq.palmid]
#' @import tidyverse
#' @import phyloseq


phyloseq.rarefy <- function(dat){
  
  # sequencing depth (spots) minimum, maximum, and each decile 
  depths <- quantile(dat %>% sample_data %>% data.frame %>% pull(spots),
                     probs = seq(0, 1, by = .1))
  
  # proportion of samples represented by rarefied depth 
  map(depths,
      
      function(x){
      prop <- dat %>% sample_data %>% data.frame %>% 
        mutate(rarefy_prop = x/spots) %>% 
        pull(rarefy_prop)
      
      # divide otu table by rarefied proportion factor
      rarefy <- (dat %>% otu_table %>% t()/prop) %>% t()
      
      # update phyloseq
      out <- dat
      otu_table(out) <- rarefy
      
      return(out)

      )
  }
  
  
  
  
}