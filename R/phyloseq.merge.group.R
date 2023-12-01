#' @title phyloseq.merge.group
#' @description
#' Merge samples by some metadata variable. Metadata will be collapsed and separated by |. 
#' The only exception is sequencing depth (spots), which will be both summed (spots) and collapse by | (spots_sep)
#'
#'By default performs merge on bio_samples. Note: it is recommended to merge only after per-run normalization and pre-processing.
#' @param dat phyloseq object 
#' @param group_by character string, metadata variable used to collapse samples  
#' @return phyloseq object with samples merged by group_by variable
#' @seealso [phyloseq.normalize(), phyloseq.palmid]
#' @import tidyverse
#' @import phyloseq
#' @export 


phyloseq.merge.group <- function(dat,
                                 group_by = "bio_sample"){
  out <- dat %>% 
    phyloseq::merge_samples(group_by) %>% 
    t() %>% 
    phyloseq::merge_phyloseq(refseq(dat)) # for some reason merge_samples removes the refseq slot
  
  
  cols <- names(dat %>% sample_data) %>% str_subset(group_by, negate = T)
  metadata <- dat %>% 
    sample_data %>% 
    data.frame %>% 
    group_by(!!!syms(group_by)) %>% 
    mutate(across(.cols  = all_of(cols),
                  .fns   = ~paste(.x %>% unique, collapse = "|")
    )
    mutate(spots_sep = spots,
           spots     = spots_sep %>% str_split("\\|") %>% map(~.x %>% as.numeric %>% sum) %>% unlist
    ) %>% distinct %>%
    mutate(rownames = get(group_by)) %>% 
    column_to_rownames(var = "rownames")
  
  
  sample_data(out) <- metadata
  
  return(out)
  
}
