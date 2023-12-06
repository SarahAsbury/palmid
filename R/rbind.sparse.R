#' @title rbind.sparse
#' @description
#' Bind rows of 2 sparse matrices. Add names.
#' @param list Two sparse matrices to bind 
#' @return Single sparse matrix derived from all row-binded sparse matrices in list.
#' @import tidyverse
#' @import Matrix
#' @export 
#' 

rbind.sparse <- function(x,y){
  
  
  # identify all unique column names from both matrices
  all_columns <- union(colnames(x), colnames(y))
  
  # create sparse matrices with missing columns
  mat_filled <- map(list(x, y),
                    
                    
                    ~cbind(
                      # bind existing matrix
                      .x, 
                      # and empty columms 
                      Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), 
                                   dims = c(nrow(.x), 
                                            setdiff(all_columns, colnames(.x)) %>% length),
                                   dimnames = list(NULL, 
                                                   setdiff(all_columns, colnames(.x))), 

                    )
  )
  )
  
  #reorder columns in matrix y to match columns in matrix x
  column_order <- colnames(mat_filled[[1]])
  mat_filled[[2]] <- mat_filled[[2]][, column_order, drop = F]

  # rbind
  out <- mat_filled %>% purrr::reduce(~rbind2(.x,.y))
 
  return(out) 
}






