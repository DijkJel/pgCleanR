#' Check Experimental Design
#'
#' This function checks the consistency and completeness of the experimental design data and protein group data.
#'
#' @importFrom assertthat assert_that
#'
#' @param ed A data frame containing the experimental design data.
#' @param pg A data frame containing the protein group data.
#'
#' @return Throws an error if the experimental design does not match the proteingroups format
#' 
#' @examples
#' ed = createExperimentalDesign(pg)
#' checkExperimentalDesign(ed, pg)
#' @export
checkExperimentalDesign = function(ed, pg){
  
  assert_that(all(colnames(ed) == c('label', 'condition',  'replicate')),
              msg = 'Wrong column names. There should be three columns labeled \'label, condition, and replicate\'')
  
  lfq_columns = paste('LFQ.intensity', ed$label, sep = '.')
  
  assert_that(all(lfq_columns %in% colnames(pg)),
              msg = 'Please check experimentalDesign.\nNot all LFQ columns detected; LFQ columns should be labeled as \'LFQ.intensity.<label>\'.')
  
}


