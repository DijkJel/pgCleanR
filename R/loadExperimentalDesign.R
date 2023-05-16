#' Load Experimental Design
#'
#' This function loads the experimental design data from a user-selected file.
#'
#' @return A data frame containing the experimental design data.
#'
#' @importFrom assertthat assert_that
#' @importFrom utils read.delim
#'
#' @examples
#' \dontrun{
#' loadExperimentalDesign()
#' }
#' 
#'
#' @keywords data
#' @export
loadExperimentalDesign = function(){
  
  
  ed = read.delim(file.choose(), stringsAsFactors = F)
  assert_that(all(colnames(ed) == c('label', 'condition',  'replicate')),
              msg = 'Wrong column names. There should be three columns labeled \'label, condition, and replicate\'')
  
  return(ed)
}

