#' Load Protein Groups
#'
#' This function loads protein group data from a user-selected file.
#'
#' @return A data frame containing the protein group data.
#'
#' @importFrom utils read.delim
#'
#'
#' @examples
#' \dontrun{loadProteinGroups()}
#' 
#' @export
loadProteinGroups = function(){
  
  pg = read.delim(file.choose(), stringsAsFactors = F)
  return(pg)
}
