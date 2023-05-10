#' Clean protein groups data
#'
#' This function takes a protein groups data frame \code{pg} and cleans it by
#' filling missing UniProt IDs using the \code{uniprot} object, removing
#' duplicate entries, and adding gene information. The resulting data frame is
#' returned.
#'
#' @param pg A data frame of protein groups, with one row per group and one or
#'   more protein IDs in each group.
#' @param uniprot An object of class \code{uniprot} containing UniProt ID
#'   mappings for the proteins in \code{pg}.
#' @param species A character string specifying the species to use for gene
#'   information. The default is \code{'Hs'} (human).
#' @return A data frame of cleaned protein groups.
#' @seealso \code{\link{fill_IDs}}, \code{\link{removeDuplicates}},
#'   \code{\link{addGeneInfo}}
#' @export
#' @examples
#' ## Load protein groups and UniProt mappings
#' data(pg)
#' data(uniprot_hs)
#' ## Clean protein groups
#' cleanPg(pg, uniprot_hs)
#' 
cleanPg = function(pg, uniprot, species = 'Hs'){
  
  message('Fill protein IDs and gene names...')
  pg2 = fill_IDs(pg, uniprot)
  message('\nRemove possible duplicate entries...')
  pg3 = removeDuplicates(pg2)
  message('\nAdd gene IDs...')
  pg4 = addGeneInfo(pg3$proteinGroups, uniprot)
  return(pg4)
}