#' Add Gene IDs to a proteinGroups Dataset
#'
#' Merge a proteomics dataset with a UniProt database to add gene IDs (Ensembl, RefSeq, Entrez)
#' for the identified proteins.
#'
#' @param pg A data frame with protein IDs under the column "Protein.IDs". Needs to be cleaned with 'Fill_IDs'
#' @param uniprot A data frame with UniProt protein information, with protein IDs under the column "Entry". Also needs to contain
#' @return A data frame with added gene IDs (Ensembl_ID, RefSeq_ID, Entrez_ID).
#'
#' @examples
#' \dontrun{addGeneIDs(pg, uniprot)}
#'
#' @export
addGeneInfo = function(pg, uniprot){

  gene_ids = merge(pg, uniprot, by.x = 'Protein.IDs', by.y = 'Entry', all.x = T, all.y = F)
  nc = ncol(gene_ids)
  pg = gene_ids[,-c((nc-9):(nc-4),(nc-2))]
  pg = pg[!duplicated(pg$Protein.IDs),]

  colnames(pg)[(ncol(pg)-2):ncol(pg)] = c('Ensembl_ID', 'RefSeq_ID', 'Entrez_ID')
  pg$RefSeq_ID = gsub('\\.\\d+', '', pg$RefSeq_ID)

  pg[is.na(pg)] = ''
  n_missing = sum(pg$Ensembl_ID == '')
  pct_missing = round(n_missing / nrow(pg) * 100, digits = 2)
  msg = paste0(n_missing, ' entries could not be matched (', pct_missing, '%).')
  message(msg)

  names(pg)[names(pg) == 'Gene.names.x'] = 'Gene.names'
  return(pg)
}
