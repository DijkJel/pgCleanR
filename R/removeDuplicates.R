#' Remove Duplicates from a Protein Group DataFrame
#'
#' This function removes duplicate entries from a given protein group DataFrame.
#'
#' @param pg A protein group DataFrame containing Protein.IDs and Razor...unique.peptides columns.
#' @return A list containing the protein group DataFrame with duplicates removed and a DataFrame of the removed duplicates.
#' @export
#' @examples
#' pg <- data.frame(Protein.IDs = c("A", "A", "B", "C"), Razor...unique.peptides = c(1, 2, 3, 4))
#' removeDuplicates(pg)
#' 
removeDuplicates = function(pg){

  if (any(duplicated(pg$Protein.IDs))){

    dup = pg[duplicated(pg$Protein.IDs), 'Protein.IDs']
    dup = dup[is.na(dup) == F]
    dup_i = lapply(dup, function(x){i = which(pg$Protein.IDs %in% x)})
    rmv = lapply(dup_i, function(i){
      pg_dup  = cbind(i, Score = pg[i, 'Razor...unique.peptides'])
      maxScore = which.max(pg_dup[,'Score'])
      return(pg_dup[-maxScore, 'i'])
    })

    pg$duplicates = F
    pg$duplicates[unlist(dup_i)] = T
    pg$removed = F
    pg$removed[unlist(rmv)] = T

    duplicate_df = pg[pg$duplicates == T,]
    duplicate_df = cbind(duplicate_df[,(ncol(duplicate_df)-1):ncol(duplicate_df)],
                         duplicate_df[,-((ncol(duplicate_df)-1):ncol(duplicate_df))])

    pg = pg[,-((ncol(pg)-1): ncol(pg))]

    pg = pg[-unlist(rmv),]

    msg = paste0(length(dup), ' non-unique entries were found.\n', sum(duplicate_df$removed), ' entries were removed.')
    message(msg)
    return(list(proteinGroups = pg, duplicates = duplicate_df))

  }
  else {
    message('No duplicates found.')
  }

  return(list(proteinGroups = pg))
}
