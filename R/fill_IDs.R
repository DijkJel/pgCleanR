#' Fill IDs and gene names based on UniProt
#'
#' This function fills in the IDs and gene names of proteins based on UniProt data.
#'It takes in a proteinGroups data frame with protein IDs in the 'Protein.IDs' column and gene names in the 'Gene.names' column.
#'It tries to assign a single, unique proteinID and corresponding official gene symbol to each entry of pg.
#'If a protein ID in the input data frame matches a reviewed protein in UniProt, the corresponding gene name is used.
#'If the protein ID matches an unreviewed protein in UniProt, the function attempts to find a reviewed protein with the same gene name, or failing that, uses the protein ID as the gene name.
#'If a protein ID does not match any protein in UniProt, the function uses the protein ID as the gene name.
#'
#'
#'
#' @param pg A maxquant proteinGroups.txt file as data frame with protein IDs in the 'Protein.IDs' column and gene names in the 'Gene.names' column
#' @param uniprot A UniProt data frame with all with reviewed and unreviewed human protein sequences. Comes with the package.
#' @param species The species to be used for UniProt data, default is 'Hs'
#'
#' @return A proteinGroups data frame with reviewed protein IDs and official gene symbols, where possible.
#' @export
#'
#' @examples
#' \dontrun{fill_IDs(my_data, uniprot_data)}
#'

fill_IDs = function(pg, uniprot, species = 'Hs'){

  line = 0
  no_match = 0
  nrow_in = nrow(pg)

  pg = pg[pg$Potential.contaminant != '+' & pg$Reverse != '+',]
  pg = pg[!is.na(pg$Razor...unique.peptides),]
  pg = pg[pg$Razor...unique.peptides > 1,]
  
  no_peptides = sapply(pg$Peptide.counts..razor.unique., function(x){
    n = as.numeric(strsplit(x, ';')[[1]])
    n = which(n == max(n))
    return(n[length(n)])
  })
  
  
  extract_proteinId = function(n, protID){
    
    pid = strsplit(protID, ';')[[1]][1:n]
    pid = paste(pid, collapse = ';')
    return(pid)
  }
  
  pg$Protein.IDs = mapply(extract_proteinId, no_peptides, pg$Protein.IDs)

  uniprot_rev = uniprot[uniprot$Status == 'reviewed',]
  uniprot_rev$Gene.names = sapply(uniprot_rev$Gene.names, function(x){strsplit(x, ';')[[1]][1]})
  uniprot_unrev = uniprot[uniprot$Status == 'unreviewed',]

  ten_pct = round(nrow(pg)/10)
  next_pct = ten_pct
  pct = 10


  rev_pid = apply(pg, 1, function(x){

    pid = x[['Protein.IDs']]
    pid = sapply(pid, function(x){strsplit(x, ';')})[[1]]
    gene_names = x[['Gene.names']]

    line <<- line + 1

    if (line > next_pct){
      msg = paste0(pct, '% done')
      message(msg)
      pct <<- pct + 10
      next_pct <<- next_pct + ten_pct
    }

    if(any(pid %in% uniprot_rev$Entry)){
      gn = uniprot_rev[uniprot_rev$Entry == pid[which(pid %in% uniprot_rev$Entry)[1]], 'Gene.names']
      gn = strsplit(gn, ' ')[[1]][1]
      id = pid[which(pid %in% uniprot_rev$Entry)[1]]

      out = list(id = id, gene_name = gn)
    }
    else if (any(pid %in% uniprot_unrev$UniProtKB.Gene.Name.ID)) {

      up = uniprot_unrev[which(up$Entry %in% pid),]
      ensembl_ids = up$Gene.stable.ID

      if (any(ensembl_ids %in% uniprot_rev$Gene.stable.ID)){
        uniprot_id = uniprot_ids[uniprot_ids %in% uniprot_rev$Entry, 'UniProtKB.Gene.Name.ID'][1]
        gn = uniprot_rev[uniprot_rev$UniProtKB.Gene.Name.ID == uniprot_id, 'Gene.names']
        gn = strsplit(gn, ' ')[[1]][1]
      }
      else if (any(ensembl_ids %in% uniprot_unrev$Gene.stable.ID)) {
        uniprot_ids = uniprot_unrev[which(ensembl_ids %in% uniprot_unrev$Gene.stable), 'UniProtKB.Gene.Name.ID']
        uniprot_id = names(which.max(table(uniprot_ids)))
        gn = uniprot_unrev[uniprot_unrev$UniProtKB.Gene.Name.ID == uniprot_id,'Gene.names']
        gn = strsplit(gn, ' ')[[1]][1]
      }
      out = list(id = uniprot_id, gene_name = gn)
      if(is.na(out$id)){print(line)}
    }
    else {
      no_match <<- no_match + 1
      uniprot_id = pid[1]
      gn = pid[1]
      out = list(id = uniprot_id, gene_name = gn)
    }
    return(out)
  })


  pg$Protein.IDs = sapply(rev_pid, function(x){x[[1]]})
  pg$Gene.names = sapply(rev_pid, function(x){x[[2]]})

  nrow_out = nrow(pg)
  difference = nrow_in-nrow_out

  msg = paste0('\n',difference, '/', nrow_in, ' entries were removed.')
  msg2 = paste0(no_match, ' gene names could not be identified. Protein IDs used as gene name in ', no_match, ' cases.')
  message(msg)
  message(msg2)

  return(pg)
}
