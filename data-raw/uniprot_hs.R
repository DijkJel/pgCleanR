## code to prepare `uniprot_hs` dataset goes here

up = read.delim('Databases/uniprot_hs_total.tsv', stringsAsFactors = F, header = T)
bm = read.delim('Databases/biomart_hs.txt', stringsAsFactors = F, header = T)

up_hs = up[up$Organism == 'Homo sapiens (Human)',]
#write.table(up_hs, 'uniprot_hs_total.tsv', quote = F, sep = ';', row.names = F)

up_hs$MANE.Select = sapply(up_hs$MANE.Select, function(x){strsplit(x, ';')[[1]][1]})
up_hs$MANE.Select = gsub(' .*', '', up_hs$MANE.Select)
up_hs$MANE.Select = gsub(';', '', up_hs$MANE.Select)
up_hs$Ensembl = sapply(up_hs$Ensembl, function(x){strsplit(x, ' ')[[1]][1]})
up_hs$Ensembl = sapply(up_hs$Ensembl, function(x){strsplit(x, ';')[[1]][1]})
up_hs$Ensembl = gsub(';', '', up_hs$Ensembl)
up_hs$Ensembl = gsub('\\.\\d+', '', up_hs$Ensembl)
up_hs$MANE.Select = gsub('\\.\\d+', '', up_hs$MANE.Select)
up_hs$Ensembl[is.na(up_hs$Ensembl)] = ''
up_hs$MANE.Select[is.na(up_hs$MANE.Select)] = ''

up_ensembl = apply(up_hs[,8:9], 1, function(x){

  if(x[[1]] == ''){
    x[[1]] = x[[2]]
    return(x)
  }
  else {return(x)}
})

up_hs[,8:9] = t(up_ensembl)
up_hs = up_hs[,c(1,2,4,3,5:8)]
names(up_hs) = c('Entry', 'Entry.name', 'Status', 'Protein.names', 'Gene.names', 'Organism', 'Length', 'Ensembl')

up_hs = merge(up_hs, bm, by.x = 'Ensembl', by.y = 'Transcript.stable.ID', all.x = T, all.y = F)
up_hs = up_hs[,c(2:9, 1, 10, 11)]
names(up_hs)[9:11] = c('Transcript.stable.ID', 'RefSeq.MANE.select.ID', 'Entrez.ID')
up_hs[is.na(up_hs)] = ''

write.table(up_hs, 'uniprot_hs_geneAnnotation.txt', sep = ';', row.names = F)


usethis::use_data(uniprot_hs, overwrite = TRUE)
