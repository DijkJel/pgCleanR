library(devtools)

pg = read.delim('C:/Data/Jelmer/NKI/pg_package/Data Sets/THRONCAT_proteinGroups.txt', stringsAsFactors = F)
uniprot_hs = read.delim('C:/Data/Jelmer/NKI/pg_package/Databases/uniprot_hs_geneAnnotation.txt', stringsAsFactors = F, header = T, sep = ';')

use_git()
use_r('fill_IDs')
use_r('data')
use_mit_license()

usethis::use_data(pg)
usethis::use_data(uniprot_hs)
usethis::use_data_raw("uniprot_hs")
usethis::use_build_ignore('setup.R')

use_r('addGeneInfo')
use_r('removeDuplicates')


pg2 = fill_IDs(pg, uniprot_hs)
pg3 = addGeneInfo(pg2, uniprot_hs)
pg4 = removeDuplicates(pg3)
