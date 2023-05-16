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
usethis::use_testthat(3)

use_r('addGeneInfo')
use_r('removeDuplicates')
use_r('createExperimentalDesign')
use_r('performImputation')
use_r('performDEP')
use_r('cleanPg')
use_r('imputeMultiple')
use_r('getImputationStats')
use_r('performMultipleImputation')
use_r('loadProteinGroups')
use_r('loadExperimentalDesign')
use_r('checkExperimentalDesign')


use_package('utils')
use_package('stats')
use_package('SummarizedExperiment')
use_package('matrixStats')
use_package('assertthat')

use_test('fill_IDs')

pg2 = fill_IDs(pg, uniprot_hs)
pg3 = removeDuplicates(pg2)
pg4 = addGeneInfo(pg3, uniprot_hs)
expDesign = createExperimentalDesign(pg4)

dep_res = performDEP(pg_clean, expDesign, thr = 0)


pg_clean = cleanPg(pg, uniprot_hs)
