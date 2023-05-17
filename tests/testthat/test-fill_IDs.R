test_that("IDs correctly filled", {

  out = fill_IDs(pg[1:50,], uniprot_hs)
  out = out[,c('Protein.IDs', 'Gene.names')]

  expect_type(out, 'list')
  
  expect_snapshot(out)
})


# test_that("IDs correctly filled", {
#   expect_snapshot(fill_IDs(pg[1:10,], uniprot_hs))
# })
