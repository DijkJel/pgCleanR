% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fill_IDs.R
\name{fill_IDs}
\alias{fill_IDs}
\title{Fill IDs and gene names based on UniProt}
\usage{
fill_IDs(pg, uniprot, species = "Hs")
}
\arguments{
\item{pg}{A maxquant proteinGroups.txt file as data frame with protein IDs in the 'Protein.IDs' column and gene names in the 'Gene.names' column}

\item{uniprot}{A UniProt data frame with all with reviewed and unreviewed human protein sequences. Comes with the package.}

\item{species}{The species to be used for UniProt data, default is 'Hs'}
}
\value{
A proteinGroups data frame with reviewed protein IDs and official gene symbols, where possible.
}
\description{
This function fills in the IDs and gene names of proteins based on UniProt data.
It takes in a proteinGroups data frame with protein IDs in the 'Protein.IDs' column and gene names in the 'Gene.names' column.
It tries to assign a single, unique proteinID and corresponding official gene symbol to each entry of pg.
If a protein ID in the input data frame matches a reviewed protein in UniProt, the corresponding gene name is used.
If the protein ID matches an unreviewed protein in UniProt, the function attempts to find a reviewed protein with the same gene name, or failing that, uses the protein ID as the gene name.
If a protein ID does not match any protein in UniProt, the function uses the protein ID as the gene name.
}
\examples{
pg = pg[1:50,]
fill_IDs(pg, uniprot_hs)

}
