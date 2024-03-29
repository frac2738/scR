% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{mouse_human_gene_conversion}
\alias{mouse_human_gene_conversion}
\title{Convert gene names from human to mouse}
\usage{
mouse_human_gene_conversion(genelist, source = "SymbolH", target = "SymbolM")
}
\arguments{
\item{genelist}{a list of human gene symbols}
}
\description{
Convert gene names from human to mouse
}
\examples{
\dontrun{mouse_human_gene_conversion(c("G6PD","CD8A","PPBP"), source = "H")}
}
