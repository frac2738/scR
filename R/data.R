#' PBMC markers
#'
#' List of markers used to annotate PBMC cell type
#' @format garnet marker file
#' @examples
#' data(pbmc_markers)
"pbmc_markers"

#' Gene markers of 29 immune cell types
#'
#' List of immune markers from Monaco et al.
#'
#' @format dataframe
#' @references Monaco et al. (2019) Cell reports vol. 26,6
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30726743/}{PubMed})
#' @source \href{https://github.com/giannimonaco/ABIS}{ABIS}
#' @examples
#' data(ABIS_markers)
#' head(ABIS_markers)
#' unique(ABIS_markers$Cell.type)
"ABIS_markers"

#' PBMC markers as data frame
#'
#' PBMC markers
#' @format dataframe
#' @examples
#' data(pbmc_markers_annotated)
#' names(pbmc_markers_annotated)
#' unique(pbmc_markers_annotated$T)
"pbmc_markers_annotated"

#' PBMC 1K test dataset
#'
#' Small pbmc dataset with around 1500 cells.
#' 
#' @format seurat
#' @examples
#' data(pbmc_1k)
#' DimPlot(pbmc_1k, label = TRUE, repel = TRUE) + NoLegend()
"pbmc_1k"

#' Mouse-Human gene name mapping table 
#'
#' Table with gene names in human and mouse format.
#' 
#' @format txt
#' @examples
#' data(MouseHumanGenes)
#' head(MouseHumanGenes)
"MouseHumanGenes"

#' ScTypeDB
#'
#' Table with marker genes
#' 
#' @format txt
#' @examples
#' data(ScTypeDB)
#' head(ScTypeDB)
"ScTypeDB"
