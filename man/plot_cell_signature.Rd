% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation_utils.R
\name{plot_cell_signature}
\alias{plot_cell_signature}
\title{Plot gene signatures}
\usage{
plot_cell_signature(
  srtObject,
  figDIR,
  cluster_col,
  srt_assay = "RNA",
  fileDir = "data/",
  colori = c("#FDE725FF", "#440154FF"),
  organism = "human",
  ncores = 8,
  pop = NULL,
  raster = TRUE
)
}
\arguments{
\item{srtObject}{Seurat object.}

\item{figDIR}{output directory}

\item{cluster_col}{Metadata column containing the clustering information.}

\item{fileDir}{path containing the ABIS_markers_mmc4_DEGmodules.csv file}

\item{colori}{colour palettes to use}
}
\description{
This function uses the gene markers defined in the file \code{ABIS_markers_mmc4_DEGmodules.csv}
to generate gene signatures and plot them using both dotplot and featureplot
}
\examples{
\dontrun{
plot_cell_signature(exp_Integrated, ".", "integrated_snn_res.0.8")
}
}
\keyword{markers}
\keyword{seurat}
