% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_pipelines.R
\name{run_qc_and_filtering}
\alias{run_qc_and_filtering}
\title{Check sample quality and filter the count matrix}
\usage{
run_qc_and_filtering(
  exp_meta,
  sample_num,
  pipeline,
  minFeatures,
  pcDims = 30,
  mtPerc = 20,
  source = "raw",
  source_type = "dir",
  organism = "human",
  isVerbose = FALSE,
  vst_flavor = NULL,
  input_matrix = NULL,
  isFixed = FALSE,
  remove_doublets = FALSE,
  mtPattern = NULL,
  riboPattern = NULL,
  seurat5 = FALSE,
  bpcells_dir = NULL,
  labeling_celline = FALSE
)
}
\arguments{
\item{exp_meta}{Metadata table (\code{data.frame} by default)}

\item{sample_num}{Row number of the metadata file containing the sample information}

\item{pipeline}{Perform SCT or lognorm normalisation}

\item{minFeatures}{Maximum number of genes to filter out}

\item{pcDims}{Number of PCs (30 by default)}

\item{source}{Use the data from the \code{raw} or \code{filtered} cellranger matrix (\code{raw} by default)}

\item{source_type}{\code{dir} by default}

\item{organism}{Set the sample organism (\code{human} by default)}

\item{isVerbose}{Print all the messages on screen (\code{FALSE} by default)}

\item{vst_flavor}{sctransform version}
}
\value{
Returns 
\item{}{a folder, named as the sample ID, containing QC figures and tables}
\item{}{a seurat object with filtering and clustering applied}
}
\description{
This function is used to reads the output from cellranger, create a seurat object, 
apply basic filtering, and normalise and cluster the filtered counts.
}
\examples{
\dontrun{
run_qc_and_filtering(exp_meta,number,".","SCT", 500, 30, source = "raw", source_type = "outs", isVerbose=FALSE)
}

}
\keyword{seurat}
