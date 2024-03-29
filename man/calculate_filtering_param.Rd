% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_utils.R
\name{calculate_filtering_param}
\alias{calculate_filtering_param}
\title{Estimate the filtering parameters (MT, nFeature, nCount)}
\usage{
calculate_filtering_param(sampleID, srtObject, outDIR, apply = FALSE)
}
\arguments{
\item{sampleID}{SampleID}

\item{srtObject}{Unfilterd seurat object}

\item{outDIR}{Save figures}

\item{apply}{Apply the calculated filters | Default \code{FALSE}}
}
\description{
Estimate the filtering parameters (MT, nFeature, nCount)
}
\details{
Details

The thresolds are set using the MAD (Median absolute deviation), and defined as:
\itemize{
\item MAX = median(log10(metric)) + 3*mad(log10(metric))
\item MIN = median(log10(metric)) - 3*mad(log10(metric))
}

Refer to the pre-processing article for a detailed usage and examples
}
\examples{
\dontrun{
filt_st <- calculate_filtering_param("45L2", unfilt_st, ".", apply = TRUE)
}
}
\keyword{soupX}
