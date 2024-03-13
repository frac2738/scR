
# scR <img src="man/figures/scR_logo.png" align="right" alt="" width="120" />

## Overview

The R package **scR** contains a collection of functions useful for the analysis 
of singel-cell datasets, with particular focus on the integration and the annotation 
with the seurat package. 

## Installation

```
download.file('https://www.dropbox.com/s/up92p4k7hew6vl7/scR_0.9.8.tar.gz?dl=1', destfile = tempfile(), method = "wget")
install.packages(paste0(tempfile(),"/scR_0.9.8.tar.gz"), repos = NULL, type = "source")
library(scR)
```

## Usage

Refer to the **Articles** section for some practical usage examples.

## References

The uses functions or data from the following packages:

<table>
  <tr><td>
     <b>Seurat package:</b>
     <i>Hao et al.</i>
     <i>,Cell 2021</i>
     <a href='https://www.sciencedirect.com/science/article/pii/S0092867421005833?via%3Dihub'>Integrated analysis of multimodal single-cell data</a>
  </td></tr>
  
  <tr><td>
     <b>SoupX package:</b>
     <i>Young, M.D., Behjati, S.</i>
     <i>GigaScience, Issue 12, 2020</i>
     <a href='https://doi.org/10.1093/gigascience/giaa151'>SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data</a>
  </td></tr>
  
  <tr><td>
     <b>Garnett package:</b>
     <i>Hannah A. Pliner, Jay Shendure*, and Cole Trapnell</i>
     <i>,Nature Methods 2019</i>
     <a href='https://cole-trapnell-lab.github.io/garnett/papers/'>Supervised Classification Enables Rapid Annotation of Cell Atlases</a>
  </td></tr>
  
  <tr><td>
     <b>STEMNET package:</b>
     <i>Velten et al.</i>
     <i>,Nature Cell Biology 1019</i>
     <a href='https://www.nature.com/articles/ncb3493'>Human haematopoietic stem cell lineage commitment is a continuous process</a>
  </td></tr>
  
  <tr><td>
     <b>ArchR package:</b>
     <i>Granja et al.</i>
     <i>,Nature Genetics 2021</i>
     <a href='https://www.archrproject.com/index.html'>ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis.</a>
  </td></tr>
  
  <tr><td>
     <b>Ucell package:</b>
     <i>Andreatta et al.</i>
     <i>,Comput Struct Biotechnol J. 2021</i>
     <a href='https://www.sciencedirect.com/science/article/pii/S2001037021002816?via%3Dihub'>Cell: Robust and scalable single-cell gene signature scoring.</a>
  </td></tr>
  
</table>


