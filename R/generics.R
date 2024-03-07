
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Export a seurat object to file
#'
#' Convenient wrapper to export a seurat object to files
#' @param srtObject seurat object to export
#' @export 
#' @examples
#' \dontrun{
#' export_seurat(pbmc_small)
#' }
export_seurat <- function(srtObject, Assay="RNA", Slot="data", bpcells = FALSE){
   
   suppressPackageStartupMessages(library(cli))
   suppressPackageStartupMessages(library(Matrix))
   
   # create output directory
   outDIR <- "seurat_export"
   if(!exists(outDIR)){
      dir.create(outDIR)
   }
   
   cli_alert_success("* Exporting metadata")
   srtObject@meta.data$barcode <- rownames(srtObject@meta.data)
   srtObject@meta.data$UMAP_1 <- srtObject@reductions$umap@cell.embeddings[,1]
   srtObject@meta.data$UMAP_2 <- srtObject@reductions$umap@cell.embeddings[,2]
   write.csv(srtObject@meta.data, file='seurat_export/metadata.csv', quote=F, row.names=F)

   cli_alert_success("* Exporting counts matrix")
   if(bpcells){
      counts_matrix <- as(object = srtObject[["RNA"]]$counts, Class = "dgCMatrix")
   } else {
      counts_matrix <- GetAssayData(srtObject, assay=Assay, slot=Slot)
   }
   writeMM(counts_matrix, file="seurat_export/counts.mtx")

   cli_alert_success("* Exporting dimensionality reduction matrix")
   write.csv(srtObject@reductions$pca@cell.embeddings, file='seurat_export/pca.csv', quote=F, row.names=F)

   cli_alert_success("* Exporting gene names")
   write.table(
      data.frame('gene'=rownames(counts_matrix)),file='seurat_export/gene_names.csv',
      quote=F,row.names=F,col.names=F
   )
   
   cli_alert_success("* The files were saved in ./seurat_export")
   
}