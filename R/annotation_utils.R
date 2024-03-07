#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Validate markers
#' 
#' @param object Seurat object.
#' @param marker_file_path output directory
#' @param noFiltering FALSE
#' @keywords markers
#' @export
#' @examples
#' \dontrun{
#' validate_markers()
#' }
validate_markers <- function(object,marker_file_path, organism = "human"){
   
   library(monocle3)
   library(SeuratWrappers)
   library(garnett)
   
   if(organism == "human"){
      library(org.Hs.eg.db)
      db_to_use=org.Hs.eg.db
      
   } else {
      library(org.Mm.eg.db)
      db_to_use=org.Mm.eg.db
      
   }
   
   message("* Extracting the raw counts ")
   rawcounts <- as(as.matrix(object@assays$RNA@counts),"sparseMatrix")
   fData <- data.frame(gene_short_name = row.names(rawcounts), row.names = row.names(rawcounts))
   
   tmp_cds <- new_cell_data_set(rawcounts,
                                 cell_metadata = as.data.frame(object@meta.data),
                                 gene_metadata = fData)
   message("  - Normalising")
   tmp_cds <- preprocess_cds(tmp_cds, num_dim = 30,verbose=FALSE)
   message("  - Batch correction")
   tmp_cds <- align_cds(tmp_cds, alignment_group = c("SampleID"),verbose = FALSE) # remove batch
   
   message("  - Check markers")
   marker_check <- check_markers(tmp_cds, marker_file_path,
                                 db=db_to_use,
                                 cds_gene_id_type = "SYMBOL",
                                 marker_file_gene_id_type = "SYMBOL")
   
   drop <- which(marker_check$summary == "Not in db")
   
   if(length(drop)>0){
      marker_check <- marker_check[-drop]
         
   }
   
   return(marker_check)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Parser
#' 
#' @param file_str File
#' @param debug FALSE
#' @keywords markers
#' @export
#' @examples
#' \dontrun{
#' parse_input()
#' }
parse_input <- function(file_str, debug = F) {
   
   library(garnett)
   
   source("R/parser_cell_type_markers.R")
   # Parse input_file
   lexer  <- rly::lex(Lexer, debug=debug)
   parser <- rly::yacc(Parser, debug=debug)
   parse_list <- parser$parse(file_str, lexer)
   
   parse_list
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Another parser
#' 
#' @param marker_file File
#' @param useDefault FALSE
#' @keywords markers
#' @export
#' @examples
#' \dontrun{
#' parse_markers_lists()
#' }
parse_markers_lists <- function(marker_file, useDefault=TRUE){
   
   if(useDefault){
      #marker_file <- "data/00_pbmc_cell_markers.txt"
      data("ABIS_markers")
      marker_file <- ABIS_markers
   }
   
   #file_str <- paste0(readChar(marker_file, file.info(marker_file)$size),"\n")
   #parse_list <- parse_input(file_str)
   #orig_name_order <- unlist(parse_list[["name_order"]])
   #rm("name_order", envir=parse_list)
   parse_list <- marker_file
   
   return(parse_list)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Yes, another parser
#' 
#' @param object File
#' @param Markers FALSE
#' @param areMarkers FALSE
#' @param listCells FALSE
#' @param markers_file_path NULL
#' @param scale_data FALSE
#' @param cluster FALSE
#' @param cut 1
#' @param ncolumns NULL
#' @param colori c("#FAC127FF","#280B54FF")
#' @keywords markers
#' @export
#' @examples
#' \dontrun{
#' markers_lookup()
#' }
markers_lookup <- function(object , Markers,
                           areMarkers = FALSE, listCells = FALSE,
                           markers_file_path = NULL, scale_data = TRUE, cluster = FALSE, cut = 1, 
                           ncolumns = NULL, colori = c("#FDE725FF","#440154FF")){
   
   library(ggplot2)
   library(patchwork)
   library(scCustomize)
   
   myColScale <- c("#0000FF","grey","#FF4C4C")
   names(myColScale) <- c("Low","Middle","High")
   
   # load markers from file
   if(!is.null(markers_file_path)){
      markers_file <- parse_markers_lists(markers_file_path, useDefault = FALSE)
   } else {
      #markers_file <- parse_markers_lists()
      #data("pbmc_markers_annotated")
      #markers_file <- pbmc_markers_annotated
   }

   if(listCells){
      # fix message to account for user provided markers file
      message(paste0("The file 00_pbmc_cell_markers.txt contains curated markers for the 
                     following cell types:"))
      message(" ")
      names(markers_file)
   } else {
      
      if(!areMarkers){
         Markers <- eval(parse(text = paste0("markers_file$`",Markers,"`")))
         Markers <- Markers@expressed
      }
      
      Markers <- unique(Markers)
      Markers <- Markers[order(Markers,decreasing = FALSE)]
      
      message(paste0("Plotting ",paste0(Markers,collapse = " ")))
      message("  ")
      feat_plot <- FeaturePlot(object, features = Markers, 
                               pt.size = 1, order = TRUE, 
                               cols = colori, ncol = ncolumns)
      
      if(!cluster){
         dot_plot <- DotPlot(object,features = Markers, scale=scale_data) + scale_color_gradientn(colors = colori)
      } else {
         dot_plot <- Clustered_DotPlot(object, features = Markers, k = cut)
            DotPlot(object,features = Markers, scale=scale_data) + scale_color_gradientn(colors = colori)
         
      }
      
      myList <- list()
      myList[["feat_plot"]] <- feat_plot
      myList[["dot_plot"]] <- dot_plot
      
      return(myList)
   }
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Plot cell markers
#'
#' This function uses the gene markers defined in the file \code{pbmc_annotation_markers.csv}
#' to plot and save their expression in both dotplot and featureplot.
#' 
#' 
#' @param srtObject Seurat object.
#' @param figDIR output directory
#' @param cluster_col Metadata column containing the clustering information.
#' @param overlapPlot something
#' @keywords seurat markers
#' @export
#' @examples
#' \dontrun{
#' plot_cells_markers(exp_Integrated, figDIR, "integrated_snn_res.0.8", overlapPlot=NULL)
#' }
plot_cells_markers <- function(srtObject, figDIR, cluster_col, srt_assay = "RNA", overlapPlot=NULL){
   
   DefaultAssay(srtObject) <- srt_assay
   Idents(srtObject) <- cluster_col
   
   data("pbmc_markers_annotated")
   pbmc_markers <- pbmc_markers_annotated
   
   names(pbmc_markers) <- gsub("\\.","_",names(pbmc_markers))
   
   for(cellType in names(pbmc_markers)){
      
      message(paste0("** Processing ",cellType))
      
      markers_genes <- unique(pbmc_markers[,cellType])
      drop <- which(markers_genes == "")
      
      if(length(drop)>0){
         markers_genes <- markers_genes[-drop]
      }
      
      tmp_markers <- markers_lookup(srtObject,markers_genes,areMarkers=TRUE, scale_data = TRUE)
      dotplot_overlap <- if (is.null(overlapPlot)) tmp_markers[[2]] else figures_overlap(tmp_markers[[2]],ref_cluster)
      
      out_file_name <- paste0("markers_",cellType)
      SaveFigure(tmp_markers[[1]], figDIR,out_file_name, type = "png", 16, 9, 300, ggplot = TRUE )
      SaveFigure(dotplot_overlap, figDIR, paste0(out_file_name,"_dot"), type = "png", 16, 9, 300, ggplot = TRUE )
      
   }
   
}

plot_cells_markers_alphab <- function(srtObject, figDIR, srt_assay = "RNA"){
   
   library(SCpubr)
   
   DefaultAssay(srtObject) <- srt_assay

   data("pbmc_markers_annotated")
   pbmc_markers <- unique(unlist(pbmc_markers_annotated))
   
   dotplot_plot <- do_DotPlot(srtObject, pbmc_markers, use_viridis = TRUE, viridis.palette = "D", viridis.direction = -1, scale = TRUE, cluster = TRUE)
   SaveFigure(dotplot_plot, figDIR, paste0("all_genes_dotplot"), type = "png", 22, 12, 300, ggplot = TRUE )
   
   for(gene in pbmc_markers){
      tryCatch({
         
         feat_plot <- do_FeaturePlot(srtObject, features = gene, use_viridis = TRUE, viridis.palette = "D", viridis.direction = -1,order = TRUE)
         dot_plot <- do_DotPlot(srtObject, features = gene, use_viridis = TRUE, viridis.palette = "D", viridis.direction = -1, scale = TRUE, cluster = TRUE)
         nebula_plot <- do_NebulosaPlot(srtObject, features = gene)
         
         final_plot <- feat_plot | dot_plot|nebula_plot
         SaveFigure(final_plot, figDIR, paste0(gene), type = "png", 22, 12, 300, ggplot = TRUE )
         
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
   }
   
}


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Plot gene signatures
#'
#' This function uses the gene markers defined in the file \code{ABIS_markers_mmc4_DEGmodules.csv}
#' to generate gene signatures and plot them using both dotplot and featureplot
#' 
#' 
#' @param srtObject Seurat object.
#' @param figDIR output directory
#' @param cluster_col Metadata column containing the clustering information.
#' @param fileDir path containing the ABIS_markers_mmc4_DEGmodules.csv file
#' @param colori colour palettes to use
#' @keywords seurat markers
#' @export
#' @examples
#' \dontrun{
#' plot_cell_signature(exp_Integrated, ".", "integrated_snn_res.0.8")
#' }
plot_cell_signature <- function(srtObject ,
                                figDIR,
                                cluster_col,
                                srt_assay = "RNA",
                                fileDir= "data/",
                                colori = c("#FDE725FF","#440154FF"),
                                organism = "human",
                                ncores = 8,
                                pop = NULL,
                                raster = TRUE
                                ){
   
   library(Seurat)
   library(UCell)
   library(wesanderson)
   
   DefaultAssay(srtObject) <- srt_assay
   
   data("ABIS_markers")
   markers <- ABIS_markers
   
   names(markers) <- c("GeneName","cellType","FDR","module")
   
   drop <- which(markers$GeneName == "") 
   markers <- markers[-drop,]
   
   markers$cellType <- gsub(pattern = " ", replacement = "_", markers$cellType)
   markers$cellType <- gsub(pattern = "\\+", replacement = "-", markers$cellType)
   markers$cellType <- gsub(pattern = "/", replacement = "-", markers$cellType)
   markers$cellType <- factor(markers$cellType)
   
   dir.create(path <- paste0(figDIR,"/cellType_signatures/"))
   
   if(organism == "mouse"){
      markers$GeneName <- gene_names_H2M(markers$GeneName)
      
   }
   
   markers <- markers[,c(1,2)]
   signatures2test <- split(markers, f = markers$cellType)
   signatures2test <- lapply(signatures2test, function(x) { x["cellType"] <- NULL; x })
   signatures2test <- lapply(signatures2test, unlist)
   
   if(!is.null(pop)){
      
      signatures2test <- signatures2test[pop]
      
   }
   
   message("** Calculating the scores...")
   srtObject <- AddModuleScore_UCell(srtObject, 
                                     features=signatures2test, name = NULL,
                                     ncores = ncores)
   
   message("** Plotting the scores...")
   for(cells in names(signatures2test)){
   
      tryCatch({
         message(paste0("  * ",cells))
         
         sigPlot <- do_FeaturePlot(srtObject, features = cells, use_viridis = TRUE, viridis.palette = "D", viridis.direction = -1,order = TRUE)
         dot_plot <- do_DotPlot(srtObject, features = cells, use_viridis = TRUE, viridis.palette = "D", viridis.direction = -1, scale = TRUE, cluster = TRUE)
         nebula_plot <- do_NebulosaPlot(srtObject, features = cells)
         
         final_plot <- sigPlot|dot_plot|nebula_plot
         ggsave(plot = final_plot, filename = paste0(figDIR,"/cellType_signatures/",cells,"_signature.png"),bg="white",width = 16,height = 9,dpi=300)
         
         # sigPlot <- FeaturePlot(srtObject, cells, order = TRUE, cols = colori, raster = raster)
         # ggsave(plot = sigPlot, filename = paste0(figDIR,"/cellType_signatures/",cells,"_signature.png"),bg="white",width = 16,height = 9,dpi=300)
         # 
         # Idents(srtObject) <- cluster_col
         # sigPlot_dot <- DotPlot(srtObject, features = cells) + scale_color_gradientn(colors = colori)
         # ggsave(plot = sigPlot_dot, filename = paste0(figDIR,"/cellType_signatures/",cells,"_signature_dot.png"),bg="white",width = 6,height = 9,dpi=300)
         
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
   }
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Calcualte genes signature
#'
#' This function 
#' 
#' @param srtObject Seurat object.
#' @keywords seurat markers
#' @export
#' @examples
#' \dontrun{
#' calcualte_genes_signature(exp_Integrated, ".", "integrated_snn_res.0.8")
#' }
calculate_genes_signature <- function(srtObject, markers_list, 
                                assay = "RNA",
                                ncores = 2
){
   
   library(Seurat)
   library(UCell)

   DefaultAssay(srtObject) <- assay
   
   message("** Calculating the scores...")
   srtObject <- AddModuleScore_UCell(srtObject, 
                                     features= markers_list, name = NULL,
                                     ncores = ncores)
   
   srtObject <- SmoothKNN(srtObject, signature.names = names(markers_list), reduction="pca")
   
   return(srtObject)
   
}
