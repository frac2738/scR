
scType_gene_sets_prepare <- function(path_to_db_file = NULL, cell_types, species = "human"){
   
   library(scR)
   library(Seurat)
   library(HGNChelper)
   library(dplyr)
   
   if(is.null(path_to_db_file)){
      data(ScTypeDB)
      cell_markers <- ScTypeDB
   } else {
      cell_markers <- openxlsx::read.xlsx(path_to_db_file)
   }
   
   cell_markers <- cell_markers[cell_markers$tissueType %in% cell_types,] 
   cell_markers$Expressed = gsub(" ","",cell_markers$Expressed)
   cell_markers$notExpressed = gsub(" ","",cell_markers$notExpressed)
   
   cell_markers$Expressed = sapply(1:nrow(cell_markers), function(i){
      
      markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$Expressed[i],",")))
      markers_all <- markers_all[markers_all != "NA" & markers_all != ""]
      markers_all <- sort(markers_all)
      
      if(length(markers_all) > 0){
         suppressMessages({markers_all <- unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
         paste0(markers_all, collapse=",")
      } else {
         ""
      }
   })
   
   cell_markers$notExpressed = sapply(1:nrow(cell_markers), function(i){
      
      markers_all <- gsub(" ", "", unlist(strsplit(cell_markers$notExpressed[i],",")))
      markers_all <- markers_all[markers_all != "NA" & markers_all != ""]
      markers_all <- sort(markers_all)
      
      if(length(markers_all) > 0){
         suppressMessages({markers_all <- unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
         paste0(markers_all, collapse=",")
      } else {
         ""
      }
   })
   
   cell_markers$Expressed <- gsub("///",",",cell_markers$Expressed)
   cell_markers$Expressed <- gsub(" ","",cell_markers$Expressed)
   
   cell_markers$notExpressed <- gsub("///",",",cell_markers$notExpressed)
   cell_markers$notExpressed <- gsub(" ","",cell_markers$notExpressed)
   
   genes_expressed <- lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$Expressed[j]),","))))
   names(genes_expressed)  <- cell_markers$cellName
   genes_not_expressed <- lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$notExpressed[j]),","))))
   names(genes_not_expressed) <- cell_markers$cellName
   
   return(list(gs_positive = genes_expressed, gs_negative = genes_not_expressed))
   
}

sctype_calculate_score <- function(scRNAseqData, assay = NULL, slot = NULL, seurat5 = FALSE, 
                         scaled = TRUE, gs_positive, gs_negative = NULL, gene_names_to_uppercase = FALSE, ...){
   
   library(scR)
   library(Seurat)
   library(HGNChelper)
   library(dplyr)
   
   # check input matrix
   if(class(scRNAseqData)[1] == "Seurat"){
      message("scRNAseqData is a seurat object")
      scRNAseqData <- eval(parse(text = paste0("scRNAseqData@assays$",assay,"@",slot)))
      
   } else {
      if(sum(dim(scRNAseqData))==0){
         warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
      }
   }
   
   message("* Ranking marker genes by sensitivity")
   marker_stat <- sort(table(unlist(gs_positive)), decreasing = T)
   marker_sensitivity <- data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs_positive),1)),
                                    gene_ = names(marker_stat), stringsAsFactors = !1)
   
   # subselect genes only found in data
   names_positive_cp <- names(gs_positive)
   names_negative_cp <- names(gs_negative);
   
   gs_positive <- lapply(1:length(gs_positive), 
                         function(d_){ GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs_positive[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
   gs_negative <- lapply(1:length(gs_negative), 
                         function(d_){ GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs_negative[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
   names(gs_positive) <- names_positive_cp
   names(gs_negative) = names_negative_cp
   
   cell_markers_genes_score <- marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs_positive)),]
   
   # z-scale if not
   if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
   
   message("* Scoring the cells")
   for(jj in 1:nrow(cell_markers_genes_score)){
      Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
   }
   
   # subselect only with marker genes
   Z = Z[unique(c(unlist(gs_positive),unlist(gs_negative))), ]
   
   message("* Combining the scores")
   es = do.call("rbind", lapply(names(gs_positive), function(gss_){ 
      sapply(1:ncol(Z), function(j) {
         gs_z = Z[gs_positive[[gss_]], j]; gz_2 = Z[gs_negative[[gss_]], j] * -1
         sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
         if(is.na(sum_t2)){
            sum_t2 = 0;
         }
         sum_t1 + sum_t2
      })
   })) 
   
   dimnames(es) = list(names(gs_positive), colnames(Z))
   es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
   
   return(es.max)
   
}

scType_assign_score <- function(srtObject, es.max, resolution = "integrated_snn_res.1.2"){
   
   resolution_column <- eval(parse(text = paste0("srtObject@meta.data$",resolution)))
   
   message("* Assigning scores")
   cL_resutls <- do.call("rbind", lapply(unique(resolution_column), 
                                         function(cl){
                                            es.max.cl <- sort(rowSums(es.max[ ,rownames(srtObject@meta.data[resolution_column == cl, ])]), decreasing = TRUE)
                                            head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(resolution_column==cl)), 10)
                                         })
   )
   
   write.csv(cL_resutls,"./scType_cL_resutls.csv")
   
   sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
   
   message("* Setting low confident clusters to 'Unknown' ")
   sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
   
   message("* Saving")
   srtObject@meta.data$ScType = ""
   for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]; 
      srtObject@meta.data$ScType[resolution_column == j] = as.character(cl_type$type[1])
   }
   
   return(srtObject)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' scType
#' 
#' Fully-automated and ultra-fast cell-type identification using specific marker 
#' combinations from single-cell transcriptomic data.
#' The functions were modified and adapted from the original repository
#' (https://github.com/IanevskiAleksandr/sc-type/tree/master).
#' 
#' * Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#' * Adapted by Francesco Carbone <francesco.carbone@institutiamgine.org>, July 2023
#' 
#' @param path_to_db_file list of markers genes. If NULL the internal object will be used.
#' @param cell_types Tissue type. Used to select the markers.
#' @param srtObject input data, either a matrix or a seurat object
#' @param assay Assay from seurat
#' @param slot Slot from seurat
#' @param scaled  TRUE if input data is already scaled
#' @param resolution  clustering resolution to pass to scType
#' @keywords annotation
#' @export
#' @examples
#' \dontrun{scType_pipeline()}
#' 
scType_pipeline <- function(path_to_db_file = NULL, cell_types, srtObject, assay = "RNA", slot = "scale.data",
                            scaled = scaled, resolution="integrated_snn_res.1.2"){
   
   gs_list <- scType_gene_sets_prepare(path_to_db_file, cell_types)
   
   es.max <- sctype_calculate_score(srtObject, assay=assay, slot=slot,
                          scaled = scaled,  gs_positive = gs_list$gs_positive, gs_negative = gs_list$gs_negative) 
   
   srtObject <- scType_assign_score(srtObject, es.max, resolution=resolution)
   
   return(srtObject)
}
   
   
   
   
   
   
