
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Setup project directory structure
#'
#' Create reproducible project directory organization when initiating a new analysis.
#'
#' @param project_location Directory in which the project will be located.  Default is "getwd()"
#' @param project_name Name of the project
#' @import cli
#' @importFrom data.table fread
#'
#' @export
#'
#' @return no return value.  Creates system folders.
#'
#' @concept organization_util
#'
#' @examples
#' \dontrun{
#' # If using built-in directory structure.
#' 
#' project_dir <- "/home/projects"
#' projectID <- "new_project
#' 
#' setup_sc_project(project_dir, projectID)
#' }
#' This will create a folder with the following structure
#'
#' new_project
#' |--- 01_metadata 
#' |--- 02_raw_data
#' |--- 03_processed_data
#'   |--- cellranger_summaries
#'   |--- cloupe_files
#'   |--- Samples_QC
#' |--- 04_additional_data
#' |--- 05_scripts
#' |--- 06_reports
#' |--- analysis
#'
setup_sc_project <- function(project_location = getwd(), project_name = "new_project") {
   
   # File paths setup
   project_path <- paste0(project_location,"/",project_name)
   
   output_dirs <- list(
      project_home <- project_path,
      metadata_path = paste0(project_path,"/01_metadata/"),
      data_path = paste0(project_path,"/02_raw_data/"),
      object_path = paste0(project_path,"/03_processed_data/"),
      summaries_path = paste0(project_path,"/03_processed_data/cellranger_summaries/"),
      cloupe_path = paste0(project_path,"/03_processed_data/cloupe_files/"),
      plots_qc_path = paste0(project_path,"/03_processed_data/Samples_QC/"),
      metadata_path = paste0(project_path,"/04_additional_data/"),
      scripts_path = paste0(project_path,"/05_scripts/"),
      reports_path = paste0(project_path,"/06_reports/"),
      analysis_path = paste0(project_path,"/analysis/"))
   
   # Check for directories and create new ones
   lapply(output_dirs, function(dir_path){
      if (!dir.exists(dir_path)){
         dir.create(path = dir_path)
         
      } else {
         cli_warn(message = "The directory {.val {dir_path}} aleady exists.  No new directory created.")
      
      }
   })
   
   # Print completion message
   cli_inform(message = "Setup complete")
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Wrapper function to save figures to disk.
#'
#' Convenient wrapper to save figures with either ggplot (default) or base R.
#' @param plot_to_save Figure to save
#' @param fig_path Where to save the figure
#' @param output_name Output file name
#' @param type Save as png or pdf (png by default)
#' @param width Set the figure width
#' @param height Set the figure height
#' @param res dpi (300 by default)
#' @param ggplot Save the figure using ggsave (TRUE by default)
#' @export 
#' @examples
#' \dontrun{
#' aPlot <- ggplot(...)
#' SaveFigure(aPlot, "./figures","figure_1","png",1400,900,300)
#' }
SaveFigure <- function(plot_to_save, fig_path = ".", output_name, type = "png", width, height, res = 300, ggplot = TRUE){
   
   if(ggplot){
      ggsave(plot_to_save, filename = paste0(fig_path,"/",output_name,".",type),
             bg="white", dpi = res, width = width, height = height)
      
   } else {
      if(type == "png") {
         png(paste0(fig_path,"/",output_name, ".", type),
             width = width, height = height, units = "in", res = 300)
         
      } else {
         pdf(paste0(fig_path,"/",output_name, ".", type), width = width, height = height)
      }
      print(plot_to_save)
      dev.off()
   }
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Function that, given a DESeq2 object, runs PCA and it returns the contributions of the  
#' genes to each PCs.
#' 
#' @param vst_object vst normalised object (e.g. DESeq2)  
#' @param ntop Number of genes to use
#' @param pc_sel Which PC to use for ordering the result table
#' @param returnPCA Return the pca object (FALSE by default)
#' @export 
#' @examples
#' \dontrun{PCA_genes(vst)}
PCA_genes <- function(vst_object, ntop= 500, pc_sel = 1, returnPCA = FALSE){
   
   rv <- rowVars(assay(vst_object))
   
   # select the ntop genes by variance
   select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
   
   # perform a PCA on the data in assay(x) for the selected genes
   pca <- prcomp(t(assay(vst_object)[select,]))
   
   # the contribution to the total variance for each component
   percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
   names(percentVar) <- paste0("PC",seq(1:length(percentVar)))
   
   pca_rotations <- pca$rotation
   pca_rotations <- pca_rotations[order(pca_rotations[,pc_sel], decreasing = TRUE),]
   
   if(returnPCA){
      return(pca)
   } else {
      return(as.data.frame(pca_rotations))
   }

}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Convert gene names from human to mouse
#' 
#' @param genelist a list of human gene symbols
#' @export 
#' @examples
#' \dontrun{mouse_human_gene_conversion(c("G6PD","CD8A","PPBP"), source = "H")}
mouse_human_gene_conversion <- function(genelist, source = "SymbolH", target = "SymbolM"){
   
   data(MouseHumanGenes)
   
   accepted_values <- c("SymbolH","SymbolM","EnsemblH")
   converted_genes <- c()
   
   if(source %in% accepted_values && source %in% accepted_values){
      
      genelist_idx <- which(MouseHumanGenes[[source]] %in% genelist)
      mapped_genes <- MouseHumanGenes[genelist_idx,]
      
      converted_genes <- mapped_genes[[target]]
      
   } else {
      message("Source/Target values not found. Accepted values are: SymbolM, SymbolH and EnsemblH")
   }

   # genelist <- tolower(genelist)
   # genelist <- paste0(toupper(substr(genelist, 1, 1)), substr(genelist, 2, nchar(genelist)))
   
   return(converted_genes)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Convert the rownames of a matrix from ensemble to symbol
#' 
#' @param genelist a matrix with rownames in ensemle format
#' @export 
#' @examples
#' \dontrun{convert_ensemble_to_symbol(matrice)}
convert_ensemble_to_symbol <- function(matrice){
   
   name_df <- data.frame(orig.id = c(rownames(matrice)))
   name_df$EnsemblH <- name_df$orig.id
   name_df$EnsemblH <- as.character(name_df$EnsemblH)
   name_df$EnsemblH <- sub("[.][0-9]*", "", name_df$EnsemblH)
   
   data(MouseHumanGenes)
   gene.df <- merge.data.frame(name_df,MouseHumanGenes,by = "EnsemblH",all.x = TRUE)
   rownames(gene.df) <- make.unique(gene.df$EnsemblH)
   gene.df <- gene.df[rownames(matrice), ]
   
   gene.df <- gene.df[!is.na(gene.df$SymbolH),]
   if(length(which(gene.df$SymbolH == ""))>1){
      drop <- which(gene.df$SymbolH == "")
      gene.df <- gene.df[-drop,]
   }
   
   matrice_filtered <- matrice[gene.df$EnsemblH, ]
   rownames(matrice_filtered) <- make.unique(gene.df[, "SymbolH"])
   return(matrice_filtered)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Calcualte the middle value
#' 
#' @param min min value
#' @param max max value
#' @export 
#' @examples
#' \dontrun{Middle_Number(min = 10, max = 20)}
Middle_Number <- function(min, max) {
   min_max <- c(min, max)
   middle <- min_max[-length(min_max)] + diff(min_max) / 2
   return(middle)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# functions to allow axis duplication with discrete values

guide_axis_label_trans <- function(label_trans = identity, ...) {
   axis_guide <- guide_axis(...)
   axis_guide$label_trans <- rlang::as_function(label_trans)
   class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
   axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
   trained <- NextMethod()
   trained$key$.label <- x$label_trans(trained$key$.label)
   trained
}

