
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Calculate the percentage of mitochondrial and ribosomal genes
#'
#' @param seurat_object seurat object
#' @param organism Human by default.
#' @param new_mito mitochondrial pattern to use if the organisms is not human or mouse
#' @param new_ribo ribosomial pattern to use if the organisms is not human or mouse
#' @details Details
#' Refer to the pre-processing article for a detailed usage and examples
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{
#' calculate_mito_ribo(seurat_object)
#' }
#' @md
calculate_mito_ribo <- function(seurat_object, organism = "human", new_mito = NULL, new_ribo = NULL){
   
   if(is.null(new_mito)|is.null(new_ribo)){
      if(organism == "human"){
         mito_pattern <- "^MT-"
         ribo_pattern <- "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"
      } else if(organism == "mouse"){
         mito_pattern <- "^Mt-|^mt-"
         ribo_pattern <- "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa"
      } else {
         
         cli_abort(message = c("There are no default values for the selected organism.",
                               "Please provide the matching patterns using {.code new_mito} and {.code new_ribo}")
         )
         
      }
   } else {
      mito_pattern <- new_mito
      ribo_pattern <-new_ribo
   }
   
   seurat_object <- PercentageFeatureSet(seurat_object, pattern = mito_pattern, col.name = "percent.mt")
   seurat_object <- PercentageFeatureSet(seurat_object, pattern = ribo_pattern, col.name = "percent.ribo")
   
   return(seurat_object)
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Calculate the complexity of the dataset
#'
#' @param seurat_object seurat object
#' @details Details
#' Refer to the pre-processing article for a detailed usage and examples
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{
#' calcualte_cell_complexity(seurat_object)
#' }
#' @md
calculate_cell_complexity <- function(seurat_object){
   seurat_object@meta.data$log10GenesPerUMI <- log10(seurat_object@meta.data$nFeature_RNA) / log10(seurat_object@meta.data$nCount_RNA)
   
   return(seurat_object)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Compare 2 seurat object (used to generate QC figures between pre-fitlered and filtered data)
#'
#' @param object1 seurat object 1
#' @param object2 seurat object 2
#' @param comp_name name to use for the figures
#' @param colori set of 2 colors to use for plotting | Default \code{c("#D3D3D3","#FFB140")}
#' @param qc_dir output directory
#' @param minFeatures number of genes
#' @param mtPerc percentage fo mt genes
#' @details Details
#' Refer to the pre-processing article for a detailed usage and examples
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{
#' compare_qc(object1, object2, "comparison_1")
#' }
#' @md
compare_qc <- function(object1, object2, comp_name, colori = NULL, qc_dir = ".", minFeatures = 500, mtPerc = 20){
   
   library(Seurat)
   library(ggplot2)
   library(patchwork)
   library(viridis)
   library(cli)
   
   if(is.null(colori)){
      color_pallet <- c("#D3D3D3","#FFB140")
      col_fill_1 <- color_pallet[1]
      col_fill_2 <- color_pallet[2]
      
   }
   
   cli_inform(message = paste0("Calculating metrics for object 1"))
   
   plot1 <- FeatureScatter(object1, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
      geom_hline(yintercept =mtPerc, linetype="dashed")  + scale_color_manual(values = col_fill_1) 
   plot2 <- FeatureScatter(object1, feature1 = "nFeature_RNA", feature2 = "percent.mt") + 
      geom_vline(xintercept = minFeatures, linetype="dashed") +
      geom_hline(yintercept = mtPerc, linetype="dashed") + scale_color_manual(values = col_fill_1) 
   plot3 <- FeatureScatter(object1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + scale_color_manual(values = col_fill_1) 
   
   nGene_1 <- VlnPlot(object1, features = c("nFeature_RNA")) + NoLegend() + scale_fill_manual(values = col_fill_1) +
   ggtitle(paste0("nFeature [median: ",round(median(object1@meta.data$nFeature_RNA))," genes]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   nCount_1 <- VlnPlot(object1, features = c("nCount_RNA")) + NoLegend() +
       scale_fill_manual(values = col_fill_1) +
   ggtitle(paste0("nCount [median: ",round(median(object1@meta.data$nCount_RNA))," counts]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   pMit_1 <- VlnPlot(object1, features = c("percent.mt")) + NoLegend() +
       scale_fill_manual(values = col_fill_1) +
   ggtitle(paste0("% Mitochondrial genes [median: ",round(median(object1@meta.data$percent.mt)),"%]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   pRib_1 <- VlnPlot(object1, features = c("percent.ribo")) + NoLegend() +
      scale_fill_manual(values = col_fill_1) +
   ggtitle(paste0("% Ribosomal genes [median: ",round(median(object1@meta.data$percent.ribo)),"%]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   complexity_1 <- ggplot(object1@meta.data, aes(x=log10GenesPerUMI)) +
      geom_density(alpha = 1, aes(color = SampleID, fill = SampleID)) + scale_fill_manual(values = col_fill_1) +
      xlim(0,1) + 
      scale_colour_manual(values="black") + 
      theme(legend.position = "none", plot.title = element_text(size = 18, hjust = 0.5)) + 
      geom_vline(xintercept = 0.8,linetype="dotted") + 
      ggtitle("Complexity/Novelty score")
   
   plots_1 <- (nGene_1 / nCount_1 / pMit_1 / pRib_1 / complexity_1) + plot_annotation(title = "Pre-filtered")
   
   cli_inform(message = paste0("Calculating metrics for object 2"))
   
   plot4 <- FeatureScatter(object2, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
      geom_hline(yintercept =mtPerc, linetype="dashed") + scale_color_manual(values = col_fill_2) 
   plot5 <- FeatureScatter(object2, feature1 = "nFeature_RNA", feature2 = "percent.mt") + 
      geom_vline(xintercept = minFeatures, linetype="dashed") +
      geom_hline(yintercept = mtPerc, linetype="dashed")  + scale_color_manual(values = col_fill_2) 
   plot6 <- FeatureScatter(object2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  + scale_color_manual(values = col_fill_2) 
   
   nGene_2 <- VlnPlot(object2, features = c("nFeature_RNA")) + NoLegend() + scale_fill_manual(values = col_fill_2) +
   ggtitle(paste0("nFeature [median: ",round(median(object2@meta.data$nFeature_RNA))," genes]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   nCount_2 <- VlnPlot(object2, features = c("nCount_RNA")) + NoLegend() +
      scale_fill_manual(values = col_fill_2) +
   ggtitle(paste0("nCount [median: ",round(median(object2@meta.data$nCount_RNA))," counts]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   pMit_2 <- VlnPlot(object2, features = c("percent.mt")) + NoLegend() +
      scale_fill_manual(values = col_fill_2) +
   ggtitle(paste0("% Mitochondrial genes [median: ",round(median(object2@meta.data$percent.mt)),"%]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   pRib_2 <- VlnPlot(object2, features = c("percent.ribo")) + NoLegend() +
       scale_fill_manual(values = col_fill_2) +
   ggtitle(paste0("% Ribosomal genes [median: ",round(median(object2@meta.data$percent.ribo)),"%]")) +
      theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0))
   
   complexity_2 <- ggplot(object2@meta.data, aes(x=log10GenesPerUMI)) +
      geom_density(alpha = 1, aes(color = SampleID, fill = SampleID)) + scale_fill_manual(values = col_fill_2) +
      xlim(0,1) + 
      scale_colour_manual(values="black") + 
      theme(legend.position = "none", plot.title = element_text(size = 18, hjust = 0.5)) + 
      geom_vline(xintercept = 0.8,linetype="dotted") + 
      ggtitle("Complexity/Novelty score")
   
   plots_2 <- (nGene_2 / nCount_2 / pMit_2 / pRib_2 / complexity_2) + plot_annotation(title = "Filtered")
   
   cli_inform(message = paste0("Combining and saving"))
   
   combined <- plots_1 | plots_2
   ggsave(combined, filename = paste0(qc_dir,"/",comp_name,".png"), dpi = 300, width = 14, height = 18)
   
   ggsave( (plot1/plot4), filename = paste0(qc_dir,"/",comp_name,"_counts-vs-MT.png"), bg = "white", width = 14, height = 12)
   ggsave( (plot2/plot5), filename = paste0(qc_dir,"/",comp_name,"_features-vs-MT.png"), bg = "white", width = 14, height = 12)
   ggsave( (plot3/plot6), filename = paste0(qc_dir,"/",comp_name,"_count-vs-features.png"), bg = "white", width = 14, height = 12)
   

}


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Estimate ambient RNA contamination using SoupX
#'
#' @param sampleid Sample ID
#' @param raw_dataDIR Directory containing the cellranger outputs (e.g. /outs)
#' @param estGenes List of genes to use to estimate the contamination | Default \code{NULL}
#' @param counts2disk Save the corrected counts as h5 | Default \code{FALSE}
#' @details Details
#' Refer to the pre-processing article for a detailed usage and examples
#' @keywords soupX
#' @export
#' @examples
#' \dontrun{
#' estimate_ambient_rna("Sample1",".", estGenes = c("PPBP","HBB","HBA1"), FALSE)
#' }
#' @md
estimate_ambient_rna <- function(sampleid, raw_dataDIR, estGenes = NULL, counts2disk = FALSE){
   
   library(SoupX)
   library(DropletUtils)
   
   message(paste0("Processing ",sampleid))
   
   sc <- load10X(paste0(raw_dataDIR,"/",sampleid))
   sc$metaData$clusters <- as.character(sc$metaData$clusters)
   
   if(is.null(estGenes)){
      sc <- autoEstCont(sc)
      
   } else {
      useToEst <- estimateNonExpressingCells(sc, nonExpressedGeneList = list(EST = estGenes))
      
      if(length(unique(useToEst))>1){
         sc <- calculateContaminationFraction(sc, list(EST = estGenes), useToEst = useToEst)
      } else {
         sc <- autoEstCont(sc)
      }
      
   }
   
   corrected_mtx <- adjustCounts(sc)
   
   estDF <- data.frame(SampleID = sampleid,
                       Contamination = sc$metaData[1,"rho"]*100)
   write.csv(estDF,paste0("./",sampleid,"_estimated_contamination.csv"), row.names = FALSE)
   
   if(counts2disk){
      message("Saving corrected counts to disk")
      DropletUtils:::write10xCounts(paste0("./",sampleid,"_raw_counts_corrected.h5"), corrected_mtx, type = "HDF5")
      
   }
   
   return(corrected_mtx)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Estimate the filtering parameters (MT, nFeature, nCount)
#'
#' @param sampleID SampleID
#' @param srtObject Unfilterd seurat object
#' @param outDIR Save figures
#' @param apply Apply the calculated filters | Default \code{FALSE}
#' @details Details
#' 
#' The thresolds are set using the MAD (Median absolute deviation), and defined as:
#' 
#' * MAX = median(log10(metric)) + 3*mad(log10(metric))
#' * MIN = median(log10(metric)) - 3*mad(log10(metric))
#'
#' Refer to the pre-processing article for a detailed usage and examples
#' 
#' @keywords soupX
#' @export
#' @examples
#' \dontrun{
#' filt_st <- calculate_filtering_param("45L2", unfilt_st, ".", apply = TRUE)
#' }
#' @md
calculate_filtering_param <- function(sampleID, srtObject, outDIR, apply = FALSE){
   
   library(Seurat)
   library(ggplot2)
   library(ggExtra)
   library(cowplot)
   library(reticulate)
   library(dplyr)
   
   qc_stats <- srtObject@meta.data
   
   message("* Calculating Mitochondrial content")
   
   max_mito_thr <- median(qc_stats$percent.mt) + 3*mad(qc_stats$percent.mt)
   min_mito_thr <- median(qc_stats$percent.mt) - 3*mad(qc_stats$percent.mt)
   if(min_mito_thr < 0){min_mito_thr <- 0}
   
   p1 <- ggplot(qc_stats, aes(x=nFeature_RNA, y=percent.mt)) +
      geom_point() +
      geom_hline(aes(yintercept = max_mito_thr), colour = "red", linetype = 2) +
      geom_hline(aes(yintercept = min_mito_thr), colour = "red", linetype = 2) +
      annotate(geom = "text", label = paste0(as.numeric(table(qc_stats$percent.mt > max_mito_thr | qc_stats$percent.mt < min_mito_thr)[2])," cells removed\n",
                                             as.numeric(table(qc_stats$percent.mt > max_mito_thr | qc_stats$percent.mt < min_mito_thr)[1])," cells remain"),
               x = max(qc_stats$nFeature_RNA)/2 , y = max(qc_stats$percent.mt)/2) +
      theme_classic()
   
   p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) 
   ggsave(p1, filename = paste0(outDIR,"/",sampleID,"/1_mt_content.png"), dpi = 300, height = 9, width = 14, bg="white")
   
   
   message("* Calculating genes and umi content")
   
   min_Genes_thr <- median(log10(qc_stats$nFeature_RNA)) - 3*mad(log10(qc_stats$nFeature_RNA))
   max_Genes_thr <- median(log10(qc_stats$nFeature_RNA)) + 3*mad(log10(qc_stats$nFeature_RNA))
   max_nUMI_thr <- median(log10(qc_stats$nCount_RNA)) + 3*mad(log10(qc_stats$nCount_RNA))
   
   p2 <- ggplot(qc_stats, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point() +
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min_Genes_thr), colour = "green", linetype = 2) +
      geom_hline(aes(yintercept = max_Genes_thr), colour = "green", linetype = 2) +
      geom_vline(aes(xintercept = max_nUMI_thr), colour = "red", linetype = 2) +
      theme_classic()
   p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
   ggsave(p2, filename = paste0(outDIR,"/",sampleID,"/2_gene-umi_content.png"), dpi = 300, height = 9, width = 14, bg="white")
   
   qc_stats <- qc_stats %>% filter(percent.mt < max_mito_thr) %>% filter(percent.mt > min_mito_thr)
   qc_stats <- qc_stats %>% filter(log10(nFeature_RNA) > min_Genes_thr) %>% filter(log10(nCount_RNA) < max_nUMI_thr)
   
   
   message("* Finalising the filtering")
   
   lm.model <- lm(data = qc_stats, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
   p3 <- ggplot(qc_stats, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point() +
      theme_classic() + 
      geom_smooth(method="lm") +
      geom_hline(aes(yintercept = min_Genes_thr), colour = "green", linetype = 2) +
      geom_hline(aes(yintercept = max_Genes_thr), colour = "green", linetype = 2) +
      geom_vline(aes(xintercept = max_nUMI_thr), colour = "red", linetype = 2) +
      geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
      annotate(geom = "text", label = paste0(dim(qc_stats)[1], " QC passed cells"),
               x = max(log10(qc_stats$nCount_RNA))-1 , y = max(log10(qc_stats$nFeature_RNA)))
   ggMarginal(p3, type = "histogram", fill="lightgrey")
   ggsave(p3, filename = paste0(outDIR,"/",sampleID,"/3_main_cells.png"), dpi = 300, height = 9, width = 14, bg="white")
   
   qc_stats$valideCells <- log10(qc_stats$nFeature_RNA) > (log10(qc_stats$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
   
   p4 <- ggplot(qc_stats, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point(aes(colour = valideCells)) +
      geom_smooth(method="lm") +
      theme_classic() + 
      geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") + 
      theme(legend.position="none") +
      annotate(geom = "text", label = paste0(as.numeric(table(qc_stats$valideCells)[2]), " QC passed cells\n",
                                             as.numeric(table(qc_stats$valideCells)[1]), " QC filtered"), x = 4, y = 3.8)
   p4 <- ggMarginal(p4, type = "histogram", fill="lightgrey")
   ggsave(p4, filename = paste0(outDIR,"/",sampleID,"/4_valid_cells.png"), dpi = 300, height = 9, width = 14, bg="white")
   
   qc_stats <- qc_stats %>% filter(valideCells)
   qc_stats$umi <- rownames(qc_stats)
   
   thresholds <- data.frame(SampleID = sampleID,
                            MT_max = round(max_mito_thr, digits = 3),
                            MT_min = round(min_mito_thr, digits = 3),
                            Genes_max_log10 = round(max_Genes_thr, digits = 3),
                            Genes_min_log10 = round(min_Genes_thr, digits = 3),
                            UMI_max_log10 = round(max_nUMI_thr, digits = 3),
                            unfiltered_cells = nrow(srtObject@meta.data),
                            filtered_cells = nrow(qc_stats),
                            removed_cells = nrow(srtObject@meta.data) - nrow(qc_stats))
   
   write.csv(thresholds, file = paste0(outDIR,"/",sampleID,"/5_calculated_thresholds.csv"), row.names = FALSE)
   
   if(apply){
      message("* Applying the filters")
      
      Idents(srtObject) <- "umi"
      srtObject <- subset(srtObject, cells = rownames(qc_stats) , invert = FALSE)
      
      return(srtObject)
   } else {
      return(thresholds) 
   }
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Load 10x GEM data
#' 
#' This function is used to reads the output from cellranger and return the matrices.
#' 
#' @param rawDIR Directory containing the data to load
#' @param sample_name Name of the sample | Optional
#' @param source Use the data from the \code{raw} or \code{filtered} cellranger matrix | Default \code{raw}
#' @param source_type | Default \code{dir}
#' @param input_matrix 
#' @keywords seurat
#' @return 
#' Returns 
#' \item{}{a matrix or a list, if multiple type are included}
#' @export
#' @examples
#' \dontrun{
#' raw_feature_matrix <- load_10x_gem("/mnt/disk2/raw_files/scanpy_test/Biopsie_17", sample_name = "Biopsie_17", source = "raw", source_type="dir")
#' }
#' 
load_10x_gem <- function(rawDIR, source = "raw", source_type="dir", input_matrix = NULL, fixed = FALSE){
   
   library(Seurat)
   
   message("** Reading the parameters ")
   message(paste0(" data directory: ",rawDIR))

   message("** Loading the data ")
   
   if(fixed){
      bc_matrix_name <- paste0("sample_",source,"_feature_bc_matrix")
   } else {
      bc_matrix_name <- paste0(source,"_feature_bc_matrix")
   }
   
   valid_sources <- c("dir","outs","h5","file","file_h5")
   
   if(!source_type %in% valid_sources){
      stop(paste0(source_type," is not a recognised 'source_type' value. Accepted values are: outs, dir, h5, and file"))
      
   } else {
      
      if(source_type == "outs"){
         message(" Loading data from ",paste0(rawDIR,"/outs/",bc_matrix_name))
         raw_feature_matrix <- Read10X(paste0(rawDIR,"/outs/",bc_matrix_name))
         
      } else if(source_type == "dir"){
         message(" Loading data from ",paste0(rawDIR,"/",bc_matrix_name))
         raw_feature_matrix <- Read10X(paste0(rawDIR,"/",bc_matrix_name))
         
      } else if(source_type == "file"){
         input_file <- paste0(rawDIR,"/",input_matrix)
         
         if(exists("input_file") && file.exists(paste(input_file))){
            
            message(" Loading data from ",paste0(rawDIR,"/",input_matrix))
            raw_feature_matrix <- read.delim(paste0(rawDIR,"/",input_matrix), header = TRUE)
            
         } else {
            stop(" File ",input_matrix," not found")
         }
         
      } else if(source_type == "file_h5"){
         input_file <- paste0(rawDIR,"/",input_matrix)
         
         if(exists("input_file") && file.exists(paste(input_file))){
            
            message(" Loading data from ",paste0(rawDIR,"/",input_matrix))
            raw_feature_matrix <- Read10X_h5(input_file)
            
         } else {
            stop(" File ",input_matrix," not found")
         }
         
      } else {
         message(" Loading data from ",paste0(rawDIR,"/",source,"_feature_bc_matrix.h5"))
         raw_feature_matrix <- Read10X_h5(paste0(rawDIR,"/",source,"_feature_bc_matrix.h5"))
      }
   }
   
   if(class(raw_feature_matrix) == "list"){
      matrices <- names(raw_feature_matrix)
      message(paste0(" Multiple matrices found. Returning a list with ",length(matrices)," elements: "))
      print(matrices[1:length(matrices)])
      
   }
   
   return(raw_feature_matrix)
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Run DoubletFinder
#' 
#' This function takes a seurat object and calculates doublet/singlet content using DoubletFinder 
#' 
#' @param srtObject An unfiltered seurat object
#' @param meta_column Column metadata containing the cell clustering
#' @param isSCT SCT normalisation | Default \code{FALSE}
#' @keywords seurat
#' @return 
#' Returns 
#' \item{}{a seurat object}
#' @export
#' @examples
#' \dontrun{
#' sample_rds <- run_doubletFinder(sample_rds)
#' }
#' 
run_doubletFinder <- function(srtObject, meta_column = "seurat_clusters", isSCT = FALSE, countMatrix = FALSE, seurat5 = FALSE){
   
   library(DoubletFinder)
   
   if(countMatrix){
      #options(Seurat.object.assay.version = "v3")
      tmp_data <- CreateSeuratObject(counts = srtObject, project = "doublet_detection", min.cells=3, min.features = 100)
      tmp_data <- NormalizeData(tmp_data)
      tmp_data <- FindVariableFeatures(tmp_data, selection.method = "vst", nfeatures = 3000)
      tmp_data <- ScaleData(tmp_data)
      tmp_data <- RunPCA(tmp_data)
      tmp_data <- RunUMAP(tmp_data, dims = 1:30)
      tmp_data <- FindNeighbors(tmp_data, dims = 1:30)
      tmp_data <- FindClusters(tmp_data)
      
   } else {
      tmp_data <- srtObject
   }
   
   sweep.res.list <- paramSweep_v3(tmp_data, PCs = 1:30, sct = isSCT)
   sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
   bcmvn <- find.pK(sweep.stats)
   
   annotations <- tmp_data@meta.data[,meta_column]
   homotypic.prop <- modelHomotypic(annotations)         ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
   nExp_poi <- round(0.075*nrow(tmp_data@meta.data))     ## Assuming 7.5% doublet formation rate - tailor for your dataset
   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
   
   tmp_data <- doubletFinder_v3(tmp_data, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = isSCT)

   return(tmp_data)
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Function that takes the annotation from one dataset and it transfers it to another one.
#'
#' @param reference_obj Annotated seurat object to use as the reference.
#' @param query_obj Query seurat object.
#' @param pipeline Normalisation method to use (SCT or LogNorm)
#' @param pcDims Number of PCs (30 by default)
#' @param annotation_col Name of the column in the metadata slot of the reference object containing the labels to transfer.
#' @param projectUMAP | Default \code{FALSE}
#' @param isVerbose Show all the messages to screen | Default \code{FALSE}
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{transfer_annotation()}
transfer_annotation <- function(referece_obj, query_obj, pipeline, annotation_col, pcDims = 30,
                                projectUMAP=FALSE, isVerbose= FALSE){
   
   object_name <- deparse(substitute(referece_obj))
   
   annotation_col <- paste0(object_name,"$",annotation_col)
   annotation_col <- eval(parse(text = annotation_col))
   
   message("* Finding transfer anchors")
   ref_anchors <- FindTransferAnchors(reference = referece_obj, query = query_obj,
                                      normalization.method = pipeline,
                                      dims = 1:30, reference.reduction = "pca",
                                      verbose = TRUE)
   
   message("* Transfering annotation")
   
   if(projectUMAP){
      
      message("** Methodology not yet validated, please do not use this results in production yet.")
      # query_obj <- TransferData(anchorset = ref_anchors, 
      #                           reference = referece_obj,
      #                           query = query_obj, 
      #                           refdata = annotation_col,
      #                           dims = 1:pcDims,
      #                           k.weight = 10)
      # 
      # query_obj <- IntegrateEmbeddings(anchorset = ref_anchors, 
      #                                  reference = referece_obj,
      #                                  query = query_obj, new.reduction.name = "ref.pca")
      # 
      # referece_obj <- RunUMAP(referece_obj, dims = 1:30,  reduction = "pca", return.model = TRUE)
      # 
      # query_obj <- ProjectUMAP(query = query_obj, 
      #                          reference = referece_obj,
      #                          query.reduction = "ref.pca", 
      #                          reference.reduction = "pca", reduction.model = "umap")
      referece_obj <- RunUMAP(referece_obj, dims = 1:30, reduction = "pca", return.model = TRUE)
      query_obj <- MapQuery(anchorset = ref_anchors, reference = referece_obj, query = query_obj,
                            refdata = list(celltype = annotation_col), reference.reduction = "pca", reduction.model = "umap")
      
   } else {
      predictions <- TransferData(anchorset = ref_anchors, 
                                  refdata = annotation_col,
                                  dims = 1:pcDims,
                                  k.weight = 10)
      
      query_obj <- AddMetaData(query_obj, metadata = predictions)
   }
   
   return(query_obj)
   
}
   

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Dump umap coordinates and metadata to load in the loupe browser
#' 
#' @param isVerbose Show all the messages to screen | Default \code{FALSE}
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{get_loupe_files()}
get_loupe_files <- function(srtObject,cellranger_cells,meta_cols, outDIR = NULL){
   
   # v.1.0
   # srtObject = Seurat object.size
   # cellranger_cells = List of cells in the loupe file
   # meta_cols = columns from srtObject@meta.data to include in the metadata file
   
   umapCoordinates <- as.data.frame(Embeddings(object = srtObject, reduction = "umap"))
   samples <- unique(srtObject@meta.data$SampleID)
   
   coordinate_fixed <- data.frame()
   for(idx in 1:length(samples)){
      
      tmp_table <- umapCoordinates[grep(samples[idx],rownames(umapCoordinates)),]
      tmp_table$Barcode <- rownames(tmp_table)
      tmp_table <- tmp_table[c("Barcode","UMAP_1","UMAP_2")]
      rownames(tmp_table) <- NULL
      tmp_table$Barcode <- gsub(paste0("^",samples[idx],"_"),"",tmp_table$Barcode)
      tmp_table$Barcode <- gsub("-1",paste0("-",idx),tmp_table$Barcode)
      tmp_table$LibraryID <- samples[idx]
      
      coordinate_fixed <- rbind(coordinate_fixed,tmp_table)
   }
   
   coordinate_fixed$Barcode_loupe <- paste0(paste0(coordinate_fixed$LibraryID,"_",coordinate_fixed$Barcode))
   coordinate_fixed$Barcode_loupe <- sub("-2$","-1",coordinate_fixed$Barcode_loupe)
   coordinate_fixed$Barcode_loupe <- sub("-3$","-1",coordinate_fixed$Barcode_loupe)
   coordinate_fixed$Barcode_loupe <- sub("-4$","-1",coordinate_fixed$Barcode_loupe)
   
   # Metadata file
   coordinate_metadata <- data.frame(Barcode = coordinate_fixed$Barcode,
                                     LibraryID = coordinate_fixed$LibraryID,
                                     Barcode_loupe = paste0(coordinate_fixed$LibraryID,"_",coordinate_fixed$Barcode))
   
   # Add columns from seurat@meta.data
   coordinate_metadata <- cbind(coordinate_metadata, srtObject@meta.data[,meta_cols])
   rownames(coordinate_metadata) <- NULL
   coordinate_metadata$Barcode_loupe <- sub("-2$","-1",coordinate_metadata$Barcode_loupe)
   coordinate_metadata$Barcode_loupe <- sub("-3$","-1",coordinate_metadata$Barcode_loupe)
   coordinate_metadata$Barcode_loupe <- sub("-4$","-1",coordinate_metadata$Barcode_loupe)
   
   # Select only barcodes in cellranger
   cellranger_cells$Barcode_loupe <- paste0(cellranger_cells$LibraryID,"_",cellranger_cells$Barcode)
   cellranger_cells$Barcode_loupe <- sub("-2$","-1",cellranger_cells$Barcode_loupe)
   cellranger_cells$Barcode_loupe <- sub("-3$","-1",cellranger_cells$Barcode_loupe)
   cellranger_cells$Barcode_loupe <- sub("-4$","-1",cellranger_cells$Barcode_loupe)
   
   coordinate_fixed <- coordinate_fixed[which(coordinate_fixed$Barcode_loupe %in% cellranger_cells$Barcode_loupe),]
   coordinate_metadata <- coordinate_metadata[which(coordinate_metadata$Barcode_loupe %in% cellranger_cells$Barcode_loupe),]
   
   if(nrow(cellranger_cells) < nrow(coordinate_fixed)){
      message("Error: The Seurat object has more cells than the loupe file. Exit.")
      
   } else {
      # Save both files
      write.csv(coordinate_fixed,paste0(outDIR,"/loupe_coordinates_umap.csv"),row.names = FALSE)
      write.csv(coordinate_metadata,paste0(outDIR,"/loupe_coordinates_metadata.csv"),row.names = FALSE)
      
   }
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Add ADT information to a seurat object
#' 
#' @param srtObject Seurat object
#' @param exp_metadata Metadata table
#' @param dataDIR Directory containing the cellranger output directory
#' @param outs Outs dir or not | Default \code{FALSE}
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{add_adt_assay()}
add_adt_assay <- function(srtObject, exp_metadata, dataDIR, outs=FALSE){
   
   # srtObject <- exp_Integrated
   # exp_metadata <- exp_meta
   # 
   adt_matrix <- data.frame()
   
   for(idx in 1:nrow(exp_metadata)){
      
      tmp_sample <- exp_metadata[idx,"SampleID"]
      tmp_sampleName <- exp_metadata[idx,"SampleID_raw"]
      
      message(paste0("Processing sample ",tmp_sample))
      
      if(outs){
         tmp_feature_matrix <- Read10X(paste0(dataDIR,"/",tmp_sampleName,"/outs/raw_feature_bc_matrix"))  
      } else {
         message(paste0(dataDIR,"/",tmp_sample,"/raw_feature_bc_matrix"))
         tmp_feature_matrix <- Read10X(paste0(dataDIR,"/",tmp_sampleName,"/raw_feature_bc_matrix"))   
      }
      
      #tmp_feature_matrix <- Read10X(paste0(dataDIR,"/",tmp_sampleDIR,"/filtered_feature_bc_matrix")) 
      tmp_feature_matrix <- tmp_feature_matrix[["Antibody Capture"]]
      colnames(tmp_feature_matrix) <- paste0(tmp_sample,"_",colnames(tmp_feature_matrix))
      
      if(idx == 1){
         adt_matrix <- rbind(adt_matrix,as.data.frame(tmp_feature_matrix))
      } else {
         adt_matrix <- cbind(adt_matrix,tmp_feature_matrix)
      }
      
   }
   
   rnaCells <- colnames(srtObject@assays$RNA)
   adtCells <- colnames(adt_matrix)
   keep_adt <- colnames(adt_matrix)[which(adtCells %in% rnaCells)]
   
   srtObject@meta.data$has_ADT <- sapply(rownames(srtObject@meta.data),
                                         function(ita) ifelse(ita %in% keep_adt,"yes","no"))
   
   print(table(srtObject@meta.data$has_ADT, srtObject@meta.data$SampleID ))
   
   Idents(srtObject) <- "has_ADT"
   adtObject <- subset(srtObject, idents = "yes")
   
   adt_matrix <- adt_matrix[,which(colnames(adt_matrix) %in% keep_adt)]
   adt_assay <- CreateAssayObject(counts = adt_matrix)
   adtObject[["ADT"]] <- adt_assay
   
   DefaultAssay(adtObject) <- "ADT"
   adtObject <- NormalizeData(adtObject, normalization.method = "CLR", margin = 2)
   DefaultAssay(adtObject) <- "RNA"
   
   return(adtObject)
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Add VDJ information to a seurat object
#'
#' @param isVerbose Show all the messages to screen | Default \code{FALSE}
#' @keywords seurat
#' @export
#' @examples 
#' \dontrun{add_vdj_assay()}
add_vdj_assay<- function(srtObject,exp_metadata,dataDIR){
   message("write something")
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Expansion of CellCycleScoring from the seurat package. 
#' 
#' In addition to the default phases, CellCycleScoring_extra supports 2 more cell phases. 
#'
#' @param object Seurat object
#' @param s.features List of genes in the S phase
#' @param g2m.features List of genes in the G2M phase
#' @param g1s.features List of genes in the G1S phase
#' @param m.features List of genes in the M phase
#' @param mg1.features List of genes in the MG1 phase
#' @param ctrl Number of control features selected from the same bin per analyzed feature supplied to AddModuleScore. Defaults to value equivalent to minimum number of features present in 's.features' and 'g2m.features'.
#' @param set.ident If true, sets identity to phase assignments Stashes old identities in 'old.ident'
#' @param ... Arguments to be passed to AddModuleScore
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{CellCycleScoring_extra()}
CellCycleScoring_extra <- function(
      object,
      s.features,
      g2m.features,
      g1s.features,
      m.features,
      mg1.features,
      ctrl = NULL,
      set.ident = FALSE,
      ...
) {
   name <- 'Cell.Cycle'
   features <- list('G1S.Score' = g1s.features,
                    'S.Score' = s.features, 
                    'G2M.Score' = g2m.features,
                    'M.Score' = m.features,
                    'MG1.Score' = mg1.features)
   
   if (is.null(x = ctrl)) {
      ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
   }
   object.cc <- AddModuleScore(object = object,features = features, name = name, ctrl = ctrl)
   
   cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), value = TRUE)
   cc.scores <- object.cc[[cc.columns]]
   
   rm(object.cc)
   CheckGC()
   assignments <- apply(
      X = cc.scores,
      MARGIN = 1,
      FUN = function(scores,first='G1S',second='S',third='G2M',fourth='M',fifth='MG1',null='G1') {
         if (all(scores < 0)) {
            return(null)
         } else {
            if (length(which(x = scores == max(scores))) > 1) {
               return('Undecided')
            } else {
               return(c(first, second,third,fourth,fifth)[which(x = scores == max(scores))])
            }
         }
      }
   )
   cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
   colnames(x = cc.scores) <- c('rownames','G1S.Score','S.Score','G2M.Score','M.Score','MG1.Score','Phase')
   rownames(x = cc.scores) <- cc.scores$rownames
   cc.scores <- cc.scores[, c('G1S.Score','S.Score','G2M.Score','M.Score','MG1.Score','Phase')]
   object[[colnames(x = cc.scores)]] <- cc.scores
   if(set.ident) {
      object[['old.ident']] <- Idents(object = object)
      Idents(object = object) <- 'Phase'
   }
   return(object)
}
