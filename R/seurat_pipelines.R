
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Check sample quality and filter the count matrix
#'
#' This function is used to reads the output from cellranger, create a seurat object, 
#' apply basic filtering, and normalise and cluster the filtered counts.
#' 
#' @param exp_meta Metadata table (\code{data.frame} by default)
#' @param sample_num Row number of the metadata file containing the sample information
#' @param pipeline Perform SCT or lognorm normalisation
#' @param minFeatures Maximum number of genes to filter out
#' @param pcDims Number of PCs (30 by default)
#' @param source Use the data from the \code{raw} or \code{filtered} cellranger matrix (\code{raw} by default)
#' @param source_type \code{dir} by default
#' @param organism Set the sample organism (\code{human} by default)
#' @param isVerbose Print all the messages on screen (\code{FALSE} by default)
#' @param vst_flavor sctransform version
#' @keywords seurat
#' @return 
#' Returns 
#' \item{}{a folder, named as the sample ID, containing QC figures and tables}
#' \item{}{a seurat object with filtering and clustering applied}
#' @export
#' @examples
#' \dontrun{
#' run_qc_and_filtering(exp_meta,number,".","SCT", 500, 30, source = "raw", source_type = "outs", isVerbose=FALSE)
#' }
#' 
run_qc_and_filtering <- function(exp_meta,sample_num,pipeline, minFeatures, pcDims = 30, mtPerc = 20, source = "raw",
                                source_type = "dir", organism = "human", isVerbose = FALSE, vst_flavor = NULL, input_matrix = NULL,
                                isFixed = FALSE, remove_doublets = FALSE, mtPattern = NULL, riboPattern = NULL, 
                                seurat5 = FALSE, bpcells_dir = NULL, labeling_celline = FALSE){
   
   library(Seurat)
   library(BPCells)
   #library(Azimuth)
   library(patchwork)
   library(wesanderson)
   library(ggplot2)
   library(cli)
   
   ggOriginal <- theme_set(theme_classic())
   
   sample_name <- as.character(exp_meta[sample_num,"SampleID"])
   dataDIR <- as.character(exp_meta[sample_num,"dataDIR"])
   
   cli_inform(message = paste0("** Processing sample ",sample_name))
   cli_inform(message = " Loading 10x raw matrix...")
   
   qc_dir <- paste0("./",sample_name,"_",pipeline)
   if(!exists(qc_dir)){
      dir.create(qc_dir)
   }
   
   raw_feature_matrix <- load_10x_gem(rawDIR = dataDIR, source = source, source_type=source_type, fixed = isFixed)
   
   if(class(raw_feature_matrix) == "list"){
      raw_st <- raw_feature_matrix[["Gene Expression"]]
      multipleAssays <- TRUE
   } else {
      raw_st <- raw_feature_matrix
      multipleAssays <- FALSE
   }
   colnames(raw_st) <- gsub("^",paste0(sample_name,"_"),colnames(raw_st))
   write.csv(data.frame(raw_umi = ncol(raw_st), raw_genes = nrow(raw_st)), file = paste0(qc_dir,"/",sample_name,"_raw_cells.csv"), row.names = FALSE)
      
   if(remove_doublets){
      
      cli_inform(message = " Setting seurat v5 to FALSE")
      options(Seurat.object.assay.version = "v5")
      doublets <- run_doubletFinder(raw_st, meta_column = "seurat_clusters", isSCT = FALSE, countMatrix = TRUE, seurat5 = FALSE)
      colnames(doublets@meta.data)[ncol(doublets@meta.data)] <- "isDoublet"
         
      #Idents(doublets) <- "isDoublet"
      colorini <- c("#FDE725FF","#440154FF")
      names(colorini) <- c("Singlet","Doublet")
      
      p0 <- DimPlot(doublets, group.by = "isDoublet", order = TRUE, cols = colorini)
      ggsave(p0, filename = paste0(qc_dir,"/",sample_name,"_doublets.png"), width = 12, height = 9)
      
      doublets <- subset(doublets, idents = "Singlet", invert = FALSE)
      if(seurat5){
         raw_st <- doublets@assays$RNA@layers$counts
      } else {
         raw_st <- doublets@assays$RNA@counts
         
      }
      
      write.csv(data.frame(singlets_umi = ncol(raw_st), singlets_genes = nrow(raw_st)), 
                file = paste0(qc_dir,"/",sample_name,"_singlets_cells.csv"), row.names = FALSE)
      
   }
   
   if(seurat5 && !is.null(bpcells_dir)){
      cli_inform(message = " Setting seurat v5 to TRUE")
      cli_inform(message = " Using BPCells capabilities")
      options(Seurat.object.assay.version = "v5")
      raw_feature_matrix <- convert_matrix_type(raw_st, type = "uint32_t")
      write_matrix_dir(raw_feature_matrix, dir = paste0(bpcells_dir,"/tmp_data/",sample_name), overwrite = TRUE)
      raw_st <- open_matrix_dir(dir = paste0(bpcells_dir,"/tmp_data/",sample_name))
      
   } else if(seurat5){
      cli_inform(message = " Setting seurat v5 to TRUE")
      options(Seurat.object.assay.version = "v5")
      
   }
   
   cli_inform(message = " Creating the unfiltered seurat object...")
   unfilt_st <- CreateSeuratObject(counts = raw_st, project = sample_name, min.cells=3, min.features = 100)
   unfilt_st <- calculate_mito_ribo(unfilt_st, organism, new_mito = mtPattern, new_ribo = riboPattern)
   unfilt_st <- calculate_cell_complexity(unfilt_st)
   unfilt_st@meta.data$SampleID <- sample_name
   
   cli_inform(message = " Creating the filtered seurat object...")
   filt_st <- CreateSeuratObject(counts = raw_st, project = sample_name, min.cells=3, min.features = minFeatures)
   filt_st <- calculate_mito_ribo(filt_st, organism, new_mito = mtPattern, new_ribo = riboPattern)
   filt_st <- calculate_cell_complexity(filt_st)
   filt_st <- subset(filt_st, subset = percent.mt < mtPerc)
   filt_st@meta.data$SampleID <- sample_name
   
   cli_inform(message = " Generating QC figures...")
   
   unfilt_st@meta.data$source <- paste0(unfilt_st@meta.data$SampleID,"_pre-filtered") 
   filt_st@meta.data$source <- paste0(filt_st@meta.data$SampleID,"_filtered") 
   
   Idents(unfilt_st) <- "source"
   Idents(filt_st) <- "source"
   
   # Generate QC figures and compare filtered/unfiltered data
   compare_qc(unfilt_st, filt_st, sample_name, colori = NULL, qc_dir = qc_dir, minFeatures = minFeatures, mtPerc = mtPerc)
   
   cli_inform(message = " Generating QC table...")
   
   unfilt_qc_table <- data.frame(SampleID = sample_name,
                                 pipeline = pipeline,
                                 Status = "Pre-filtered",
                                 nCells = nrow(unfilt_st@meta.data),
                                 nFeature_median = round(median(unfilt_st@meta.data$nFeature_RNA)),
                                 nCount_median = round(median(unfilt_st@meta.data$nCount_RNA)),
                                 MT_median = round(median(unfilt_st@meta.data$percent.mt)),
                                 Ribo_median = round(median(unfilt_st@meta.data$percent.ribo)), 
                                 Complexity_median = median(unfilt_st@meta.data$log10GenesPerUMI))
   
   filt_qc_table <- data.frame(SampleID = sample_name,
                               pipeline = pipeline,
                               Status = "Filtered",
                               nCells = nrow(filt_st@meta.data),
                               nFeature_median = round(median(filt_st@meta.data$nFeature_RNA)),
                               nCount_median = round(median(filt_st@meta.data$nCount_RNA)),
                               MT_median = round(median(filt_st@meta.data$percent.mt)),
                               Ribo_median = round(median(filt_st@meta.data$percent.ribo)),
                               Complexity_median = median(filt_st@meta.data$log10GenesPerUMI))

   
   qc_table <- rbind(unfilt_qc_table,filt_qc_table)
   write.csv(qc_table,paste0(qc_dir,"/",sample_name,"_",pipeline,"_qc_summary.csv"), row.names = FALSE)
   
   cli_inform(message = " Normalising filtered data ...")
   
   if(pipeline == "SCT"){
      cli_inform(message = " NOTE: Using SCTransform ...")
      filt_st <- SCTransform(filt_st, vst.flavor = vst_flavor, vars.to.regress = c("percent.mt"), verbose = isVerbose)

   } else {
      cli_inform(message = "  NOTE: SCTransform NOT in use ...")
      filt_st <- NormalizeData(filt_st, verbose = isVerbose)
      filt_st <- FindVariableFeatures(filt_st, selection.method = "vst", nfeatures = 3000, verbose= isVerbose)
      filt_st <- ScaleData(filt_st, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = isVerbose)

   }

   cli_inform(message = " Running PCA and uMAP ..")
   filt_st <- RunPCA(filt_st, verbose = isVerbose)
   filt_st <- RunUMAP(filt_st, dims = 1:pcDims, verbose = isVerbose)

   cli_inform(message = " Clustering ...")
   filt_st <- FindNeighbors(filt_st, dims = 1:pcDims, verbose = isVerbose)
   filt_st <- FindClusters(filt_st, verbose = isVerbose)
   
   sample_umap <- DimPlot(filt_st, label = TRUE) + NoLegend()
   ggsave(sample_umap, filename = paste0(qc_dir,"/",sample_name,"_",pipeline,"_umap.png"), dpi = 300, width = 14, height = 9)

   if(labeling_celline){
      
      library(scGate)
      
      cli_inform(message = " Finding the celline UMIs")
      
      celline_sig <- list(K562 = c("MAGEA1+","MAGEA12+","MAGEA3+","MAGEA4+","PRAME+","MKI67+","EPCAM+","PTPRC-"))
      
      filt_st <- calculate_genes_signature(filt_st,celline_sig, assay = "RNA")
      
      DefaultAssay(filt_st) <- "RNA"
      Idents(filt_st) <- "seurat_clusters"
      
      cli_inform(message = " Gating the celline UMIs")
      
      my_scGate_model <- gating_model(name = "K562", signature = celline_sig$K562)  
      #filt_st <- scGate(filt_st, model = my_scGate_model)
      
      c1 <- DimPlot(filt_st, group.by = "seurat_clusters", label = TRUE) + NoLegend()
      c2 <- FeaturePlot(filt_st, "K562_kNN", order = TRUE, cols = c("#FDE725FF","#440154FF"))
      c3 <- DotPlot(filt_st, group.by = "seurat_clusters",features = "K562_kNN", cols = c("#FDE725FF","#440154FF"))
      #c4 <- DimPlot(filt_st, group.by = "is.pure", label = TRUE) + NoLegend()
      combined_plot <- (c1|c2|c3)
      ggsave(combined_plot, filename = paste0(qc_dir,"/K562_cells.png"), bg="white", dpi = 300, width = 18, height = 9)
      
   }
   
   cli_inform(paste0(" Saving seurat object: ./filtered_",sample_name,"_",pipeline,".rds"))

   if(seurat5 && !is.null(bpcells_dir)){
      options(Seurat.object.assay.version = "v5")
      SaveSeuratRds(filt_st, file = paste0(bpcells_dir,"/filtered_",sample_name,"_",pipeline,".rds"), destdir = bpcells_dir, relative = TRUE)
   } else {
      saveRDS(filt_st,paste0("./filtered_",sample_name,"_",pipeline,".rds"))
   }
   
   cli_inform(message = " Clean up...")
   rm(unfilt_st)
   rm(filt_st)
   
   cli_inform(message = " Done!")
   message(message = "..............................")

}


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Samples integration using seurat v5 (logNorm only)
#'
#' @param dataDIR Path to filtered rds files. 
#' @param intMethod d
#' @param outname Name of the intergated object.
#' @param nfeatures Number of genes to use when integrating (3000 by default)
#' @param pcDims Number of PCs (30 by default)
#' @param isVerbose Show all the messages to screen.FALSE by default
#' @param vst_flavor sctransform version
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{
#' run_seurat_integration(".","SCT", "cca", "integration1", 3000, 30, FALSE , "v2")
#' }
run_seurat_integration_v5 <- function(dataDIR, samples_vect = NULL, intMethod, outname, nfeatures = 3000, 
                                      pcDims = 30, isVerbose = FALSE, isSCT = FALSE, vst_flavor = NULL, ref_samples = NULL,
                                      merge_layers = TRUE, seurat5 = TRUE, bpcells_dir=NULL, checkpoint = FALSE){
   
   options(Seurat.object.assay.version = "v5")
   message("** Loading samples data: ")
   
   if(is.null(samples_vect)){
      samples_vect <- list.files(dataDIR, full.names = FALSE, recursive = FALSE, pattern = "rds") 
   } 
   
   sample_list <- c()
   metadata_list <- data.frame(matrix(ncol = 9, nrow = 0))
   
   if(isSCT){
      colonnine <- c("nCount_SCT","nFeature_SCT")
   } else {
      colonnine <- c("nCount_RNA","nFeature_RNA")
   }
   
   names(metadata_list) <- c("orig.ident",colonnine, "percent.mt","percent.ribo",
                             "log10GenesPerUMI","SampleID","source","seurat_clusters")
   
   for(srt_object in samples_vect){
      
      message(paste0(" ",srt_object))
      
      rds_file <- srt_object
      sample_name <- gsub("filtered_","",rds_file)
      sample_name <- gsub("_logNorm.rds","",sample_name)
      sample_name <- gsub("_SCT_K562_removed.rds","",sample_name)
      sample_name <- gsub("_SCT.rds","",sample_name)
      sample_name <- gsub("_logNorm_K562_removed.rds","",sample_name)
      
      rds_file <- paste0(dataDIR,"/",rds_file)
      mySample <- readRDS(rds_file)
      
      # # remove celline columns      
      # if(length(grep("K562|K562_kNN", names(mySample@meta.data))) >0){
      #    mySample@meta.data <- mySample@meta.data[,names(metadata_list)]
      # }
      
      mySample@meta.data <- mySample@meta.data[,names(metadata_list)]
      
      metadata_list <- rbind(metadata_list,mySample@meta.data)
      sample_list[[sample_name]] <- mySample@assays$RNA$counts
      rm(mySample)
      
      gc()
   }
   
   message("** Create layered object")
   integrated_obj <- CreateSeuratObject(counts = sample_list,
                                        meta.data = metadata_list)
   
   if(isSCT){
      norm_method <- "SCT"
   } else {
      norm_method <- "LogNormalize"
   }
   
   if(norm_method == "SCT"){
      message("** Normalising the layered object")
      integrated_obj <- SCTransform(integrated_obj)
      
   } else {
      message("** Normalising the layered object")
      integrated_obj <- NormalizeData(integrated_obj, verbose = isVerbose)
      message("** Finding HVG")
      integrated_obj <- FindVariableFeatures(integrated_obj, nfeatures = nfeatures, verbose = isVerbose)
      message("** Scaling the layered object")
      integrated_obj <- ScaleData(integrated_obj, vars.to.regress = c("SampleID","nCount_RNA"), verbose = isVerbose)
   }
   
   message("** Running PCA")
   integrated_obj <- RunPCA(integrated_obj, verbose = isVerbose)
   gc()
   
   if(checkpoint){
      message("** Saving checkpoint")
      outDIR <- paste0(dataDIR,"/",outname)
      if(seurat5 && !is.null(bpcells_dir)){
         SaveSeuratRds(integrated_obj, file = paste0(outname,"_tmp.Rds"), destdir = outDIR,  relative = TRUE)
      } else if(seurat5){
         saveRDS(integrated_obj, file = paste0(outname,"_tmp.Rds"))
      } 
   }
   
   message("** Integrating")
   
   if(!is.null(ref_samples)){
      
      if(!isSCT){
         ref_samples_idx <- which(Layers(integrated_obj, search = 'data') %in% paste0("data.",ref_samples))
      } else {
         ref_samples_idx <- which(names(integrated_obj@assays$SCT@SCTModel.list) %in% ref_samples)
      }
      
   } else {
      ref_samples_idx <- NULL
   }
   
   print(ref_samples_idx)
   
   if(intMethod == "rpca"){
      new_reduction <- paste0("red_",intMethod)
      new_umap <- paste0("umap_",intMethod)
      clust_name <- "rpca_"
      integrated_obj <- IntegrateLayers(integrated_obj, method = RPCAIntegration, normalization.method =norm_method,
                                        orig.reduction = "pca", new.reduction = new_reduction,
                                        verbose = isVerbose, reference = ref_samples_idx)
      
   } else if(intMethod == "cca"){
      new_reduction <- paste0("red_",intMethod)
      new_umap <- paste0("umap_",intMethod)
      clust_name <- "cca_"
      integrated_obj <- IntegrateLayers(integrated_obj, method = CCAIntegration, normalization.method =norm_method,
                                        orig.reduction = "pca", new.reduction = new_reduction,
                                        verbose = isVerbose, reference = ref_samples_idx)
   } else {
      new_reduction <- paste0("red_",intMethod)
      new_umap <- paste0("umap_",intMethod)
      clust_name <- "harmony_"
      integrated_obj <- IntegrateLayers(integrated_obj, method = HarmonyIntegration,
                                        orig.reduction = "pca", new.reduction = new_reduction,
                                        verbose = isVerbose, reference = ref_samples_idx)
   }
   
   message("** Reclustering")
   integrated_obj <- FindNeighbors(integrated_obj, reduction = new_reduction, dims = 1:pcDims, verbose = isVerbose)
   integrated_obj <- FindClusters(integrated_obj, resolution = c(0.8), verbose = isVerbose, cluster.name = paste0(clust_name,"0.8"))
   integrated_obj <- FindClusters(integrated_obj, resolution = c(1), verbose = isVerbose, cluster.name = paste0(clust_name,"1"))
   integrated_obj <- FindClusters(integrated_obj, resolution = c(1.2), verbose = isVerbose, cluster.name = paste0(clust_name,"1.2"))
   integrated_obj <- FindClusters(integrated_obj, resolution = c(1.4), verbose = isVerbose, cluster.name = paste0(clust_name,"1.4"))
   integrated_obj <- FindClusters(integrated_obj, resolution = c(1.6), verbose = isVerbose, cluster.name = paste0(clust_name,"1.6"))
   integrated_obj <- FindClusters(integrated_obj, resolution = c(1.8), verbose = isVerbose, cluster.name = paste0(clust_name,"1.8"))
   
   integrated_obj <- RunUMAP(integrated_obj, reduction = new_reduction, dims = 1:pcDims, reduction.name = new_umap)
   
   gc()
   
   if(merge_layers){
      if(!isSCT){
         message("** Merging the layers")
         integrated_obj <- JoinLayers(integrated_obj)
      }
   }
   
   message("** Saving the integrated object")
   outDIR <- paste0(dataDIR,"/",outname)
   if(seurat5 && !is.null(bpcells_dir)){
      SaveSeuratRds(integrated_obj, file = paste0(outname,".Rds"), destdir = outDIR,  relative = TRUE)
   } else if(seurat5){
      saveRDS(integrated_obj, file = paste0(outname,".Rds"))
   }
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Samples integration using seurat 
#'
#' This function takes a list of rds files (generated by `run_qc_and_filtering`) and integrates them together.
#' @param dataDIR Path to filtered rds files. 
#' @param pipeline Use SCT or lognorm normalisation.
#' @param srt_reduction Use cca or rpca.
#' @param outname Name of the intergated object.
#' @param nfeatures Number of genes to use when integrating (3000 by default)
#' @param pcDims Number of PCs (30 by default)
#' @param kvalue Set the number of neighbors used when picking and weighting anchors (100 by default) 
#' @param isVerbose Show all the messages to screen.FALSE by default
#' @param vst_flavor sctransform version
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{
#' run_seurat_integration(".","SCT", "cca", "integration1", 3000, 30, FALSE , "v2")
#' }
run_seurat_integration <- function(dataDIR, pipeline, srt_reduction, outname,nfeatures = 3000, 
                                   pcDims = 30, kvalue =100, isVerbose = FALSE, vst_flavor = NULL){
   
   message("** Loading samples data: ")
   
   samples_data <- list.files(dataDIR, pattern = "^filtered")
   samples_data <- samples_data[grep(pipeline,samples_data)]
   
   sample_list <- list()
   for(srt_object in samples_data){
      
      message(paste0(" ",srt_object))
      mySample <- readRDS(paste0(dataDIR,"/",srt_object))
      
      sampleid <- unlist(strsplit(srt_object,split = "\\."))[[1]]
      sampleid <- gsub("filtered_","",sampleid)
      sampleid <- gsub("_standard","",sampleid)
      sampleid <- gsub("_SCT","",sampleid)
      
      sample_list[[sampleid]] <- mySample
      rm(mySample)
      
      gc()
   }
   
   message("** Starting the integraton")
   
   message(" Selecting integration anchors")
   exp_Features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = nfeatures)
   
   if(pipeline == "SCT"){
      
      message(" Preparing integration objects")
      sample_list <- PrepSCTIntegration(object.list = sample_list,
                                        anchor.features = exp_Features, verbose = isVerbose)
      
      normalization_method <- "SCT"
      
   } else {
      
      normalization_method <- "LogNormalize"
      
   }
   
   message(" Finding integration anchors")
   exp_Anchors <- FindIntegrationAnchors(object.list = sample_list,
                                         normalization.method = normalization_method,
                                         reduction = srt_reduction,
                                         k.anchor = kvalue,
                                         dims = 1:pcDims,
                                         scale = TRUE,
                                         anchor.features = exp_Features, verbose = isVerbose)
   message(" Integrating")
   integrated_obj <- IntegrateData(anchorset = exp_Anchors, 
                                   normalization.method = normalization_method, 
                                   k.weight = kvalue, 
                                   dims = 1:pcDims, 
                                   verbose = isVerbose)
   
   message("** Running PCA and uMAP on the integrated object")
   
   DefaultAssay(integrated_obj) <- "integrated"
   
   if(pipeline != "SCT"){
      message("** Scaling")
      integrated_obj <- ScaleData(integrated_obj, vars.to.regress = c("SampleID","nCount_RNA") , verbose = isVerbose)
   }
   
   #integrated_obj <- ScaleData(integrated_obj)
   integrated_obj <- RunPCA(integrated_obj, verbose = isVerbose)
   integrated_obj <- RunUMAP(integrated_obj, dims = 1:pcDims, verbose = isVerbose)
   
   message("** Reclustering")
   integrated_obj <- FindNeighbors(integrated_obj, dims = 1:pcDims, verbose = isVerbose)
   integrated_obj <- FindClusters(integrated_obj, resolution = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2), verbose = isVerbose)
   
   if(pipeline == "SCT" && vst_flavor == "v2"){
      
      message("** SCT v2 selected: Recorrecting the counts")
      integrated_obj <- PrepSCTFindMarkers(integrated_obj, assay = "SCT", verbose = TRUE)
      
   }
   
   message("** Saving the integrated object")
   saveRDS(integrated_obj, paste0(outname,".rds"))
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Samples reference based integration using seurat 
#'
#' This function takes a list of rds files (generated by `run_qc_and_filtering`) and integrates them together.
#' @param dataDIR Path to filtered rds files. 
#' @param pipeline Use SCT or lognorm normalisation.
#' @param srt_reduction Use cca or rpca.
#' @param pcDims Number of PCs (30 by default)
#' @param outname Name of the intergated object.
#' @param ref_samples List of samples to use as reference.
#' @param ref_as_rds Boolean. FALSE by default
#' \itemize{
#' \item If TRUE 
#' \item If FALSE
#'}
#' @param isVerbose Show all the messages to screen.FALSE by default
#' @keywords seurat
#' @export 
#' @examples 
#' \dontrun{
#' run_integration_with_Ref(".", "SCT", "cca", pcDims = 30, "intergation_WRef", c("sample1","sample12"), ref_as_rds = FALSE, isVerbose = FALSE)
#' }
run_integration_with_Ref <- function(dataDIR, pipeline, srt_reduction, outname, ref_samples, ref_as_rds = FALSE, pcDims = 30, isVerbose = FALSE, compress_rds = TRUE){
   
   message("** Loading samples data")
   
   samples_data <- list.files(dataDIR, pattern = "^filtered")
   samples_data <- samples_data[grep(pipeline,samples_data)]
   message("** Processing samples: ")
   
   if(pipeline == "SCT"){
      normalization_method <- "SCT"
      default_assay <- "SCT"
   } else {
      normalization_method <- "LogNormalize"
      default_assay <- "RNA"
   }
   
   sample_list <- list()
   for(srt_object in samples_data){
      
      message(paste0(" ",srt_object))
      mySample <- readRDS(paste0(dataDIR,"/",srt_object))
      
      sampleid <- unlist(strsplit(srt_object,split = "\\."))[[1]]
      sampleid <- gsub("filtered_","",sampleid)
      sampleid <- gsub(paste0("_",pipeline),"",sampleid)
      sampleid <- gsub("_K562_removed","",sampleid)
      
      DefaultAssay(mySample) <- default_assay
      sample_list[[sampleid]] <- mySample
      rm(mySample)
      
      gc()
   }
   
   if(ref_as_rds){
      sample_list[["reference"]] <- ref_samples
      ref_idx <- which(names(sample_list) %in% "reference")
      
   } else {
      ref_idx <- which(names(sample_list) %in% ref_samples)
   }
   
   message("** Starting the integraton")
   
   message(" Selecting integration anchors")
   exp_Features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)
   
   if(pipeline == "SCT"){   
      message(" Preparing integration objects")
      sample_list <- PrepSCTIntegration(object.list = sample_list,
                                        anchor.features = exp_Features, verbose = isVerbose)
      
      normalization_method <- "SCT"
   } else {
      normalization_method <- "LogNormalize"
   }
   
   message(" Finding integration anchors")
   
   exp_Anchors <- FindIntegrationAnchors(object.list = sample_list,
                                         normalization.method = normalization_method,
                                         reduction = srt_reduction,
                                         #k.anchor = 10, 
                                         reference = ref_idx,
                                         anchor.features = exp_Features, verbose = isVerbose)
   
   #saveRDS(exp_Anchors,"./tmp_anchors.rds")
   
   message(" Integrating")
   integrated_obj <- IntegrateData(anchorset = exp_Anchors, 
                                   normalization.method = normalization_method, verbose = isVerbose)
   
   message("** Running PCA and uMAP on the integrated object")
   
   DefaultAssay(integrated_obj) <- "integrated"
   
   if(normalization_method == "LogNormalize"){
      integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
   }
   
   integrated_obj <- RunPCA(integrated_obj, verbose = FALSE)
   integrated_obj <- RunUMAP(integrated_obj, dims = 1:pcDims, verbose = FALSE)
   
   message("** Reclustering")
   integrated_obj <- FindNeighbors(integrated_obj, dims = 1:pcDims, verbose = isVerbose)
   integrated_obj <- FindClusters(integrated_obj, resolution = c(0.6,0.8,0.9,1,1.2,1.4,1.6,1.8), verbose = isVerbose)
   
   message("** Saving the integrated object")
   saveRDS(integrated_obj, paste0(outname,".rds"), compress = compress_rds)
   
}

