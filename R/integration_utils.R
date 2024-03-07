
sketch_integration <- function(dataDIR, samples_vect = NULL, intMethod, outname, nfeatures = 3000, 
                               pcDims = 30, isVerbose = FALSE, merge_layers = TRUE, seurat5 = TRUE,
                               checkpoint = FALSE){
   
   # dataDIR <- "/mnt/disk1/projects/scATLAS/03_processed_data/full_atlas_v1"
   # samples_vect <- NULL
   # intMethod <- "rpca"
   # nfeatures <- 3000
   # pcDims <- 30
   # isVerbose <- TRUE
   
   library(cli)
   
   options(Seurat.object.assay.version = "v5")
   cli_inform(message = "** Loading samples data: ")
   
   if(is.null(samples_vect)){
      samples_vect <- list.files(dataDIR, full.names = FALSE, recursive = FALSE, pattern = "rds") 
   } 
   
   sample_list <- c()
   metadata_list <- data.frame(matrix(ncol = 10, nrow = 0))
   names(metadata_list) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt",
                             "percent.ribo","log10GenesPerUMI","SampleID","source",
                             "RNA_snn_res.0.8","seurat_clusters")
   
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
      
      # remove celline columns      
      if(length(grep("K562|K562_kNN", names(mySample@meta.data))) >0){
         mySample@meta.data <- mySample@meta.data[,names(metadata_list)]
      }
      
      metadata_list <- rbind(metadata_list,mySample@meta.data)
      sample_list[[sample_name]] <- mySample@assays$RNA$counts
      rm(mySample)
      
      gc()
   }
   
   message("** Create layered object")
   integrated_obj <- CreateSeuratObject(counts = sample_list,
                                        meta.data = metadata_list)
   
   integrated_obj <- NormalizeData(integrated_obj)
   integrated_obj <- FindVariableFeatures(integrated_obj)
   
   sketch_cells <- round(nrow(integrated_obj@meta.data)/4, digits = 0)
   
   integrated_obj <- SketchData(
      object = integrated_obj,
      ncells = sketch_cells,
      method = "LeverageScore",
      sketched.assay = "sketch",
      seed = 666,
   )
   
   DefaultAssay(integrated_obj) <- "sketch"
   integrated_obj <- FindVariableFeatures(integrated_obj)
   integrated_obj <- ScaleData(integrated_obj)
   integrated_obj <- RunPCA(integrated_obj)
   
   ref_samples_idx <- which(Layers(integrated_obj, search = "data") %in% paste0("data.",ref_samples))
   #ref_samples_idx <- which(Layers(integrated_obj, search = 'data') %in% paste0("counts.",ref_samples))
   
   integrated_obj <- IntegrateLayers(integrated_obj, method = RPCAIntegration, 
                                     orig = "pca", new.reduction = "integrated.rpca",
                                     dims = 1:30, k.anchor = 20,
                                     reference = ref_samples_idx)
   
   integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30)
   integrated_obj <- FindClusters(integrated_obj, resolution = 2)
   integrated_obj <- RunUMAP(integrated_obj, dims = 1:30, return.model = T)

   integrated_obj <- ProjectIntegration(integrated_obj, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")

   integrated_obj <- ProjectData(integrated_obj, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                         full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(cluster_full = "seurat_clusters"))
   
   integrated_obj <- RunUMAP(integrated_obj, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
   
   DefaultAssay(integrated_obj) <- "RNA"
   DimPlot(integrated_obj, label = T, label.size = 3, reduction = "umap.full", group.by = "cluster_full") + NoLegend()
   
   return(integrated_obj)

}

