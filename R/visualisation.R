#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Improved violin plot
#' 
#' @param srtObject Seurat object
#' @param features Vector of features to plot
#' @param group_by Name of one or more metadata columns to group cells
#' @param split_by Name of a metadata column to split plot by
#' @param idents Identity classes to include in plot (default is all)
#' @param invert Invert the Identity selected (idents) 
#' @param assay Name of assay to use, defaults to the active assay
#' @param slot Name of the slot to use
#' @param boxplot Overlay a boxplot
#' @param group_comparison vector of colours to use to colour the feature names
#' @param colors_pallet Color palette to use. viridis by default
#' @param title Figure title
#' @keywords seurat
#' @export
#' @examples
#' \donttest{
#' library(Seurat)
#' 
#' features <- c("MS4A1","CD3E","CD8A")
#' slot <- "data"
#' assay <- "RNA"
#' group_by <- "groups"
#' group_comparison  <- list(c("g1","g2"))
#' VlnPlot_modified(pbmc_small, features, group_by, split_by = NULL, idents = NULL, 
#'                   invert = FALSE, assay, slot, boxplot = TRUE, group_comparison = group_comparison,"a title")
#'}
VlnPlot_modified <- function(srtObject, features, group_by, split_by = NULL, idents = NULL, 
                             invert = FALSE, assay, slot, log = FALSE, boxplot= FALSE, group_comparison = NULL, alpha = 0.6,
                             colors_pallet = "viridis", title, right_labels_col = NULL, dup_x_axis = FALSE,
                             y_scale = "free_y"){
   
   suppressPackageStartupMessages(library(ggplot2))
   suppressPackageStartupMessages(library(ggpubr))
   suppressPackageStartupMessages(library(wesanderson))
   suppressPackageStartupMessages(library(viridis))
   suppressPackageStartupMessages(library(data.table))
   
   DefaultAssay(srtObject) <- assay
   
   if(!is.null(idents)){
      srtObject <- subset(srtObject, idents = idents, invert = invert)
   }
   
   # check if features are in the meta.data or exp matrix
   inMetadata <- ifelse(features %in% names(srtObject@meta.data),TRUE,FALSE)
   
   if(all(inMetadata)){
      plot_data <- srtObject@meta.data[,features, drop=FALSE]
      
   } else {
      plot_data <- FetchData(object = srtObject, vars = features, slot = slot)
   }
   
   #plot_data$group <- srtObject@meta.data$group
   plot_data$group_by <- eval(parse(text = paste0("srtObject@meta.data$",group_by) ))
   
   if(!is.null(split_by)){
      plot_data$split_by <- eval(parse(text = paste0("srtObject@meta.data$",split_by) ))
      
   }
   
   plot_data <- setDT(plot_data)
   melt_by <- names(plot_data)
   plot_data <- melt(plot_data,id.vars = setdiff(names(plot_data), features))
   names(plot_data) <- sub("group_by","group",names(plot_data))
   names(plot_data) <- sub("split_by","split",names(plot_data))
   
   if(log){
      plot_data$value <- log10(plot_data$value+1)
      y_label_text <- paste0("log(",slot,"+1)")
   } else {
      y_label_text <- paste0(assay)
   }
   
   # set color palette based on the value of group_by
   col_length <- length(unique(plot_data$group))
   
   if(length(colors_pallet) >1){
      colors_use <- colors_pallet
   } else {
      
      if(colors_pallet %in% c("Zissou1","BottleRocket2","Rushmore1","Darjeeling1")){
         colors_use <- wes_palette(colors_pallet, col_length, type = "continuous")
         
      } else if(colors_pallet %in% c("magma","inferno","plasma","viridis","cividis","rocket","mako","mako")) {
         colors_use <- viridis(col_length, alpha = 1, begin = 0, end = 1, direction = -1, option = colors_pallet)
      } 
   }
   
   #colors_use <- col_set1
   
   vln_return <- ggplot(plot_data, aes(group, value, fill = group , color = group)) + 
      geom_violin( alpha = alpha, trim = TRUE, draw_quantiles = FALSE, scale = "width") +
      facet_grid(rows = vars(variable), scales = y_scale, space = "free_x") +
      xlab("") +
      ylab(y_label_text) +
      ggtitle(title) +
      scale_fill_manual( values = colors_use, na.value = "grey") +
      scale_color_manual(values = colors_use, na.value = "grey") +
      theme(legend.position = "none",
         axis.ticks.x = element_blank(),
         axis.title.y = element_text(size = rel(1), angle = 90),
         axis.text.y = element_text(size = rel(1)),
         axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
         strip.text.x = element_text(size = 14),
         strip.text.y = element_text(size = 14),
         panel.background = element_rect(fill = "white", colour = "grey50")) 
   
   if(dup_x_axis){
      vln_return + guides(x.sec = guide_axis_label_trans(~ .x))
      
   }
   
   if(boxplot){
      vln_return <- vln_return + geom_boxplot(outlier.shape = NA, width = 0.1, color="black") 
   }
   
   if(!is.null(split_by)){
      vln_return <- vln_return + facet_grid(rows = vars(variable), cols = vars(split) ,scales = y_scale, space = "free_x")
      #vln_return <- vln_return + facet_grid(rows = vars(split), cols = vars(variable) ,scales = y_scale, space = "free_x")
      
   }
   
   if(!is.null(group_comparison)){
      vln_return <- vln_return + stat_compare_means(comparisons = group_comparison, label = "p.signif",  hide.ns = TRUE) 
   }
   
   #vln_return + guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
   #return(vln_return)
   
   if(!is.null(right_labels_col)){
      
      library(tidyverse)
      library(grid)
      library(gridExtra)
      
      g <- ggplot_gtable(ggplot_build(vln_return))
      stripr <- which(grepl('strip-r', g$layout$name))
      fills <- right_labels_col
      k <- 1
      for (i in stripr) {
         j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
         g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
         k <- k+1
      }
      
      message("Figure saved in the working directory")
      pdf(file = "./vlnPlot_grid_arrange.pdf",width = 16, height = 22)
      print(grid.draw(g)) 
      dev.off()
      
      return("./vlnPlot_grid_arrange.pdf")
      
   } else {
      return(vln_return)
      
   }
   
}

#' Improved Dotplot
#'
#' This function is a modified version of the Clustered DotPlots function from the scCustomize package.
#' The changes made are:
#' 
#' * Removed the colomn annotation bar
#' * Added the option of clustering the identities (rows)
#'
#' @param seurat_object Seurat object name.
#' @param features Features to plot.
#' @param colors_use_exp Color palette to use for plotting expression scale.  Default is `viridis::plasma(n = 20, direction = -1)`.
#' @param exp_color_min Minimum scaled average expression threshold (everything smaller will be set to this).
#' Default is -2.
#' @param exp_color_middle What scaled expression value to use for the middle of the provided `colors_use_exp`.
#' By default will be set to value in middle of `exp_color_min` and `exp_color_max`.
#' @param exp_color_max Minimum scaled average expression threshold (everything smaller will be set to this).
#' Default is 2.
#' @param print_exp_quantiles Whether to print the quantiles of expression data in addition to plots.
#' Default is FALSE.  NOTE: These values will be altered by choices of `exp_color_min` and `exp_color_min`
#' if there are values below or above those cutoffs, respectively.
#' @param colors_use_idents specify color palette to used for identity labels.  By default if
#' number of levels plotted is less than or equal to 36 it will use "polychrome" and if greater than 36
#' will use "varibow" with shuffle = TRUE both from `DiscretePalette_scCustomize`.
#' @param x_lab_rotate How to rotate column labels.  By default set to `TRUE` which rotates labels 45 degrees.
#' If set `FALSE` rotation is set to 0 degrees.  Users can also supply custom angle for text rotation.
#' @param flip logical, whether to flip the axes of final plot.  Default is FALSE; rows = features and
#' columns = idents.It coul dbe because of a lot of You should add some details (e.g. figures and code used).
#' @param k Value to use for k-means clustering on features  Sets (km) parameter in `ComplexHeatmap::Heatmap()`.
#' From `ComplexHeatmap::Heatmap()`: Apply k-means clustering on rows. If the value is larger than 1, the
#' heatmap will be split by rows according to the k-means clustering. For each row slice, hierarchical
#' clustering is still applied with parameters above.
#' @param feature_km_repeats Number of k-means runs to get a consensus k-means clustering for features.
#' Note if `feature_km_repeats` is set to value greater than one, the final number of groups might be
#' smaller than row_km, but this might mean the original row_km is not a good choice.  Default is 1000.
#' @param row_km_repeats `r lifecycle::badge("deprecated")` soft-deprecated.  See `feature_km_repeats`
#' @param ident_km_repeats Number of k-means runs to get a consensus k-means clustering. Similar to
#' `feature_km_repeats`.  Default is 1000.
#' @param column_km_repeats `r lifecycle::badge("deprecated")` soft-deprecated.  See `ident_km_repeats`
#' @param row_label_size Size of the feature labels.  Provided to `row_names_gp` in Heatmap call.
#' @param raster Logical, whether to render in raster format (faster plotting, smaller files).  Default is FALSE.
#' @param plot_km_elbow Logical, whether or not to return the Sum Squared Error Elbow Plot for k-means clustering.
#' Estimating elbow of this plot is one way to determine "optimal" value for `k`.
#' Based on: \url{https://stackoverflow.com/a/15376462/15568251}.
#' @param elbow_kmax The maximum value of k to use for `plot_km_elbow`.  Suggest setting larger value so the
#' true shape of plot can be observed.  Value must be 1 less than number of features provided.  If NULL parameter
#' will be set dependent on length of feature list up to `elbow_kmax = 20`.
#' @param assay Name of assay to use, defaults to the active assay.
#' @param group.by Group (color) cells in different ways (for example, orig.ident).
#' @param idents Which classes to include in the plot (default is all).
#' @param show_parent_dend_line Logical, Sets parameter of same name in `ComplexHeatmap::Heatmap()`.
#' From `ComplexHeatmap::Heatmap()`: When heatmap is split, whether to add a dashed line to mark parent
#' dendrogram and children dendrograms.  Default is TRUE.
#' @param ggplot_default_colors logical.  If `colors_use = NULL`, Whether or not to return plot using
#' default ggplot2 "hue" palette instead of default "polychrome" or "varibow" palettes.
#' @param color_seed random seed for the "varibow" palette shuffle if `colors_use = NULL` and number of
#' groups plotted is greater than 36.  Default = 123.
#' @param seed Sets seed for reproducible plotting (ComplexHeatmap plot).
#' @param col_dendr Logical, cluster the genes (default is FALSE).
#' @param row_dendr Logical, cluster the clusters (default is FALSE).
#' @param figure_title Figure title column
#' @return A ComplexHeatmap or if plot_km_elbow = TRUE a list containing ggplot2 object and ComplexHeatmap.
#'
#' @import cli
#' @import ggplot2
#' @importFrom circlize colorRamp2
#' @importFrom dplyr any_of filter select
#' @importFrom grid grid.circle grid.rect gpar
#' @importFrom magrittr "%>%"
#' @importFrom Seurat DotPlot
#' @importFrom SeuratObject PackageCheck
#' @importFrom stats quantile
#' @importFrom tidyr pivot_wider
#' 
#' @export
#'
#' @concept seurat_plotting
#'
#' @author Ming Tang (Original Code), Sam Marsh (Wrap single function, added/modified functionality) Francesco Carbone (Final tweakings)
#' @references \url{https://divingintogeneticsandgenomics.rbind.io/post/clustered-dotplot-for-single-cell-rnaseq/}
#'
#' @examples
#' \donttest{
#' library(Seurat)
#' Clustered_DotPlot_nl(seurat_object = pbmc_small, features = c("CD3E", "CD8", "GZMB", "MS4A1"),kcut = 2)
#'}
#'
#'
Clustered_DotPlot_nl <- function(
   seurat_object,
   features,
   scale_dot = TRUE,
   colors_use_exp = c("#0000FF","grey","#FF4C4C"),
   exp_color_min = -2,
   exp_color_middle = NULL,
   exp_color_max = 2,
   print_exp_quantiles = FALSE,
   colors_use_idents = NULL,
   x_lab_rotate = TRUE,
   kcut = 1,
   border = "black",
   row_km_repeats = 1000,
   column_km_repeats = 1000,
   row_label_size = 8,
   col_label_size = 8,
   raster = FALSE,
   plot_km_elbow = TRUE,
   elbow_kmax = NULL,
   assay = NULL,
   group.by = NULL,
   idents = NULL,
   show_parent_dend_line = TRUE,
   ggplot_default_colors = FALSE,
   seed = 666,
   col_dendr = TRUE,
   row_dendr = TRUE,
   figure_title = ""
) {
   
   library(scCustomize)
   library(ComplexHeatmap)
   library(circlize)
   library(tidyr)
   library(dplyr)
   
   # Check for packages
   ComplexHeatmap_check <- PackageCheck("ComplexHeatmap", error = FALSE)
   if (!ComplexHeatmap_check[1]) {
      stop(
         "Please install the ComplexHeatmap package to use Clustered_DotPlot_nl",
         call. = FALSE
      )
   }
   
   # Check Seurat
   #Is_Seurat(seurat_object = seurat_object)
   
   # Check unique features
   features_unique <- unique(x = features)
   
   if (length(x = features_unique) != length(x = features)) {
      warning("Feature list contains duplicates, making unique.")
   }
   
   # Check exp min/max set correctly
   if (!exp_color_min < exp_color_max) {
      stop("The value for 'exp_color_min': ", exp_color_min, ", must be less than the value for 'exp_color_max': ", exp_color_max, ".")
   }
   
   # Get DotPlot data
   seurat_plot <- DotPlot(object = seurat_object, features = features_unique, assay = assay, group.by = group.by, 
                          scale = scale_dot, idents = idents, col.min = NULL, col.max = NULL)
   #seurat_plot <- DotPlot(object = seurat_object, features = features_unique, scale = TRUE,  col.min = NULL, col.max = NULL)
   
   data <- seurat_plot$data
   
   # remove NaN
   drop <- which(is.na(data$avg.exp.scaled))
   if(length(drop)>0){
      data <- data[-drop,]
   }
   
   # Get expression data
   exp_mat <- data %>%
      select(-pct.exp, -avg.exp) %>%
      pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
      as.data.frame()
   
   drop <- which(is.na(exp_mat$features.plot))
   if(length(drop)>0){
      exp_mat <- exp_mat[-drop,]
   }
   
   row.names(x = exp_mat)  <- exp_mat$features.plot
   exp_mat <- as.matrix(exp_mat[,-1])
   
   excluded_features <- exp_mat[which(rowSums(is.na(x = exp_mat)) > 0),, drop=FALSE] %>% rownames()
   
   # Extract good features
   good_features <- rownames(exp_mat)[-which(rownames(exp_mat) %in% excluded_features)]
   #    
   # # Remove rows with NAs
   # exp_mat <- exp_mat[which(rownames(exp_mat) %in% good_features),]
   
   exp_mat <- t(exp_mat)
   
   # Get percent expressed data
   percent_mat <- data %>%
      select(-avg.exp, -avg.exp.scaled) %>%
      pivot_wider(names_from = id, values_from = pct.exp) %>%
      as.data.frame()
   
   drop <- which(is.na(percent_mat$features.plot))
   if(length(drop)>0){
      percent_mat <- percent_mat[-drop,]
   }
   row.names(x = percent_mat) <- percent_mat$features.plot
   
   # Subset dataframe for NAs if idents so that exp_mat and percent_mat match
   #percent_mat <- percent_mat[which(rownames(percent_mat) %in% good_features),]
   percent_mat <- as.matrix(percent_mat[,-1])
   percent_mat <- t(percent_mat)
   
   # set assay (if null set to active assay)
   assay <- assay %||% DefaultAssay(object = seurat_object)
   
   # Set middle of color scale if not specified
   if (is.null(x = exp_color_middle)) {
      exp_color_middle <- Middle_Number(min = exp_color_min, max = exp_color_max)
   }
   palette_length <- length(colors_use_exp)
   palette_middle <- Middle_Number(min = 0, max = palette_length)
   
   # Set default color palpercent_matette based on number of levels being plotted
   if (is.null(x = group.by)) {
      group_by_length <- length(x = unique(x = seurat_object@active.ident))
   } else {
      group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
   } 
   
   cell_fun <- function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = NA, fill = NA))
      grid.circle(x=x,y=y,r= sqrt(percent_mat[i, j]/100) * unit(2, "mm"),
                  gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
   }
   
   # Create palette
   col_fun = colorRamp2(c(exp_color_min, exp_color_middle, exp_color_max), colors_use_exp[c(1,palette_middle, palette_length)])
   
   # Create legend for point size
   lgd_list = list(
      ComplexHeatmap::Legend(labels = c(0.25,0.5,0.75,1), title = "Percent Expressing",
                             graphics = list(
                                function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(2, "mm"),
                                                                 gp = gpar(fill = "black")),
                                function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(2, "mm"),
                                                                 gp = gpar(fill = "black")),
                                function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(2, "mm"),
                                                                 gp = gpar(fill = "black")),
                                function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                                                 gp = gpar(fill = "black")))
      )
   )
   
   # Set x label roration
   if(is.numeric(x = x_lab_rotate)) {
      x_lab_rotate <- x_lab_rotate
   } else if (isTRUE(x = x_lab_rotate)) {
      x_lab_rotate <- 60
   } else {
      x_lab_rotate <- 0
   }
   
   # Create Plot
   set.seed(seed = seed)
   exp_mat[is.na(exp_mat)] <- 0
   cluster_dot_plot <- ComplexHeatmap::Heatmap(exp_mat,
                                               col=colors_use_exp,
                                               rect_gp = gpar(type = "none"),
                                               cell_fun = cell_fun,
                                               row_names_gp = gpar(fontsize = row_label_size),
                                               column_names_gp = gpar(fontsize = col_label_size),
                                               row_km = kcut,
                                               row_km_repeats = row_km_repeats,
                                               border = border,
                                               cluster_columns = col_dendr, 
                                               cluster_rows = row_dendr,
                                               column_km_repeats = column_km_repeats,
                                               show_parent_dend_line = show_parent_dend_line,
                                               column_names_rot = x_lab_rotate,
                                               column_title = figure_title)
   cluster_dot_plot
   return(ComplexHeatmap::draw(cluster_dot_plot, annotation_legend_list = lgd_list))
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Calculate cell proportions
#'
#' This function calculate the cell proportions.
#' 
#' @param srtObject Seurat object.
#' @param cluster_col Metadata column containing the clustering information.
#' @param groups_annotation data.frame containing sample information. 
#' @param mycolors something
#' @param plotType something
#' @param plotTitle  something
#' @param sampleOrder  something
#' @param groupOrder something
#' @param noClusters something
#' @param stats_compare something
#' @param png something
#' @keywords seurat
#' @export
#' @examples
#' \dontrun{calculate_cell_proportions()}
#' 
calculate_cell_proportions <- function(srtObject, cluster_col, groups_annotation, 
                                       mycolors, plotType, plotTitle,
                                       sampleOrder=NULL, groupOrder = NULL,
                                       noClusters= FALSE, stats_compare = FALSE, png=TRUE){
   
   old <- theme_set(theme_bw())
   
   if(png){
      out_format <- "png"
   } else {
      out_format <- "pdf"
   }
   
   metadata <- srtObject@meta.data
   #metadata$source <- metadata$orig.ident
   
   cluster_col_val <- eval(parse(text = paste0("metadata$`",cluster_col,"`")))
   defined_cellTypes <- unique(metadata$cellType_final)
   
   mappedCells <- list()
   for(df_cells in defined_cellTypes){
      tmp_df <- metadata[which(metadata$cellType_final == df_cells),]
      tmp_clusters <- unique(eval(parse(text = paste0("tmp_df$`",cluster_col,"`"))))
      mappedCells[[df_cells]] <- tmp_clusters
      
   }
   
   # add cell-cluster annotation
   if(!noClusters){
      
      metadata$cellType_wclust <- sapply(metadata$cellType_final, 
                                         function(t){
                                            list_idx <- which(sapply(names(mappedCells), FUN=function(X) t %in% X))
                                            tmp_values <- paste0(mappedCells[[list_idx]],collapse = "-")
                                            tmp_name <- paste0(names(list_idx)," [",tmp_values,"]")
                                            return(tmp_name)
                                         } 
      )
      
      metadata$cluster_wType <- paste0(cluster_col_val," [",metadata$cellType_final,"]")
      metadata$cellType_wclust <- factor(metadata$cellType_wclust, 
                                         levels = unique(metadata$cellType_wclust)[order(unique(metadata$cellType_wclust),decreasing = FALSE)])
      metadata$cluster_wType <- factor(metadata$cluster_wType, 
                                       levels = unique(metadata$cluster_wType)[order(unique(metadata$cluster_wType),decreasing = FALSE)])
      
      
   }
   
   # metadata$cluster
   metadata$clusters <- cluster_col_val
   
   finalTable <- data.frame()
   clusterTable <- data.frame()
   
   if(noClusters){
      for(mySample in unique(metadata$SampleID)){
         tmp_meta <- metadata[which(metadata$SampleID == mySample ),]
         
         tmp_meta_cell <- as.data.frame(table(tmp_meta$cellType_final)/sum(table(tmp_meta$cellType_final))*100)
         tmp_meta_cell$sample <- mySample
         finalTable <- rbind(finalTable,tmp_meta_cell)
         
         tmp_meta_clust <- as.data.frame(table(tmp_meta$clusters)/sum(table(tmp_meta$clusters))*100)
         tmp_meta_clust$sample <- mySample
         clusterTable <- rbind(clusterTable,tmp_meta_clust)
         
      }   
      
   } else {
      
      for(mySample in unique(metadata$SampleID)){
         tmp_meta <- metadata[which(metadata$SampleID == mySample ),]
         
         tmp_meta_cell <- as.data.frame(table(tmp_meta$cellType_wclust)/sum(table(tmp_meta$cellType_wclust))*100)
         tmp_meta_cell$sample <- mySample
         finalTable <- rbind(finalTable,tmp_meta_cell)
         
         tmp_meta_clust <- as.data.frame(table(tmp_meta$cluster_wType)/sum(table(tmp_meta$cluster_wType))*100)
         tmp_meta_clust$sample <- mySample
         clusterTable <- rbind(clusterTable,tmp_meta_clust)
         
      }   
   }
   
   names(finalTable) <- c("cellType","percentage","Sample")
   #finalTable$Sample <- factor(finalTable$Sample)
   finalTable$group <- sapply(finalTable$Sample, function(ita) groups_annotation[match(ita,groups_annotation$SampleID),"group"])
   
   names(clusterTable) <- c("cluster","percentage","Sample")
   #clusterTable$Sample <- factor(clusterTable$Sample)
   clusterTable$group <- sapply(clusterTable$Sample, function(ita) groups_annotation$group[match(ita,groups_annotation$SampleID)])
   
   if(!is.null(sampleOrder)){
      clusterTable$Sample <- factor(clusterTable$Sample, levels = sampleOrder)
      finalTable$Sample <- factor(finalTable$Sample, levels = sampleOrder)
   }
   
   if(!is.null(groupOrder)){
      clusterTable$group <- factor(clusterTable$group, levels = groupOrder)
      finalTable$group <- factor(finalTable$group, levels = groupOrder)
   }
   
   finalTable$percentage <- round(finalTable$percentage,digits = 1)
   clusterTable$percentage <- round(clusterTable$percentage,digits = 1)
   
   write.csv(finalTable,"./cell_proportions_byCellType_table.csv",row.names = FALSE)
   write.csv(clusterTable,"./cell_proportions_byCluster_table.csv",row.names = FALSE)
   
   if(plotType == "boxplot"){
      cb_byCellType_box <- ggplot(finalTable, aes(group,percentage,fill=group)) + 
         geom_boxplot(outlier.shape = NA) +
         facet_wrap("cellType",scales = "free_y") +
         geom_jitter() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ggtitle(plotTitle) +
         scale_fill_manual(values=mycolors) 
      
      cb_byCluster_box <- ggplot(clusterTable, aes(group,percentage,fill=group)) + 
         geom_boxplot(outlier.shape = NA) +
         facet_wrap("cluster",scales = "free_y") +
         geom_jitter() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ggtitle(plotTitle) +
         scale_fill_manual(values=mycolors) 
      
      outname_boxplot_celltypes <- paste0("./cell_proportions_byCellType.",out_format)
      outname_boxplot_clusters <- paste0("./cell_proportions_byCluster.",out_format)
      
      if(stats_compare){
         
         if(exists("groupsComparisons")){
            
            cb_byCellType_box <- cb_byCellType_box + 
               stat_compare_means(size = 4, test = "wilcox.test", label = "p.signif",
                                  hide.ns = TRUE , comparisons = groupsComparisons) 
            
            cb_byCluster_box <- cb_byCluster_box + 
               stat_compare_means(size = 4, test = "wilcox.test", label = "p.signif",
                                  hide.ns = TRUE , comparisons = groupsComparisons) 
            
            outname_boxplot_celltypes <- paste0("./cell_proportions_byCellType_wStats.",out_format)
            outname_boxplot_clusters <- paste0("./cell_proportions_byCluster_wStats.",out_format)
            
         } else {
            
            message("The groups for the wilcox.test are not defined. Statistical test WILL NOT BE PERFORMED.
Please define 'groupsComparisons' as a list containing the comparisons to test.
               e.g. groupsComparisons <- list(c('Healthy','AGS'))")
         }
         
      }
      
      ggsave(plot = cb_byCellType_box, filename = outname_boxplot_celltypes,
             dpi=300, width=22, height = 18)
      
      ggsave(plot = cb_byCluster_box, filename = outname_boxplot_clusters,
             dpi=300, width=22, height = 18)
   }
   
   if(plotType == "barplot"){
      
      cb_byCellType_bar <- ggplot(finalTable, aes(Sample,percentage,fill=group)) + 
         geom_bar(stat="identity") +
         facet_wrap(~cellType,scales = "free_y") +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ggtitle(plotTitle) +
         scale_fill_manual(values=mycolors) 
      
      ggsave(plot = cb_byCellType_bar, 
             filename = paste0("./cell_proportions_byCellType_barplot.",out_format), dpi=300, width=22, height = 18)
      
      cb_byCluster_bar <- ggplot(clusterTable, aes(Sample,percentage,fill=group)) + 
         geom_bar(stat="identity") +
         facet_wrap("cluster",scales = "free_y") +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
         ggtitle(plotTitle) +
         scale_fill_manual(values=mycolors) 
      
      ggsave(plot = cb_byCluster_bar, 
             filename = paste0("./cell_proportions_byCluster_barplot.",out_format), dpi=300, width=22, height = 18)
      
   }
   Clustered_DotPlot_nl
   theme_set(old)
   
}


#'Plot a STEMNET object by Circular A Posteriori Projection
#'
#'@param x An object of class \code{\link{stemnet}}
#'@param arrangement The order of developmental endpoints on the CAP plot. If \code{NULL} (default), the optimal arrangement is determined using a travelling salesman algorithm, see below.
#'@param top If \code{arrangement} is \code{NULL}, show the top n solutions for arranging the developmental endpoints. Defaults to 1, i.e. show only the best solution.
#'@param color Specify the variable color-coded on the plot. The following options exist:
#'\itemize{
#' \item One of "direction" or "degree": Highlight the degree or direction of lineage priming.
#' \item A numeric vector of length equal the number of cells included in the STEMNET plot
#' \item A boolean vector to impose a density plot on the STEMNET plot
#'}
#'@return Invisibly returns a list of two elements: First, the data underlying the plot; Second, the arrangement of endpoints on the cirlce.
#'@details The posterior probability matrix computed by STEMNET places each cell into a n-simplex with as many corners as there a developmental endpoints. To reduce this high-dimensional space to 2 dimensions for plotting, the endpoints are arranged on the edge of a circle and all cells are placed in between.
#'@details In brief, each endpoint is assigned with an angle \eqn{\alpha_k}. The class probabilties \eqn{p_{ik}} are then transformed to cartesian coordinates by \deqn{x_i = \sum_k p_{ik} \cos (\alpha_k) } and  \deqn{y_i = \sum_k p_{ik} \sin (\alpha_k) }
#'@details To find the optimal arrangement of the developmental endpoints on the circle, lineages with common precursor stages are placed next to each other. The proximity between lineages l and k computed by \deqn{D_{kl} = \sum_i p_{il} * p_{ik}}
#'@details The arrangement with the highest proximity is chosen.
#'@references D. Jaitin, E. Kenigsberg, et al.: Massively parallel single-cell RNA-seq for marker-free decomposition of tissues into cell types. *Science* 343: 776-779
#'@export
stemnet_plot <- function(x, arrangement = NULL, top = 1,  color = "direction") {
   
   x <- stemnet_result_cb
   x <- stemnet_result
   center <- apply(x@posteriors[is.na(x@populations),],2,mean)
   p2norm <- t(apply(x@posteriors,1,function(x) {out <- x/center; out/sum(out)}))
   
   uik <- as.matrix(p2norm)
   
   #if no arrangement is given, determine proximity of lineages
   if (is.null(arrangement)) {
      
      bipotent <- matrix(0,ncol(p2norm),ncol(p2norm),dimnames = list(colnames(p2norm),colnames(p2norm)))
      for (i in colnames(p2norm)) {
         for (j in colnames(p2norm)) {
            bipotent[i,j] <- sum(p2norm[,i] * p2norm[,j])
         }
      }
      
      #now identify the optimal arrangement
      getLength <- function(d, permutation, l = 0) {
         if (length(permutation) < 2) return(l) else {
            current <- permutation[1]
            n <- permutation[2]
            l <- getLength(d, permutation[-1], l = l + d[current, n])
            return(l)
         }
         
      }
      
      all_tours <- combinat::permn(2:ncol(uik))
      #all_tours <- combinat::permn(2:9)
      lengths <- sapply(all_tours, function(x) getLength(bipotent, c(1,x,1)))
      optimals <- all_tours[order(lengths, decreasing = T)]
      arrangement <- lapply(optimals[1:top], function(x) c(1,x))
   }
   if(!is.list(arrangement)) arrangement <- list(arrangement)
   
   #for each arrangement, create plot.
   plots <- list()
   CAPframe <- list()
   display_legend = length(arrangement) == 1
   
   for (k in 1:length(arrangement)) {
      
      optimal <- arrangement[[k]]
      #and determine the positions on the unit circle.
      alpha <- seq(0, 2* pi, length.out= ncol(uik)+1)[1:ncol(uik)]
      names(alpha) <- colnames(p2norm)[optimal]
      alpha <- alpha[colnames(p2norm)]
      #for each cell, determine x & y position
      xx <-c(); yy<-c();
      for (i in 1:nrow(uik)){
         xx <- c(xx, sum(uik[i,] * cos(alpha)))
         yy <- c(yy, sum(uik[i,] * sin(alpha)))
      }
      
      CAPframe[[k]] <- data.frame(x = xx, y =yy, direction = primingDirection(x), class = is.na(x@populations), amount = primingDegree(x))
      
      if (length(color) ==1 ){
         if (color == "direction") {
            plots[[k]] <- ggplot2::qplot(x = x, y=y, color=direction, shape = as.factor(class), data = CAPframe[[k]]) + ggplot2::theme_minimal()  + ggplot2::scale_color_discrete(name="Direction d\nof priming", guide = ifelse(display_legend, "legend",F)) + ggplot2::scale_shape_manual(name = "", labels = c("TRUE" = "Stem Cell", "FALSE" = "Endpoint"), values = c("TRUE" = 18, "FALSE" =  17),guide=ifelse(display_legend, "legend",F))
         } else if (color == "degree") {
            plots[[k]] <- ggplot2::qplot(x = x, y=y, color=log10(amount), shape = as.factor(class), data = CAPframe[[k]])   + ggplot2::theme_minimal() + ggplot2::scale_color_gradientn(colours = c("grey","black","blue","red"),name = "Degree\nof priming",guide=ifelse(display_legend, "colorbar",F))+ ggplot2::theme_minimal() + ggplot2::scale_shape_manual(name = "", labels = c("TRUE" = "Endpoint", "FALSE" = "Stem Cell"), values = c("TRUE" = 17, "FALSE" =  18),guide=ifelse(display_legend, "legend",F))
         } else {
            stop("Invalid string constant supplied to color, choose one of direction or degree")
         }
      } else {
         if (is.logical(color)) {
            plots[[k]] <-ggplot2::ggplot(ggplot2::aes(x = x, y=y),  data = CAPframe[[k]][color,])+   ggplot2::theme_minimal() + ggplot2::geom_point(color="black",alpha=0.1,size=0.1,data=CAPframe[[k]]) + ggplot2::stat_density2d(fill = "green",color=NA, bins = 10,size=0.5,alpha=0.2, geom="polygon")  + ggplot2::coord_cartesian(xlim=c(-1,1),ylim=c(-1,1))
         } else if (is.numeric(color)) {
            plots[[k]] <- ggplot2::qplot(x = x, y=y, color=color, shape = as.factor(class), data = CAPframe[[k]])   + ggplot2::theme_minimal() + ggplot2::scale_color_gradientn(colours = c("grey","black","blue","red"),name = "color",guide=ifelse(display_legend, "colorbar",F))+ ggplot2::theme_minimal() + ggplot2::scale_shape_manual(name = "", labels = c("TRUE" = "Endpoint", "FALSE" = "Stem Cell"), values = c("TRUE" = 17, "FALSE" =  18),guide=ifelse(display_legend, "legend",F))
         }
      }
      
   }
   do.call(gridExtra::grid.arrange, plots)
   if (length(arrangement) == 1) invisible(list ( PlotData = CAPframe[[1]], Arrangement = arrangement[[1]])) else invisible(list ( PlotData = CAPframe, Arrangement = arrangement))
}
