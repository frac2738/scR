

liana_get_freq <- function(liana_res){
   liana_res %>%
      group_by(source, target) %>%
      summarise(freq = n(), .groups = 'keep') %>%
      pivot_wider(id_cols = source,
                  names_from = target,
                  values_from = freq,
                  values_fill = 0) %>%
      arrange(source) %>%
      ungroup() %>%
      as.data.frame() %>%
      column_to_rownames('source') %>%
      as.matrix()
   
}

liana_anno_barplot <- function(x,
                          cell_anno,
                          axis.font.size){
   anno_barplot(x,
                gp = gpar(fill = cell_anno,
                          col = cell_anno,
                          font.size=axis.font.size),
                axis_param = list(gp=gpar(font.size=axis.font.size)),
                title="",
                border = FALSE)
}


liana_heatmap_fixed <- function(mat,
                          font_size = 12,
                          grid_text = FALSE,
                          name = 'Frequency',
                          pallette = c("white", "violetred2"),
                          row_title = "Sender (Cell types)",
                          column_title = "Receiver (Cell types)",
                          ...){
   
   library(rlang)
   library(ComplexHeatmap)
   
   if(grid_text){
      grid_text <- function(j, i, x, y, width, height, fill) {
         grid_text <- grid.text(sprintf("%d", mat[i, j]),
                                x, y, gp = gpar(fontsize = font_size*0.83))
      }
   } else{
      grid_text <- NULL
   }
   
   # define Annotations and Barplots
   cell_anno <- unique(rownames(mat))
   cell_anno <- grDevices::colorRampPalette(
      (RColorBrewer::brewer.pal(n = 8, name = 'Dark2'))
   )(length(cell_anno)) %>%
      setNames(cell_anno)
   
   ## Annotations
   ha_opts <- list(show_legend = FALSE,
                   show_annotation_name = FALSE,
                   col = list("anno"=cell_anno),
                   simple_anno_size = grid::unit(0.25, "cm"))
   column_ha <- exec("HeatmapAnnotation", anno = names(cell_anno), !!!ha_opts)
   row_ha <- exec("rowAnnotation", anno = names(cell_anno), !!!ha_opts)
   
   # Barplots
   column_bar <- ComplexHeatmap::HeatmapAnnotation(
      bar = liana_anno_barplot(colSums(mat),
                          cell_anno,
                          axis.font.size = font_size*0.5
      ),
      annotation_name_gp = gpar(fontsize = font_size*0.5),
      show_legend = FALSE,
      show_annotation_name = FALSE)
   
   row_bar <- ComplexHeatmap::rowAnnotation(
      bar2 = liana_anno_barplot(rowSums(mat),
                           cell_anno,
                           font_size*0.5
      ),
      gp = gpar(fill = cell_anno,
                col = cell_anno),
      show_legend = FALSE,
      show_annotation_name = FALSE)
   
   # Heatmap
   ComplexHeatmap::Heatmap(mat,
                           col=colorRampPalette(pallette)(10),
                           cluster_rows = FALSE,
                           cluster_columns = FALSE,
                           row_names_side = "left",
                           top_annotation = column_bar,
                           #bottom_annotation = column_ha,
                           right_annotation = row_bar,
                           #left_annotation = row_ha,
                           row_title = row_title,
                           row_names_gp = gpar(fontsize = font_size),
                           row_title_gp = gpar(fontsize = font_size*1.2),
                           column_names_gp = gpar(fontsize = font_size),
                           column_title = column_title,
                           column_title_gp = gpar(fontsize = font_size*1.2),
                           column_title_side = "bottom",
                           heatmap_legend_param = list(title_gp = gpar(fontsize = font_size*0.9,
                                                                       fontface = 'bold'),
                                                       border = NA,
                                                       labels_gp = gpar(fontsize = font_size*0.9),
                                                       grid_width = unit(2, "mm")),
                           name = name,
                           cell_fun = grid_text
   )
   
}