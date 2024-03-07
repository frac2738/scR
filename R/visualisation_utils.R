#' Modify scales of a Feature Plot.
#' 
#' This function transforms all the values in a FeaturePlot that are below min.cutoff to min.cutoff and above max.cutoff to max.cutoff, effectively
#' achieving the same behaviour as in Seurat::FeaturePlot. All credits for this idea goes to the Seurat developers.
#' 
#' @param p Plot resulting from calling SCpubr::do_FeaturePlot() with a single feature and no extra parameters.
#' @param feature Feature used in p.
#' @param min.cutoff Minimum range of the scale.
#' @param max.cutoff Maximum range of the scale.
#'
#' @return A ggplot2 object with the scales modified.
#' @export
#'
#' @examples .
modify_scales <- function(p,
                          feature,
                          min.cutoff = NULL,
                          max.cutoff = NULL){
   `%>%` <- magrittr::`%>%`
   `:=` <- rlang::`:=`
   
   # Apply min.cutoff.
   if (!is.null(min.cutoff)){
      p$data <- p$data %>% 
         dplyr::mutate("{feature}" := ifelse(.data[[feature]] < min.cutoff, min.cutoff, .data[[feature]]))
   }
   
   # Apply max.cutoff.
   if (!is.null(max.cutoff)){
      p$data <- p$data %>% 
         dplyr::mutate("{feature}" := ifelse(.data[[feature]] > max.cutoff, max.cutoff, .data[[feature]]))
   }
   
   return(p)
}
