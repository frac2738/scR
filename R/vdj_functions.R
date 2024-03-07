#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Build and visualise phylogenetic trees on AIRR formatted file 
#' 
#' @param germ_file AIRR formatted table
#' @param clonesID List of clones to build the tree for. If \code{NULL} all the clones will be used | Default \code{NULL} 
#' @param method Method used to estimate tree topology and branch lengths. Available methods are: pratchet, pml, dnapars, dnaml and igphyml | Default \code{pratchet} 
#' @param exec Path to dnapars executable (only used if \code{method = dnapars}) | Default \code{NULL} 
#' @param scale With of the scale bar | Default \code{0.01} 
#' @export 
#' @examples
#' \dontrun{
#' bcr_lineage_reconstruction("IGHV-genotyped_M2_germ-pass.tab",c(379,375,455,291),scale = 1)
#' }
bcr_lineage_reconstruction <- function(germ_file,clonesID = NULL,method="pratchet", exec = NULL, scale = 0.01){
   
   suppressPackageStartupMessages(library(ggtree))
   library(dowser)
   
   if(method == "dnapars" && is.null(exec)){
      exec <- "/home/fcarbone/softwares/phylip-3.697/exe/dnapars"
   }
   
   airr_table <- read.delim(germ_file)
   names(airr_table) <- tolower(names(airr_table))
   airr_table$clone_id <- airr_table$clone
   
   if(is.null(clonesID)){
      available_clones <- as.data.frame(table(airr_table$clone_id))
      available_clones <- available_clones[available_clones$Freq>1,]
      available_clones <- available_clones[order(available_clones$Freq, decreasing = TRUE),]
      clonesID <- available_clones$Var1
      
   }
   
   clones_table <- airr_table[which(airr_table$clone_id %in% clonesID),]
   
   clones <- formatClones(clones_table, traits="c_call", seq = "sequence_imgt", germ = "germline_imgt_d_mask")
   trees <- getTrees(clones, build = method, exec = exec)
   trees <- scaleBranches(trees, edge_type="mutations")
   plots <- plotTrees(trees, tips="c_call",tipsize=2, scale = scale)
   
   out_dir <- paste0("Lineage_clones_",paste(as.character(clonesID), collapse = "-"))
   dir.create(out_dir)
   
   for(plottino in 1:length(plots)){
      
      message(paste0("Plotting Clone ",clonesID[plottino]))
      gg_tree <- plots[[plottino]]
      
      tree_plot_name <- paste0(out_dir,"/",clonesID[plottino],".png")
      ggsave(gg_tree, filename= tree_plot_name)
   }
   
}

