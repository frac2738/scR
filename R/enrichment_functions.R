#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Run enrichR
#' 
#' @param genes2enrich .
#' @param databases .
#' @param n_top .
#' @param direction .
#' @param out_dir .
#' @keywords enrichR
#' @export
#' @examples
#' \donttest{
#' 
#'}
enrichR_analysis <- function(genes2enrich, databases= NULL, n_top = NULL, direction = "up", out_dir = "."){
   
   library(enrichR)
   library(RColorBrewer)
   
   if(direction == "up") colori <- rev(brewer.pal(9, "YlOrRd")[1:7]) else colori <- rev(brewer.pal(9, "YlGnBu")[1:7])
   
   if(is.null(databases)){
      enrichR_db <- c("Reactome_2022","MSigDB_Hallmark_2020","WikiPathway_2021_Human",
                      "KEGG_2021_Human","BioCarta_2016") 
   } else {
      enrichR_db <- databases
   }
   
   ora_results <- enrichr(genes2enrich, enrichR_db)
   
   return(ora_results)

}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Generate the enrichment figure for enrichR
#' 
#' @param genes2enrich .
#' @param databases .
#' @param n_top .
#' @param direction .
#' @param out_dir .
#' @keywords enrichR
#' @export
#' @examples
#' \donttest{
#' 
#'}
#' @examples
#' \donttest{
#' 
#'}
plot_enrichR_result <- function(enrichR_table, nPaths, colori = NULL, title = ""){
   
   library(ggplot2)
   library(cli)
   
   nPaths <- 20
   
   if(is.null(colori)){
      colori <- rev(brewer.pal(9, "YlOrRd")[1:7])
   }
   
   if(nrow(enrichR_table)>0){
      
      if(nrow(enrichR_table) < nPaths){
         enrichR_table
         titolo <- title
      } else {
         enrichR_table <- enrichR_table[1:nPaths,]
         titolo <- paste0(title," (top ",nrow(enrichR_table)," paths)")
      }
      
      enrichR_table <- enrichR_table[order(enrichR_table$Combined.Score,decreasing = T),]
      enrichR_table$Signif <- ""
      enrichR_table$Signif[which(enrichR_table$Adjusted.P.value<=0.05)] <- "* "
      enrichR_table$Signif[which(enrichR_table$Adjusted.P.value<=0.01)] <- "**"
      enrichR_table$Signif[which(enrichR_table$Adjusted.P.value<=0.001)] <- "***"
      enrichR_table$Signif[which(enrichR_table$Adjusted.P.value<=0.0001)] <- "****"
      
      ggEnriched <- ggplot(enrichR_table, aes(x = Term, y = Combined.Score, fill = Adjusted.P.value)) +
         geom_bar(stat="identity", width=0.8) + coord_flip() + 
         scale_x_discrete(limits = rev(enrichR_table$Term)) +
         theme_bw() + 
         scale_fill_gradientn(colours = colori) + 
         geom_text(aes(label=Signif),position ="stack", hjust = -0.3) +
         geom_text(aes(label=Overlap), position = position_stack(0.5)) +
         ylim(0, max(enrichR_table$Combined.Score)+0.2*max(enrichR_table$Combined.Score)) +
         xlab("") + ggtitle(titolo)
      
      
   } else {
      cli_alert_warning("The enrichment table provided is empty. Nothing to do!")
      
   }
   
   return(ggEnriched)
   
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#' Rank genes for FGSEA analysis
#' 
#' @param deg_output .
#' @param isHuman .
#' @param entrez .
#' @keywords FGSEA
#' @export
#' @examples
#' \donttest{
#' 
#'}
fgsea_rank_genes <- function(deg_output, isHuman=TRUE, entrez=TRUE){
   
   ranked_genes <- deg_output$avg_log2FC
   names(ranked_genes) <- deg_output$GeneName
   
   if(isHuman){
      source <- org.Hs.eg.db
      
   } else {
      source <- org.Mm.eg.db
   }
   
   if(entrez){
      symbols <- mapIds(source, keys = names(ranked_genes), 
                        keytype = "SYMBOL", column=c("ENTREZID"))
      
      fgsea_input <- ranked_genes
      names(fgsea_input) <- symbols
      
   } else {
      fgsea_input <- ranked_genes
      
   }
   fgsea_input <- fgsea_input[!is.na(names(fgsea_input))]
   
   return(fgsea_input)
   
}
