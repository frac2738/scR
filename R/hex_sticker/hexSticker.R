
#--- generate a hexSticker

library(Seurat)
library(ggplot2)
library(hexSticker)
library(cowplot)
library(viridis)
library(wesanderson)

ciao <- readRDS("/mnt/disk1/projects/big-integration/rnaSeq_v2/rds_data/atlas_v5-bis.rds")

# Old sticker
# p1 <- DimPlot(ciao, label = FALSE) + theme_void()+ NoLegend() +
#    theme(axis.text.x=element_blank(),
#          axis.ticks.x=element_blank(),
#          axis.title.x=element_blank(),
#          axis.text.y=element_blank(),
#          axis.ticks.y=element_blank(),
#          axis.title.y=element_blank(),
#          panel.border = element_blank(),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.background = element_rect(fill = '#160042', color = '#160042'),
#          plot.background = element_rect(color = "#160042"))
# p1
# ggsave(p1, filename = "/mnt/disk1/sc-pipelines/scR/man/figures/uMAP_reclustered.png", dpi = 300, height = 10, width = 14, bg = "transparent")
# 
# rm(ciao)
# gc()
# imgurl <- "/mnt/disk1/sc-pipelines/scR/man/figures/uMAP_reclustered.png"
# sticker(imgurl, package="scR", s_x=1, s_y=0.7, s_width=0.5, s_height = 0.4,
#         p_size=40,
#         h_fill= "#160042", h_color = "#BBE116", white_around_sticker = FALSE,dpi = 300, spotlight=FALSE,
#         filename="/mnt/disk1/sc-pipelines/scR/man/figures/scR_logo.png")
# 
# 


umap_3d <- FetchData(ciao, vars = c("UMAP_1", "UMAP_2", "seurat_clusters"))

zissou_palette <- wes_palette("Zissou1", 56, type = "continuous")

p2 <- ggplot(umap_3d, aes(UMAP_1,UMAP_2, colour= seurat_clusters)) + geom_point() + 
   scale_color_manual(values = as.character(zissou_palette)) +
   theme(legend.position="none",
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.title.y=element_blank(),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = '#ffffff', color = '#ffffff'),
         plot.background = element_rect(color = "#ffffff"))
p2 + theme_transparent() + theme_void() 
p2

ggsave(p2, filename = "/mnt/disk1/code_repository/scR/man/figures/atlas_umap.png", dpi = 300, height = 10, width = 14, bg = "transparent")
imgurl <- "/mnt/disk1/code_repository/scR/man/figures/atlas_umap.png"
sticker(imgurl, package="scR", s_x=1, s_y=0.7, s_width=0.6, s_height = 0.6,
        p_size=40, p_color = "#ff9a00",
        h_fill= "#ffffff", h_color = "#ff9a00", white_around_sticker = FALSE, dpi = 300, spotlight=FALSE,
        filename="/mnt/disk1/code_repository/scR/man/figures/scR_logo_v4.png")

# #ff9a00 = light
# #2C002C = dark