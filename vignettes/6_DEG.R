
library(DiagrammeR)

#---- 1

DiagrammeR::grViz("digraph flowchart {

      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']

      # edge definitions with the node IDs
      tab1 -> tab2;
      tab1 -> tab3;
      tab1 -> tab4;
      tab2 -> tab5; 
      tab3 -> tab5;
      tab4 -> tab6;
      tab6 -> tab7;
      tab4 -> tab8;
      }

      [1]: 'Seurat Object'
      [2]: 'log normalisation'
      [3]: 'sctransform (v1)'
      [4]: 'sctransform (v2)'
      [5]: 'RNA'
      [6]: 'PrepSCT'
      [7]: 'SCT'
      [8]: 'RNA'
      
      ")



VizDimLoadings(epithelials, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)


ProjectDim(epithelials,  reduction = "umap", dims = 1:2)
JackStraw(epithelials,  reduction = "pca", dims = 1:2)

VlnPlot(epithelials,  reduction = "pca",features = "PCA_1")



DefaultAssay(epithelials) <- "integrated"
DimHeatmap(epithelials, reduction = "umap", dims = 1, balanced = TRUE)

