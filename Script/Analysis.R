merged_seu <- readRDS("/media/rna/VAQUITA/snRNAseq_potato/03_ANALYSIS/merged_seu.rds")
###RUN DImension reduction
merged_seu <- NormalizeData(merged_seu)
merged_seu <- FindVariableFeatures(merged_seu)
merged_seu <- ScaleData(merged_seu)
merged_seu <- RunPCA(merged_seu, npcs = 30)

ElbowPlot(merged_seu, ndims = 30)

###find local neighbourhood clusters with just 10 pcs

merged_seu <- RunUMAP(merged_seu, dims = 1:10)
merged_seu <- FindNeighbors(object = merged_seu, dims = 1:10)              
merged_seu <- FindClusters(object = merged_seu, resolution = 0.5)

DimPlot(merged_seu, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
DimPlot(merged_seu, reduction = "umap", split.by = "orig.ident")





####Now we integrate the dataset from two conditions 
merged_seu <- IntegrateLayers(object = merged_seu, method = RPCAIntegration, orig.reduction = "pca", 
                              new.reduction = "integrated.RPCA",
                              verbose = TRUE)
Reductions(merged_seu)

# re-join layers after integration
merged_seu[["RNA"]] <- JoinLayers(merged_seu[["RNA"]])


merged_seu <- FindNeighbors(merged_seu, reduction = "integrated.RPCA", dims = 1:10)
merged_seu <- FindClusters(merged_seu, resolution = 0.3)

merged_seu <- RunUMAP(merged_seu, dims = 1:10, reduction = "integrated.RPCA")

DimPlot(merged_seu, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = TRUE)



##find conserved markers for each cell clusters
Idents(merged_seu) <- "seurat_clusters"
ConservedMarkers_list <-list() 
for (Cell_clusters in levels(Idents(merged_seu))) {
  
  ConservedMarkers_list[[as.character(Cell_clusters)]] <- FindConservedMarkers(merged_seu, ident.1 = Cell_clusters, 
                                                                               grouping.var = "orig.ident", verbose = TRUE)
  
}

##markers test from paper
Markers<- c("Soltu.DM.01G029600.v6.1",
            "Soltu.DM.05G011870.v6.1",
            "Soltu.DM.08G015830.v6.1",
            "Soltu.DM.09G028830.v6.1",
            "Soltu.DM.02G013580.v6.1",
            "Soltu.DM.08G001280.v6.1",
            "Soltu.DM.05G013680.v6.1",
            "Soltu.DM.10G016620.v6.1",
            "Soltu.DM.08G001720.v6.1",
            "Soltu.DM.06G008810.v6.1",
            "Soltu.DM.11G009620.v6.1",
            "Soltu.DM.09G024150.v6.1",
            "Soltu.DM.04G036750.v6.1",
            "Soltu.DM.02G011940.v6.1",
            "Soltu.DM.03G019570.v6.1",
            "Soltu.DM.03G020090.v6.1",
            "Soltu.DM.05G012430.v6.1",
            "Soltu.DM.04G031670.v6.1",
            "Soltu.DM.11G010180.v6.1",
            "Soltu.DM.12G020640.v6.1",
            "Soltu.DM.07G014110.v6.1",
            "Soltu.DM.04G028940.v6.1",
            "Soltu.DM.12G007170.v6.1",
            "Soltu.DM.07G012350.v6.1",
            "Soltu.DM.08G012380.v6.1",
            "Soltu.DM.09G024210.v6.1",
            "Soltu.DM.12G023640.v6.1",
            "Soltu.DM.02G010430.v6.1",
            "Soltu.DM.02G032930.v6.1",
            "Soltu.DM.01G031660.v6.1",
            "Soltu.DM.04G036350.v6.1",
            "Soltu.DM.03G014810.v6.1",
            "Soltu.DM.05G021590.v6.1",
            "Soltu.DM.04G033400.v6.1",
            "Soltu.DM.07G012640.v6.1",
            "Soltu.DM.02G007010.v6.1",
            "Soltu.DM.05G023840.v6.1",
            "Soltu.DM.09G028830.v6.1",
            "Soltu.DM.06G024930.v6.1",
            "Soltu.DM.04G024870.v6.1",
            "Soltu.DM.04G027770.v6.1",
            "Soltu.DM.06G014260.v6.1")


S.tuberosum.v6.1_latest$locusName <- paste0(S.tuberosum.v6.1_latest$locusName,".v6.1")

###based on markers identified from paper
DotPlot(
  merged_seu,
  features = c("Soltu.DM.01G029600.v6.1",
               "Soltu.DM.05G011870.v6.1",
               "Soltu.DM.08G015830.v6.1",
               "Soltu.DM.09G028830.v6.1",
               "Soltu.DM.02G013580.v6.1",
               "Soltu.DM.08G001280.v6.1",
               "Soltu.DM.05G013680.v6.1",
               "Soltu.DM.10G016620.v6.1",
               "Soltu.DM.08G001720.v6.1",
               "Soltu.DM.06G008810.v6.1",
               "Soltu.DM.11G009620.v6.1",
               "Soltu.DM.09G024150.v6.1",
               "Soltu.DM.04G036750.v6.1",
               "Soltu.DM.02G011940.v6.1",
               "Soltu.DM.03G014810.v6.1")
) + RotatedAxis() + coord_flip()



cluster_ids <- c(
  "0" = "Mesophyll",
  "1" = "Epidermis",
  "2" = "Mesophyll",
  "3" = "Mesophyll",
  "4" = "Mesophyll",
  "5" = "unknown1",
  "6" = "Mesophyll",
  "7" = "Epidermis",
  "8" = "Epidermis",
  "9" = "Epidermis",
  "10" = "unknown4",
  "11" = "Guard cells",
  "12" = "unknown2",
  "13" = "unknown3",
  "14" = "Primoridium cells")


##mapped the cell identities to the seurat clusters
merged_seu$celltype <- plyr::mapvalues(
  x = merged_seu$seurat_clusters,
  from = names(cluster_ids),
  to = cluster_ids
)

DimPlot(merged_seu, reduction = "umap", group.by = c("orig.ident", "celltype"))
DimPlot(merged_seu, reduction = "umap", group.by = c("celltype"), )


#### redo the clustering with original pca reduction with just 10 PCs
## to cluster label by cell type
merged_seu <- FindNeighbors(merged_seu, reduction = "pca", dims = 1:10)
merged_seu <- FindClusters(merged_seu, resolution = 0.3)

merged_seu <- RunUMAP(merged_seu, dims = 1:10, reduction = "pca")

DimPlot(merged_seu, reduction = "umap", group.by = c( "celltype"), split.by = "orig.ident", label = TRUE)
