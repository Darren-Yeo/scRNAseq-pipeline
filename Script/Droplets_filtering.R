
## this script will analyse the snRNAseq data from
##Wen, S., Wang, K., Liang, W., Liu, R., Li, Z., Chen, X., ... & Jian, H. (2026). 
#Transcription Profiling of Potato Leaves in Response to Heat Stress at Single‐Cell Resolution. Plant Biotechnology Journal.
BiocManager::install("DropletUtils")

library(Seurat)
library(SeuratData)
library(patchwork)
library(Matrix)
library(DropletUtils)
library(tidyverse)
# Read matrix into R
path="/media/rna/VAQUITA/snRNAseq_potato/02_KB_count/counts_unfiltered/"

read_matrix <- ReadMtx(mtx = paste0(path,"cells_x_genes.mtx"),
                       features = paste0(path,"cells_x_genes.genes.txt"),
                       cells = paste0(path,"cells_x_genes.barcodes.txt"),
                       feature.column = 1) #because no gene name, cell ranger has two columns ID, gene_name
?ReadMtx
mtx <- readMM(paste0(path,"cells_x_genes.mtx"))
dim(mtx)
?read_count_output

read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(path, "/", name, ge)))
  barcodes <- readLines(file(paste0(path, "/", name, ".prefix.barcodes_together.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}


res_mat <- read_count_output(paste0(path), name = "cells_x_genes")
as.matrix(res_mat[1:10, 1:10])
as.matrix(res_mat[(nrow(res_mat)-9):nrow(res_mat), (ncol(res_mat)-9):ncol(res_mat)])

cell_barcodes<- colnames(res_mat)

##summary of the matrix
tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

##separate cols by sample identifier

AAAAAAAAAAAAAAAC
AAAAAAAAAAAAAAAG
AAAAAAAAAAAAAAAT
AAAAAAAAAAAAAACA
AAAAAAAAAAAAAACC

res_mat_Sep_list<- list(
  CK1 = res_mat[, grepl("AAAAAAAAAAAAAAAA", colnames(res_mat))],
  CK2 = res_mat[, grepl("AAAAAAAAAAAAAAAC", colnames(res_mat))],
  CK3 = res_mat[, grepl("AAAAAAAAAAAAAAAG", colnames(res_mat))],
  HS1 = res_mat[, grepl("AAAAAAAAAAAAAAAT", colnames(res_mat))],
  HS2  = res_mat[, grepl("AAAAAAAAAAAAAACA", colnames(res_mat))],
  HS3  = res_mat[, grepl("AAAAAAAAAAAAAACC", colnames(res_mat))])

bc_rank_list<- lapply(res_mat_Sep_list, barcodeRanks, lower=100)

bc_rank <- barcodeRanks(res_mat, lower = 10)

knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}

##kneeplot of the droplets
gridExtra::grid.arrange(
  knee_plot(bc_rank_list$CK1),
  knee_plot(bc_rank_list$CK2),
  knee_plot(bc_rank_list$CK3),
  knee_plot(bc_rank_list$HS1),
  knee_plot(bc_rank_list$HS2),
  knee_plot(bc_rank_list$HS3)
)

run_emptydrops <- function(mat, fdr_cutoff =  0.001) {
  
  # EmptyDrops expects genes x cells OR cells x genes depending on version
  # kb-count is usually genes x cells, so we keep as-is
  
  ed <- emptyDrops(mat, lower = 500)
  
  # keep only significant barcodes
  keep <- which(ed$FDR < fdr_cutoff)
  
  mat[, keep, drop = FALSE]
}

s <- "CK1"

mat <- res_mat_Sep_list[[s]]

# remove empty genes
mat <- mat[Matrix::rowSums(mat) > 0, ]

# run EmptyDrops
ed <- DropletUtils::emptyDrops(mat, lower = 500)

# filter cells
mat_filtered <- mat[, !is.na(ed$FDR) & ed$FDR < 0.001]

hist(Matrix::colSums(mat_filtered), breaks = 200)
summary(Matrix::colSums(mat_filtered))
res_mat_filtered <- lapply(names(res_mat_Sep_list), function(s) {
  
  
  mat <- res_mat_Sep_list[[s]]
  
  # optional but recommended: remove all-zero genes early
  mat <- mat[Matrix::rowSums(mat) > 0, ]
  
  # EmptyDrops filtering
  run_emptydrops(mat, fdr_cutoff = 0.001)
})

## check number of cells left after UMI filtering (empty drops)
lapply(res_mat_filtered, dim)


hist(Matrix::colSums(res_mat_filtered[[1]]), breaks = 200,log = "x")
hist(Matrix::colSums(res_mat_filtered[[2]]), breaks = 200,log = "x")
hist(Matrix::colSums(res_mat_filtered[[3]]), breaks = 200,log = "x")
hist(Matrix::colSums(res_mat_filtered[[4]]), breaks = 200,log = "x")
hist(Matrix::colSums(res_mat_filtered[[5]]), breaks = 200,log = "x")
hist(Matrix::colSums(res_mat_filtered[[6]]), breaks = 200,log = "x")


##summary
summary(Matrix::colSums(res_mat_filtered[[1]]))
summary(Matrix::colSums(res_mat_filtered[[2]]))
summary(Matrix::colSums(res_mat_filtered[[3]]))
summary(Matrix::colSums(res_mat_filtered[[4]]))
summary(Matrix::colSums(res_mat_filtered[[5]]))
summary(Matrix::colSums(res_mat_filtered[[6]]))


##create seuratobject
seuObj_list<- lapply(res_mat_filtered, CreateSeuratObject)
saveRDS(seuObj_list, file = "/media/rna/VAQUITA/snRNAseq_potato/03_ANALYSIS/seuObj_list.rds")
seuObj_list <- readRDS("/media/rna/VAQUITA/snRNAseq_potato/03_ANALYSIS/seuObj_list.rds")

###Now find doublets...
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
library(DoubletFinder)

PreprocessingStep<- function(seuratObj, number_PCs_PCA=30, number_PCs_UMAP=10){
  
  ### classic seurat normalization, 30 pcs
  print("Running PCA, Normalization")
  seuratObj <- NormalizeData(seuratObj)
  seuratObj <- FindVariableFeatures(seuratObj)
  seuratObj <- ScaleData(seuratObj)
  seuratObj <- RunPCA(seuratObj, npcs = number_PCs_PCA)
  
  ###find local neighbourhood clusters with just 10 pcs
  
  print("Running UMAP, finding clusters")
  seuratObj <- RunUMAP(seuratObj, dims = 1:number_PCs_UMAP)
  seuratObj <- FindNeighbors(object = seuratObj, dims = 1:number_PCs_UMAP)              
  seuratObj <- FindClusters(object = seuratObj, resolution = 0.1)
  
  return(seuratObj)
  
}

####preprocesss, running normalization, PCA, UMAP before doublet finding
seuObj_list<- lapply(seuObj_list, PreprocessingStep)
ElbowPlot(seuObj_list[[1]], ndims = 30)
ElbowPlot(seuObj_list[[2]], ndims = 30)
ElbowPlot(seuObj_list[[3]], ndims = 30)
ElbowPlot(seuObj_list[[4]], ndims = 30)
ElbowPlot(seuObj_list[[5]], ndims = 30)
ElbowPlot(seuObj_list[[6]], ndims = 30)

#
#test_seu_CK1 
#test_seu_CK1 <- NormalizeData(seuObj_list[[1]])
#test_seu_CK1 <- FindVariableFeatures(test_seu_CK1)
#test_seu_CK1 <- ScaleData(test_seu_CK1)
#test_seu_CK1 <- RunPCA(test_seu_CK1, npcs = 30)


#stdv <- test_seu_CK1[["pca"]]@stdev
#
#percent_var <- (stdv^2/sum(stdv^2)) * 100
#cumulative_var <- cumsum(percent_var)
#co1 <- which(cumulative_var > 90)[1]
#co2 <- which(diff(percent_var) < 1)[1] + 1
#min_pc <- min(co1, co2)

# Finish pre-processing with 
#test_seu_CK1 <- RunUMAP(test_seu_CK1, dims = 1:10)
#test_seu_CK1 <- FindNeighbors(object = test_seu_CK1, dims = 1:10)              
#test_seu_CK1 <- FindClusters(object = test_seu_CK1, resolution = 0.1)


###identify pK, no ground truth
###Ground truth = experimentally known labels, eg you already know which cells are doublets and singlets
##can happen like genetic demultiplexing (e.g. SNP-based),cell hashing (antibody barcodes)
pK_identification<- function(seuObj, number_PCs=10){
  
  ##no ground truth
  seuObj_sweep.res <- paramSweep(seuObj, PCs = 1:number_PCs, sct = FALSE)
  seuObj_sweep.stats <- summarizeSweep(seuObj_sweep.res, GT = FALSE)
  seuObj_CK1_bcmvn <- find.pK(seuObj_sweep.stats)
  seuObj_CK1_bcmvn$pK <- as.numeric(as.character(seuObj_CK1_bcmvn$pK))
  
  return(seuObj_CK1_bcmvn)
}


CK1_bcmvn<- pK_identification(seuObj_list[[1]])
CK2_bcmvn<- pK_identification(seuObj_list[[2]])
CK3_bcmvn<- pK_identification(seuObj_list[[3]])
HS1_bcmvn<- pK_identification(seuObj_list[[4]])
HS2_bcmvn<- pK_identification(seuObj_list[[5]])
HS3_bcmvn<- pK_identification(seuObj_list[[6]])

bcmvn_list<- list(
  CK1 = CK1_bcmvn,
  CK2 = CK2_bcmvn,
  CK3 = CK3_bcmvn,
  HS1 = HS1_bcmvn,
  HS2 = HS2_bcmvn,
  HS3 = HS3_bcmvn
)


## pK Identification (no ground-truth) ------
#test_seu_CK1_sweep.res <- paramSweep(test_seu_CK1, PCs = 1:10, sct = FALSE,)
#test_seu_CK1_sweep.stats <- summarizeSweep(test_seu_CK1_sweep.res, GT = FALSE)
#test_seu_CK1_bcmvn <- find.pK(test_seu_CK1_sweep.stats)
#test_seu_CK1_bcmvn$pK <- as.numeric(as.character(test_seu_CK1_bcmvn$pK))
#
#test_seu_CK1_bcmvn$pK[which.max(test_seu_CK1_bcmvn$BCmetric)]

##plot the pKs

gridExtra::grid.arrange(
  ggplot(bcmvn_list[["CK1"]], aes(x = pK, y = BCmetric)) + geom_point() + geom_line() + theme_classic(),
  ggplot(bcmvn_list[["CK2"]], aes(x = pK, y = BCmetric)) + geom_point() + geom_line() + theme_classic(),
  ggplot(bcmvn_list[["CK3"]], aes(x = pK, y = BCmetric)) + geom_point() + geom_line() + theme_classic(),
  ggplot(bcmvn_list[["HS1"]], aes(x = pK, y = BCmetric)) + geom_point() + geom_line() + theme_classic(),
  ggplot(bcmvn_list[["HS2"]], aes(x = pK, y = BCmetric)) + geom_point() + geom_line() + theme_classic(),
  ggplot(bcmvn_list[["HS3"]], aes(x = pK, y = BCmetric)) + geom_point() + geom_line() + theme_classic()
)

bcmvn_list$CK1
bcmvn_list$CK1[which.max(bcmvn_list$CK1$BCmetric), ]
bcmvn_list$CK2[which.max(bcmvn_list$CK2$BCmetric), ]
bcmvn_list$CK3[which.max(bcmvn_list$CK3$BCmetric), ]
bcmvn_list$HS1[which.max(bcmvn_list$HS1$BCmetric), ]
bcmvn_list$HS2[which.max(bcmvn_list$HS2$BCmetric), ]
bcmvn_list$HS3[which.max(bcmvn_list$HS3$BCmetric), ]


calculate_nExp<- function(seuObj){
  
  ###### Homotypic Doublet Proportion Estimate -------
  homotypic.prop <- modelHomotypic(seuObj$seurat_clusters)
  nExp <- round(0.075 * ncol(seuObj))
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  
  return(nExp.adj)
}

nExp_list<- lapply(seuObj_list, calculate_nExp)

seuObj_list[[1]] <- doubletFinder(seuObj_list[[1]], PCs = 1:10, pN = 0.25, pK = 0.03, nExp = 700, reuse.pANN = NULL, sct = FALSE)
seuObj_list[[2]] <- doubletFinder(seuObj_list[[2]], PCs = 1:10, pN = 0.25, pK = 0.1, nExp = 1166, reuse.pANN = NULL, sct = FALSE)
seuObj_list[[3]] <- doubletFinder(seuObj_list[[3]], PCs = 1:10, pN = 0.25, pK = 0.055, nExp = 852, reuse.pANN = NULL, sct = FALSE)
seuObj_list[[4]] <- doubletFinder(seuObj_list[[4]], PCs = 1:10, pN = 0.25, pK = 0.23, nExp = 386, reuse.pANN = NULL, sct = FALSE)
seuObj_list[[5]] <- doubletFinder(seuObj_list[[5]], PCs = 1:10, pN = 0.25, pK = 0.02, nExp = 889, reuse.pANN = NULL, sct = FALSE)
seuObj_list[[6]] <- doubletFinder(seuObj_list[[6]], PCs = 1:10, pN = 0.25, pK = 0.11, nExp = 799, reuse.pANN = NULL, sct = FALSE)




table(seuObj_list[[1]]$DF.classifications_0.25_0.03_700)


DimPlot(seuObj_list[[1]], group.by = "DF.classifications_0.25_0.03_700")
DimPlot(seuObj_list[[2]], group.by = "DF.classifications_0.25_0.1_1166")
DimPlot(seuObj_list[[3]], group.by = "DF.classifications_0.25_0.055_852")
DimPlot(seuObj_list[[4]], group.by = "DF.classifications_0.25_0.23_386")
DimPlot(seuObj_list[[5]], group.by = "DF.classifications_0.25_0.02_889")
DimPlot(seuObj_list[[6]], group.by = "DF.classifications_0.25_0.11_799") )

gridExtra::grid.arrange(DimPlot(seuObj_list[[2]], group.by = "DF.classifications_0.25_0.1_1166"),
                        DimPlot(seuObj_list[[2]], group.by = "seurat_clusters"))


table(seuObj_list[[2]]$seurat_clusters, seuObj_list[[2]]$DF.classifications_0.25_0.1_1166)

aggregate(nCount_RNA ~ DF.classifications_0.25_0.03_700, data = seuObj_list[[1]]@meta.data, FUN = summary)
aggregate(nCount_RNA ~ DF.classifications_0.25_0.1_1166, data = seuObj_list[[2]]@meta.data, FUN = summary)
aggregate(nCount_RNA ~ DF.classifications_0.25_0.055_852, data = seuObj_list[[3]]@meta.data, FUN = summary)
aggregate(nCount_RNA ~ DF.classifications_0.25_0.23_386, data = seuObj_list[[4]]@meta.data, FUN = summary)
aggregate(nCount_RNA ~ DF.classifications_0.25_0.02_889, data = seuObj_list[[5]]@meta.data, FUN = summary)
aggregate(nCount_RNA ~ DF.classifications_0.25_0.11_799, data = seuObj_list[[6]]@meta.data, FUN = summary)

###get mitochondrial genes, determine cells that have undergo apoptosis
Mitochrondrial_Genes<- S.tuberosum.v6.1_latest %>% filter(str_detect(MAPMan.functional.group.basieredn.auf.v3.1,"mitochondrial electron transport "))
Mitochrondrial_Genes$locusName <- paste0(Mitochrondrial_Genes$locusName, ".v6.1")


DetectMTgenesAnd_addPercentFeature<- function(seuObj){
  
  ##Mitochrondrial_Genes$locusName is the total MT genes annotated to potato reference based on mapmap
  
  
  ##determine which are annotated MT genes in the seurat object
  mt_genes <- intersect(Mitochrondrial_Genes$locusName,
                        rownames(seuObj))
  
  #get the expression of MT genes
  seuObj$percent_mt <- PercentageFeatureSet(seuObj, features = mt_genes
  )
  
  return(seuObj)
}

seuObj_list<- lapply(seuObj_list, DetectMTgenesAnd_addPercentFeature)



VlnPlot(seuObj_list[[1]],group.by = "DF.classifications_0.25_0.03_700",features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),ncol = 3,pt.size = 0)
VlnPlot(seuObj_list[[2]],group.by = "DF.classifications_0.25_0.1_1166",features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),ncol = 3,pt.size = 0)
VlnPlot(seuObj_list[[3]],group.by = "DF.classifications_0.25_0.055_852",features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),ncol = 3,pt.size = 0)
VlnPlot(seuObj_list[[4]],group.by = "DF.classifications_0.25_0.23_386",features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),ncol = 3,pt.size = 0)
VlnPlot(seuObj_list[[5]],group.by = "DF.classifications_0.25_0.02_889",features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),ncol = 3,pt.size = 0)
VlnPlot(seuObj_list[[6]],group.by = "DF.classifications_0.25_0.11_799",features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),ncol = 3,pt.size = 0)


###sumarize total doublet found compared to singlet  
prop.table(table(seuObj_list[[1]][["DF.classifications_0.25_0.03_700"]])) * 100
prop.table(table(seuObj_list[[2]][["DF.classifications_0.25_0.1_1166"]])) * 100
prop.table(table(seuObj_list[[3]][["DF.classifications_0.25_0.055_852"]])) * 100
prop.table(table(seuObj_list[[4]][["DF.classifications_0.25_0.23_386"]])) * 100
prop.table(table(seuObj_list[[5]][["DF.classifications_0.25_0.02_889"]])) * 100
prop.table(table(seuObj_list[[6]][["DF.classifications_0.25_0.11_799"]])) * 100

##filter seuratobject
seuObj_filtered_list <- lapply(seuObj_list, function(seu) {
  
  ### because every df.classication has a different name for each sample
  Doublet_column<-  grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
  
  subset(seu,
         subset =
           seu[[Doublet_column]] == "Singlet" &
           nFeature_RNA > 500 &
           percent_mt < 25)
})

###check initial table
lapply(seuObj_list, function(seu) {
  ### because every df.classication has a different name for each sample
  Doublet_column<-  grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
  
  table(seu[[Doublet_column]])
})

##check final filter table...
lapply(seuObj_filtered_list, function(seu) {
  ### because every df.classication has a different name for each sample
  Doublet_column<-  grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
  
  table(seu[[Doublet_column]])
})

### add unique identitites to each seu object
seuObj_filtered_list[[1]]$orig.ident <- paste0("CK")
seuObj_filtered_list[[2]]$orig.ident <- paste0("CK")
seuObj_filtered_list[[3]]$orig.ident <- paste0("CK")
seuObj_filtered_list[[4]]$orig.ident <- paste0("HS")
seuObj_filtered_list[[5]]$orig.ident <- paste0("HS")
seuObj_filtered_list[[6]]$orig.ident <- paste0("HS")



###merge all seurat objects together.... , all cells from all samples
merged_seu <- merge(
  x = seuObj_filtered_list[[1]],
  y = seuObj_filtered_list[2:6],
  # add.cell.ids = paste0("Sample", 1:6)
)

saveRDS(merged_seu, file = "/media/rna/VAQUITA/snRNAseq_potato/03_ANALYSIS/merged_seu.rds")