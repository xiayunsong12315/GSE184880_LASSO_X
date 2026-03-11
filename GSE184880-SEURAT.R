library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(DoubletFinder)
library(SingleR)
library(GSVA)
library(GSEA)
library(GSEABase)
library(clusterProfiler)
library(limma)
library(pheatmap)



setwd('/data/xys/20250428_c4_pd1/')


GSM5599220_Norm1 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599220_Norm1",
                               gene.column = 1)
GSM5599221_Norm2 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599221_Norm2",
                               gene.column = 1)
GSM5599222_Norm3 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599222_Norm3",
                               gene.column = 1)
GSM5599223_Norm4 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599223_Norm4",
                               gene.column = 1)
GSM5599224_Norm5 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599224_Norm5",
                               gene.column = 1)
GSM5599225_cancer1 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599225_Cancer1",
                               gene.column = 1)
GSM5599226_cancer2 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599226_Cancer2",
                               gene.column = 1)
GSM5599227_cancer3 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599227_Cancer3",
                               gene.column = 1)
GSM5599228_cancer4 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599228_Cancer4",
                               gene.column = 1)
GSM5599229_cancer5 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599229_Cancer5",
                               gene.column = 1)
GSM5599230_cancer6 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599230_Cancer6",
                                gene.column = 1)
GSM5599231_cancer7 <- Read10X(data.dir = "/home/data/xys/work2/GSM184880_RAW/GSM5599231_Cancer7",
                              gene.column = 1)


GSM5599220_Norm1 <- CreateSeuratObject(counts = GSM5599220_Norm1, project = "GSM5599220_Norm1", min.cells = 3, min.features = 200)
GSM5599221_Norm2 <- CreateSeuratObject(counts = GSM5599221_Norm2, project = "GSM5599221_Norm2", min.cells = 3, min.features = 200)
GSM5599222_Norm3 <- CreateSeuratObject(counts = GSM5599222_Norm3, project = "GSM5599222_Norm3", min.cells = 3, min.features = 200)
GSM5599223_Norm4 <- CreateSeuratObject(counts = GSM5599223_Norm4, project = "GSM5599223_Norm4", min.cells = 3, min.features = 200)
GSM5599224_Norm5 <- CreateSeuratObject(counts = GSM5599224_Norm5, project = "GSM5599224_Norm5", min.cells = 3, min.features = 200)
GSM5599225_cancer1 <- CreateSeuratObject(counts = GSM5599225_cancer1, project = "GSM5599225_cancer1", min.cells = 3, min.features = 200)
GSM5599226_cancer2 <- CreateSeuratObject(counts = GSM5599226_cancer2, project = "GSM5599226_cancer2", min.cells = 3, min.features = 200)
GSM5599227_cancer3 <- CreateSeuratObject(counts = GSM5599227_cancer3, project = "GSM5599227_cancer3", min.cells = 3, min.features = 200)
GSM5599228_cancer4 <- CreateSeuratObject(counts = GSM5599228_cancer4, project = "GSM5599228_cancer4", min.cells = 3, min.features = 200)
GSM5599229_cancer5 <- CreateSeuratObject(counts = GSM5599229_cancer5, project = "GSM5599229_cancer5", min.cells = 3, min.features = 200)
GSM5599230_cancer6 <- CreateSeuratObject(counts = GSM5599230_cancer6, project = "GSM5599230_cancer6", min.cells = 3, min.features = 200)
GSM5599231_cancer7 <- CreateSeuratObject(counts = GSM5599231_cancer7, project = "GSM5599231_cancer7", min.cells = 3, min.features = 200)


GSM5599220_Norm1[["percent.mt"]] <- PercentageFeatureSet(GSM5599220_Norm1, pattern = "^mt-")
GSM5599221_Norm2[["percent.mt"]] <- PercentageFeatureSet(GSM5599221_Norm2, pattern = "^mt-")
GSM5599222_Norm3[["percent.mt"]] <- PercentageFeatureSet(GSM5599222_Norm3, pattern = "^mt-")
GSM5599223_Norm4[["percent.mt"]] <- PercentageFeatureSet(GSM5599223_Norm4, pattern = "^mt-")
GSM5599224_Norm5[["percent.mt"]] <- PercentageFeatureSet(GSM5599224_Norm5, pattern = "^mt-")
GSM5599225_cancer1[["percent.mt"]] <- PercentageFeatureSet(GSM5599225_cancer1, pattern = "^mt-")
GSM5599226_cancer2[["percent.mt"]] <- PercentageFeatureSet(GSM5599226_cancer2, pattern = "^mt-")
GSM5599227_cancer3[["percent.mt"]] <- PercentageFeatureSet(GSM5599227_cancer3, pattern = "^mt-")
GSM5599228_cancer4[["percent.mt"]] <- PercentageFeatureSet(GSM5599228_cancer4, pattern = "^mt-")
GSM5599229_cancer5[["percent.mt"]] <- PercentageFeatureSet(GSM5599229_cancer5, pattern = "^mt-")
GSM5599230_cancer6[["percent.mt"]] <- PercentageFeatureSet(GSM5599230_cancer6, pattern = "^mt-")
GSM5599231_cancer7[["percent.mt"]] <- PercentageFeatureSet(GSM5599231_cancer7, pattern = "^mt-")


HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599220_Norm1@assays$RNA)) 
HB.genes <- rownames(GSM5599220_Norm1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599220_Norm1[["percent.HB"]]<-PercentageFeatureSet(GSM5599220_Norm1, features=HB.genes) 


HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599221_Norm2@assays$RNA)) 
HB.genes <- rownames(GSM5599221_Norm2@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599221_Norm2[["percent.HB"]]<-PercentageFeatureSet(GSM5599221_Norm2, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599222_Norm3@assays$RNA)) 
HB.genes <- rownames(GSM5599222_Norm3@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599222_Norm3[["percent.HB"]]<-PercentageFeatureSet(GSM5599222_Norm3, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599223_Norm4@assays$RNA)) 
HB.genes <- rownames(GSM5599223_Norm4@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599223_Norm4[["percent.HB"]]<-PercentageFeatureSet(GSM5599223_Norm4, features=HB.genes) 

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599224_Norm5@assays$RNA)) 
HB.genes <- rownames(GSM5599224_Norm5@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599224_Norm5[["percent.HB"]]<-PercentageFeatureSet(GSM5599224_Norm5, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599225_cancer1@assays$RNA)) 
HB.genes <- rownames(GSM5599225_cancer1@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599225_cancer1[["percent.HB"]]<-PercentageFeatureSet(GSM5599225_cancer1, features=HB.genes) 

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599226_cancer2@assays$RNA)) 
HB.genes <- rownames(GSM5599226_cancer2@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599226_cancer2[["percent.HB"]]<-PercentageFeatureSet(GSM5599226_cancer2, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599227_cancer3@assays$RNA)) 
HB.genes <- rownames(GSM5599227_cancer3@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599227_cancer3[["percent.HB"]]<-PercentageFeatureSet(GSM5599227_cancer3, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599228_cancer4@assays$RNA)) 
HB.genes <- rownames(GSM5599228_cancer4@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599228_cancer4[["percent.HB"]]<-PercentageFeatureSet(GSM5599228_cancer4, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599229_cancer5@assays$RNA)) 
HB.genes <- rownames(GSM5599229_cancer5@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599229_cancer5[["percent.HB"]]<-PercentageFeatureSet(GSM5599229_cancer5, features=HB.genes)

HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599230_cancer6@assays$RNA)) 
HB.genes <- rownames(GSM5599230_cancer6@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599230_cancer6[["percent.HB"]]<-PercentageFeatureSet(GSM5599230_cancer6, features=HB.genes) 


HB.genes <- c("Hba-a1","Hba-a2","Hba-x","Hbb","Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y")
HB_m <- match(HB.genes, rownames(GSM5599231_cancer7@assays$RNA)) 
HB.genes <- rownames(GSM5599231_cancer7@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
GSM5599231_cancer7[["percent.HB"]]<-PercentageFeatureSet(GSM5599231_cancer7, features=HB.genes) 




GSM5599220_Norm1 <- subset(GSM5599220_Norm1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599221_Norm2 <- subset(GSM5599221_Norm2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599222_Norm3 <- subset(GSM5599222_Norm3, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599223_Norm4 <- subset(GSM5599223_Norm4, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599224_Norm5 <- subset(GSM5599224_Norm5, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599225_cancer1 <- subset(GSM5599225_cancer1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599226_cancer2 <- subset(GSM5599226_cancer2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599227_cancer3 <- subset(GSM5599227_cancer3, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599228_cancer4 <- subset(GSM5599228_cancer4, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599229_cancer5 <- subset(GSM5599229_cancer5, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599230_cancer6 <- subset(GSM5599230_cancer6, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)
GSM5599231_cancer7 <- subset(GSM5599231_cancer7, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10 & percent.HB < 10)



GSM5599220_Norm1 <- NormalizeData(GSM5599220_Norm1)
GSM5599220_Norm1 <- FindVariableFeatures(GSM5599220_Norm1, selection.method = "vst", nfeatures = 2000)
GSM5599220_Norm1 <- ScaleData(GSM5599220_Norm1, features = VariableFeatures(object = GSM5599220_Norm1))
GSM5599220_Norm1 <- RunPCA(GSM5599220_Norm1, verbose = FALSE)
GSM5599220_Norm1 <- RunUMAP(GSM5599220_Norm1, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599220_Norm1, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599220_Norm1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599220_Norm1@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599220_Norm1 <- doubletFinder_v3(GSM5599220_Norm1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599220_Norm1 <- doubletFinder_v3(GSM5599220_Norm1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599220_Norm1@meta.data)[6], sct = FALSE)
colnames(GSM5599220_Norm1@meta.data)
GSM5599220_Norm1 <- subset(GSM5599220_Norm1, subset = DF.classifications_0.25_0.3_476 == "Singlet")

GSM5599221_Norm2 <- NormalizeData(GSM5599221_Norm2)
GSM5599221_Norm2 <- FindVariableFeatures(GSM5599221_Norm2, selection.method = "vst", nfeatures = 2000)
GSM5599221_Norm2 <- ScaleData(GSM5599221_Norm2, features = VariableFeatures(object = GSM5599221_Norm2))
GSM5599221_Norm2 <- RunPCA(GSM5599221_Norm2, verbose = FALSE)
GSM5599221_Norm2 <- RunUMAP(GSM5599221_Norm2, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599221_Norm2, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599221_Norm2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599221_Norm2@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599221_Norm2 <- doubletFinder_v3(GSM5599221_Norm2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599221_Norm2 <- doubletFinder_v3(GSM5599221_Norm2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599221_Norm2@meta.data)[6], sct = FALSE)
colnames(GSM5599221_Norm2@meta.data)
GSM5599221_Norm2 <- subset(GSM5599221_Norm2, subset = DF.classifications_0.25_0.18_386 == "Singlet")



GSM5599222_Norm3 <- NormalizeData(GSM5599222_Norm3)
GSM5599222_Norm3 <- FindVariableFeatures(GSM5599222_Norm3, selection.method = "vst", nfeatures = 2000)
GSM5599222_Norm3 <- ScaleData(GSM5599222_Norm3, features = VariableFeatures(object = GSM5599222_Norm3))
GSM5599222_Norm3 <- RunPCA(GSM5599222_Norm3, verbose = FALSE)
GSM5599222_Norm3 <- RunUMAP(GSM5599222_Norm3, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599222_Norm3, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599222_Norm3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599222_Norm3@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599222_Norm3 <- doubletFinder_v3(GSM5599222_Norm3, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599222_Norm3 <- doubletFinder_v3(GSM5599222_Norm3, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599222_Norm3@meta.data)[6], sct = FALSE)
colnames(GSM5599222_Norm3@meta.data)
GSM5599222_Norm3 <- subset(GSM5599222_Norm3, subset = DF.classifications_0.25_0.29_333 == "Singlet")

GSM5599223_Norm4 <- NormalizeData(GSM5599223_Norm4)
GSM5599223_Norm4 <- FindVariableFeatures(GSM5599223_Norm4, selection.method = "vst", nfeatures = 2000)
GSM5599223_Norm4 <- ScaleData(GSM5599223_Norm4, features = VariableFeatures(object = GSM5599223_Norm4))
GSM5599223_Norm4 <- RunPCA(GSM5599223_Norm4, verbose = FALSE)
GSM5599223_Norm4 <- RunUMAP(GSM5599223_Norm4, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599223_Norm4, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599223_Norm4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599223_Norm4@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599223_Norm4 <- doubletFinder_v3(GSM5599223_Norm4, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599223_Norm4 <- doubletFinder_v3(GSM5599223_Norm4, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599223_Norm4@meta.data)[6], sct = FALSE)
colnames(GSM5599223_Norm4@meta.data)
GSM5599223_Norm4 <- subset(GSM5599223_Norm4, subset = DF.classifications_0.25_0.21_497 == "Singlet")



GSM5599224_Norm5 <- NormalizeData(GSM5599224_Norm5)
GSM5599224_Norm5 <- FindVariableFeatures(GSM5599224_Norm5, selection.method = "vst", nfeatures = 2000)
GSM5599224_Norm5 <- ScaleData(GSM5599224_Norm5, features = VariableFeatures(object = GSM5599224_Norm5))
GSM5599224_Norm5 <- RunPCA(GSM5599224_Norm5, verbose = FALSE)
GSM5599224_Norm5 <- RunUMAP(GSM5599224_Norm5, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599224_Norm5, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599224_Norm5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599224_Norm5@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599224_Norm5 <- doubletFinder_v3(GSM5599224_Norm5, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599224_Norm5 <- doubletFinder_v3(GSM5599224_Norm5, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599224_Norm5@meta.data)[6], sct = FALSE)
colnames(GSM5599224_Norm5@meta.data)
GSM5599224_Norm5 <- subset(GSM5599224_Norm5, subset = DF.classifications_0.25_0.005_321 == "Singlet")




GSM5599225_cancer1 <- NormalizeData(GSM5599225_cancer1)
GSM5599225_cancer1 <- FindVariableFeatures(GSM5599225_cancer1, selection.method = "vst", nfeatures = 2000)
GSM5599225_cancer1 <- ScaleData(GSM5599225_cancer1, features = VariableFeatures(object = GSM5599225_cancer1))
GSM5599225_cancer1 <- RunPCA(GSM5599225_cancer1, verbose = FALSE)
GSM5599225_cancer1 <- RunUMAP(GSM5599225_cancer1, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599225_cancer1, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599225_cancer1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599225_cancer1@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599225_cancer1 <- doubletFinder_v3(GSM5599225_cancer1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599225_cancer1 <- doubletFinder_v3(GSM5599225_cancer1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599225_cancer1@meta.data)[6], sct = FALSE)
colnames(GSM5599225_cancer1@meta.data)
GSM5599225_cancer1 <- subset(GSM5599225_cancer1, subset = DF.classifications_0.25_0.21_651 == "Singlet")


GSM5599226_cancer2 <- NormalizeData(GSM5599226_cancer2)
GSM5599226_cancer2 <- FindVariableFeatures(GSM5599226_cancer2, selection.method = "vst", nfeatures = 2000)
GSM5599226_cancer2 <- ScaleData(GSM5599226_cancer2, features = VariableFeatures(object = GSM5599226_cancer2))
GSM5599226_cancer2 <- RunPCA(GSM5599226_cancer2, verbose = FALSE)
GSM5599226_cancer2 <- RunUMAP(GSM5599226_cancer2, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599226_cancer2, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599226_cancer2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599226_cancer2@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599226_cancer2 <- doubletFinder_v3(GSM5599226_cancer2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599226_cancer2 <- doubletFinder_v3(GSM5599226_cancer2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599226_cancer2@meta.data)[6], sct = FALSE)
colnames(GSM5599226_cancer2@meta.data)
GSM5599226_cancer2 <- subset(GSM5599226_cancer2, subset = DF.classifications_0.25_0.005_315 == "Singlet")

GSM5599227_cancer3 <- NormalizeData(GSM5599227_cancer3)
GSM5599227_cancer3 <- FindVariableFeatures(GSM5599227_cancer3, selection.method = "vst", nfeatures = 2000)
GSM5599227_cancer3 <- ScaleData(GSM5599227_cancer3, features = VariableFeatures(object = GSM5599227_cancer3))
GSM5599227_cancer3 <- RunPCA(GSM5599227_cancer3, verbose = FALSE)
GSM5599227_cancer3 <- RunUMAP(GSM5599227_cancer3, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599227_cancer3, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599227_cancer3@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599227_cancer3@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599227_cancer3 <- doubletFinder_v3(GSM5599227_cancer3, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599227_cancer3 <- doubletFinder_v3(GSM5599227_cancer3, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599227_cancer3@meta.data)[6], sct = FALSE)
colnames(GSM5599227_cancer3@meta.data)
GSM5599227_cancer3 <- subset(GSM5599227_cancer3, subset = DF.classifications_0.25_0.005_372 == "Singlet")

GSM5599228_cancer4 <- NormalizeData(GSM5599228_cancer4)
GSM5599228_cancer4 <- FindVariableFeatures(GSM5599228_cancer4, selection.method = "vst", nfeatures = 2000)
GSM5599228_cancer4 <- ScaleData(GSM5599228_cancer4, features = VariableFeatures(object = GSM5599228_cancer4))
GSM5599228_cancer4 <- RunPCA(GSM5599228_cancer4, verbose = FALSE)
GSM5599228_cancer4 <- RunUMAP(GSM5599228_cancer4, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599228_cancer4, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599228_cancer4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599228_cancer4@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599228_cancer4 <- doubletFinder_v3(GSM5599228_cancer4, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599228_cancer4 <- doubletFinder_v3(GSM5599228_cancer4, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599228_cancer4@meta.data)[6], sct = FALSE)
colnames(GSM5599228_cancer4@meta.data)
GSM5599228_cancer4 <- subset(GSM5599228_cancer4, subset = DF.classifications_0.25_0.01_295 == "Singlet")



GSM5599229_cancer5 <- NormalizeData(GSM5599229_cancer5)
GSM5599229_cancer5 <- FindVariableFeatures(GSM5599229_cancer5, selection.method = "vst", nfeatures = 2000)
GSM5599229_cancer5 <- ScaleData(GSM5599229_cancer5, features = VariableFeatures(object = GSM5599229_cancer5))
GSM5599229_cancer5 <- RunPCA(GSM5599229_cancer5, verbose = FALSE)
GSM5599229_cancer5 <- RunUMAP(GSM5599229_cancer5, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599229_cancer5, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599229_cancer5@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599229_cancer5@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599229_cancer5 <- doubletFinder_v3(GSM5599229_cancer5, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599229_cancer5 <- doubletFinder_v3(GSM5599229_cancer5, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599229_cancer5@meta.data)[6], sct = FALSE)
colnames(GSM5599229_cancer5@meta.data)
GSM5599229_cancer5 <- subset(GSM5599229_cancer5, subset = DF.classifications_0.25_0.29_431 == "Singlet")


GSM5599230_cancer6 <- NormalizeData(GSM5599230_cancer6)
GSM5599230_cancer6 <- FindVariableFeatures(GSM5599230_cancer6, selection.method = "vst", nfeatures = 2000)
GSM5599230_cancer6 <- ScaleData(GSM5599230_cancer6, features = VariableFeatures(object = GSM5599230_cancer6))
GSM5599230_cancer6 <- RunPCA(GSM5599230_cancer6, verbose = FALSE)
GSM5599230_cancer6 <- RunUMAP(GSM5599230_cancer6, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599230_cancer6, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599230_cancer6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599230_cancer6@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599230_cancer6 <- doubletFinder_v3(GSM5599230_cancer6, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599230_cancer6 <- doubletFinder_v3(GSM5599230_cancer6, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599230_cancer6@meta.data)[6], sct = FALSE)
colnames(GSM5599230_cancer6@meta.data)
GSM5599230_cancer6 <- subset(GSM5599230_cancer6, subset = DF.classifications_0.25_0.07_380 == "Singlet")


GSM5599231_cancer7 <- NormalizeData(GSM5599231_cancer7)
GSM5599231_cancer7 <- FindVariableFeatures(GSM5599231_cancer7, selection.method = "vst", nfeatures = 2000)
GSM5599231_cancer7 <- ScaleData(GSM5599231_cancer7, features = VariableFeatures(object = GSM5599231_cancer7))
GSM5599231_cancer7 <- RunPCA(GSM5599231_cancer7, verbose = FALSE)
GSM5599231_cancer7 <- RunUMAP(GSM5599231_cancer7, dims = 1:30)
sweep.res.list <- paramSweep_v3(GSM5599231_cancer7, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
annotations <- GSM5599231_cancer7@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(GSM5599231_cancer7@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
GSM5599231_cancer7 <- doubletFinder_v3(GSM5599231_cancer7, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
GSM5599231_cancer7 <- doubletFinder_v3(GSM5599231_cancer7, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(GSM5599231_cancer7@meta.data)[6], sct = FALSE)
colnames(GSM5599231_cancer7@meta.data)
GSM5599231_cancer7 <- subset(GSM5599231_cancer7, subset = DF.classifications_0.25_0.3_379 == "Singlet")





GSM184880 <- merge(x = GSM5599222_Norm3, y = list (GSM5599220_Norm1,GSM5599221_Norm2,GSM5599223_Norm4,GSM5599224_Norm5,GSM5599225_cancer1,GSM5599226_cancer2,GSM5599227_cancer3,GSM5599228_cancer4,GSM5599229_cancer5,GSM5599230_cancer6,GSM5599231_cancer7))
GSM184880
head(GSM184880@meta.data)
library("harmony")
GSM184880 <- NormalizeData(GSM184880, normalization.method = "LogNormalize", scale.factor = 10000)
GSM184880 <- FindVariableFeatures(GSM184880, selection.method = "vst", nfeatures = 2000)
GSM184880<-ScaleData(GSM184880)%>%RunPCA(verbose=FALSE)
system.time({GSM184880 <- RunHarmony(GSM184880, group.by.vars = "orig.ident")})

GSM184880<-FindNeighbors(GSM184880, reduction = "harmony",dims = 1:30)
GSM184880<-FindClusters(GSM184880, reduction = "harmony",resolution = 0.2)
GSM184880<-RunUMAP(GSM184880,dims = 1:30,reduction = "harmony")
#GSM184880<-RunTSNE(GSM184880,dims = 1:30,reduction = "harmony")

DimPlot(GSM184880, reduction = "umap",label = T, pt.size = 0.5)
DimPlot(GSM184880,reduction = "umap",label = T,pt.size = 0.5,split.by = "orig.ident")



library(Seurat)
library(dplyr)
library(ggplot2)

# 假设整合后的 Seurat 对象是 GSM184880

# -------------------------
# 1️⃣ 标准化 & 找高变基因
# -------------------------
GSM184880 <- NormalizeData(GSM184880)
GSM184880 <- FindVariableFeatures(GSM184880, selection.method = "vst", nfeatures = 2000)
GSM184880 <- ScaleData(GSM184880)

# -------------------------
# 2️⃣ PCA 降维
# -------------------------
GSM184880 <- RunPCA(GSM184880, features = VariableFeatures(GSM184880))
ElbowPlot(GSM184880)  # 查看主成分数量

# -------------------------
# 3️⃣ 邻近图 & 聚类
# -------------------------
GSM184880 <- FindNeighbors(GSM184880, dims = 1:20)  # 取前20个PCs，可根据ElbowPlot调整
GSM184880 <- FindClusters(GSM184880, resolution = 0.5)  # 分辨率可调

# -------------------------
# 4️⃣ UMAP 可视化
# -------------------------
GSM184880 <- RunUMAP(GSM184880, dims = 1:20)
DimPlot(GSM184880, reduction = "umap", label = TRUE)

# -------------------------
# 5️⃣ 细胞类型注释（根据已知标记基因）
# -------------------------
# 常见标记基因示例：
# EPCAM / KRT8 / KRT18 - 上皮细胞
# PTPRC (CD45) - 免疫细胞
# CD3D / CD3E - T 细胞
# MS4A1 (CD20) - B 细胞
# CD14 / LYZ - 单核细胞/巨噬细胞
# COL1A1 / ACTA2 - 基质细胞 / 成纤维细胞
# PECAM1 (CD31) - 内皮细胞



head(rownames(GSM184880))

library(biomaRt)
library(Seurat)
library(dplyr)

# 1️⃣ 连接 human Ensembl
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# 2️⃣ 获取 Seurat 对象中的 ENSEMBL ID
ensembl_ids <- rownames(GSM184880)  # 替换成你的对象名字

library(org.Hs.eg.db)
library(AnnotationDbi)


# 卸载旧包
remove.packages("org.Hs.eg.db")
remove.packages("AnnotationDbi")

# URL
url_gtf <- "ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"

# 下载
download.file(url_gtf, destfile = "Homo_sapiens.GRCh38.109.gtf.gz")

library(data.table)

# 读取 GTF
gtf <- fread("zcat Homo_sapiens.GRCh38.109.gtf.gz", sep="\t", header=FALSE, data.table=FALSE)

# GTF 第9列是 info 字段，包含 gene_id 和 gene_name
gtf_genes <- gtf[gtf$V3=="gene", ]  # 只保留 gene 行

# 提取 gene_id 和 gene_name
library(stringr)
gtf_genes$gene_id <- str_match(gtf_genes$V9, "gene_id \"(ENSG[0-9]+)\";")[,2]
gtf_genes$gene_name <- str_match(gtf_genes$V9, "gene_name \"([^\"]+)\";")[,2]

# 创建映射表
gene_map <- gtf_genes[, c("gene_id", "gene_name")]
head(gene_map)


library(Seurat)

# 假设 GSM184880 是 raw counts matrix
GSE184880 <- CreateSeuratObject(counts = GSM184880, project = "Ovarian_Cancer")


# 假设 gene_symbols 是和行顺序对应的向量
rownames(GSM184880) <- gene_symbols

# 去掉没有对应 symbol 的行
GSM184880 <- GSM184880[!is.na(rownames(GSM184880)), ]

# 检查前几行
head(rownames(GSM184880))

# counts
rownames(GSM184880[["RNA"]]@counts) <- gene_symbols

# normalized data
rownames(GSM184880[["RNA"]]@data) <- gene_symbols

# 不修改 scale.data，Seurat 会自动处理

genes_in_scaled <- rownames(GSM184880[["RNA"]]@scale.data)
# 匹配 gene_symbols
matched_symbols <- gene_symbols[match(genes_in_scaled, rownames(GSM184880[["RNA"]]@counts))]

rownames(GSM184880[["RNA"]]@scale.data) <- matched_symbols

# FeaturePlot
FeaturePlot(GSM184880, features = c("EPCAM","PTPRC","CD3D","MS4A1","CD14","COL1A1","PECAM1"))



FeaturePlot(GSM184880, features = c("CD3D","CD3E","CD2","CD24","KRT19","EPCAM","COL1A1","DCN","C1R","VWF", "CDH5","LYZ", "CD68", "TYROBP","CD79A", "CD79B", "MZB1"))

# PCA / 聚类 / UMAP
GSM184880 <- FindVariableFeatures(GSM184880)
GSM184880 <- ScaleData(GSM184880)
GSM184880 <- RunPCA(GSM184880)
GSM184880 <- FindNeighbors(GSM184880, dims = 1:20)
GSM184880 <- FindClusters(GSM184880, resolution = 0.5)
GSM184880 <- RunUMAP(GSM184880, dims = 1:20)
DimPlot(GSM184880, reduction = "umap", label = TRUE)

# 找每个 cluster 标志基因
celltype_markers <- FindAllMarkers(GSM184880, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

table(GSM184880$seurat_clusters)  # 查看 cluster 数量和分布

cluster_markers <- FindAllMarkers(
  GSM184880,
  only.pos = TRUE,
  min.pct = 0,
  logfc.threshold = 0
)

# 查看前几行
head(cluster_markers)

# 保存为 CSV 文件
write.csv(cluster_markers, file = "GSM184880_cluster_markers.csv", row.names = FALSE)



counts <- GSM184880[["RNA"]]@counts
rownames(counts) <- gene_symbols
counts <- counts[!is.na(rownames(counts)), ]

library(Seurat)

GSM184880 <- CreateSeuratObject(counts = counts)

GSM184880 <- NormalizeData(GSM184880)

GSM184880 <- FindVariableFeatures(GSM184880)

GSM184880 <- ScaleData(GSM184880)

GSM184880 <- RunPCA(GSM184880)

GSM184880 <- FindNeighbors(GSM184880, dims = 1:20)

GSM184880 <- FindClusters(GSM184880, resolution = 0.5)

GSM184880 <- RunUMAP(GSM184880, dims = 1:20)

cluster_markers <- FindAllMarkers(
  GSM184880,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1
)

write.csv(cluster_markers, "cluster_markers.csv", row.names = FALSE)

































# 1️⃣ 提取 Seurat 对象行名（ENSEMBL ID）
ensembl_ids <- rownames(GSM184880)

# 2️⃣ 对应 gene_name
gene_symbols <- gene_map$gene_name[match(ensembl_ids, gene_map$gene_id)]

# 3️⃣ 替换 Seurat 行名
rownames(GSM184880) <- gene_symbols

# 4️⃣ 去掉没有对应 symbol 的行
GSM184880 <- GSM184880[!is.na(rownames(GSM184880)), ]

# 检查前几行
head(rownames(GSM184880))





FeaturePlot(GSM184880, features = c("EPCAM","PTPRC","CD3D","MS4A1","CD14","COL1A1","PECAM1"))

# 根据聚类和标记基因为每个 cluster 赋予细胞类型
new.cluster.ids <- c(
  "T_cells", "B_cells", "Macrophages", "Epithelial", 
  "Fibroblasts", "Endothelial", "NK_cells", "Other"  # 根据实际cluster数量调整
)
names(new.cluster.ids) <- levels(GSM184880)
GSM184880 <- RenameIdents(GSM184880, new.cluster.ids)

# 可视化标注后的细胞类型
DimPlot(GSM184880, reduction = "umap", label = TRUE, label.size = 5)

# -------------------------
# 6️⃣ 每个细胞类型的标志基因
# -------------------------
# FindAllMarkers 会找每个 cluster vs 其他所有 cluster 的差异表达基因
celltype_markers <- FindAllMarkers(
  GSM184880,
  only.pos = TRUE,  # 只取上调基因
  min.pct = 0.25,   # 至少25%细胞表达
  logfc.threshold = 0.25
)

# 查看 top 5 标志基因
celltype_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 可以保存到 CSV
write.csv(celltype_markers, "celltype_markers_GSM184880.csv", row.names = FALSE)





















FeaturePlot(GSM184880, c("Ptprc","Pdgfra","Vwf"))
FeaturePlot(GSM184880, c("Krt14","Epcam","Krt10","Krt5","Krt17","Krt18","CLDN4","Krt6","Krt7","Krt20"),label = T)
FeaturePlot(GSM184880, "CD200")
FeaturePlot(bca_d.integrated, "Epcam")
FeaturePlot(bca_d.integrated, "Ptprc")#删除这群细胞
FeaturePlot(bca_d.integrated, "Vwf")#删除
FeaturePlot(bca_d.integrated, "Col1a1")
FeaturePlot(bca_d.integrated, "Sox2")
FeaturePlot(GSM184880, "Rhox5")
table(bca_d.integrated$seurat_clusters,bca_d.integrated$orig.ident)




bca_pd1_x.integrated <- subset(GSM184880, idents = c(3))
DimPlot(bca_pd1_x.integrated,reduction = "umap",label = T,pt.size = 0.5,split.by = "orig.ident")
bca_pd1_x.integrated <- NormalizeData(bca_pd1_x.integrated)
bca_pd1_x.integrated <- FindVariableFeatures(bca_pd1_x.integrated, selection.method = "vst", nfeatures = 2000)
bca_pd1_x.integrated <- ScaleData(bca_pd1_x.integrated, features = VariableFeatures(object = bca_pd1_x.integrated))
bca_pd1_x.integrated <- RunPCA(bca_pd1_x.integrated, features = VariableFeatures(object = bca_pd1_x.integrated))
ElbowPlot(bca_pd1_x.integrated,ndims = 50)
bca_pd1_x.integrated <- FindNeighbors(bca_pd1_x.integrated, dims = 1:20,k.param = 20)
bca_pd1_x.integrated <- FindClusters(bca_pd1_x.integrated, resolution = 0.3)
bca_pd1_x.integrated <- RunUMAP(bca_pd1_x.integrated, dims = 1:20)

DimPlot(bca_pd1_x.integrated,reduction = "umap",label = TRUE,pt.size = 0.5)
DimPlot(bca_pd1_x.integrated,reduction = "umap",label = T,pt.size = 0.5,split.by = "orig.ident")

#统计一下数量
table(combined_integrated$orig.ident)


# 给 bca.integrated 加一个 metadata 标记，用于高亮
cells_dWT_day14_1 <- colnames(subset_seurat)[subset_seurat$orig.ident == "dWT_day14_1"]

# 添加一个新列 highlight_group，标记是否为目标细胞
bca.integrated$highlight_group <- ifelse(colnames(bca.integrated) %in% cells_dWT_day14_1,
                                         "dWT_day14_1", "Other")

# 绘图：按 orig.ident 分组，按 highlight_group 上色
DimPlot(bca.integrated, reduction = "umap",
        group.by = "highlight_group", 
        split.by = "orig.ident",
        cols = c("dWT_day14_1" = "red", "Other" = "gray80")) +
  ggtitle("Highlighting dWT_day14_1 in split UMAPs (by orig.ident)")


#CCA处理新批次和旧批次数据
old_samples <- SplitObject(subset_seurat_x, split.by = "orig.ident")

# 设置 orig.ident（如果还没设置）
GSM5599222_Norm3$orig.ident <- " GSM5599222_Norm3"
GSM5599221_Norm2$orig.ident <- "GSM5599221_Norm2"
GSM5599223_Norm4$orig.ident <- "GSM5599223_Norm4"
GSM5599224_Norm5$orig.ident <- "GSM5599224_Norm5"
GSM5599225_cancer1$orig.ident <- "GSM5599225_cancer1"
GSM5599226_cancer2$orig.ident <- "GSM5599226_cancer2"

# 把新样本打包成列表
new_samples <- list(GSM5599222_Norm3,GSM5599223_Norm4,GSM5599227_cancer3,GSM5599226_cancer2)
all_samples <- c(old_samples, bca_pd1_x.integrated)
all_samples <- lapply(all_samples, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})
anchors <- FindIntegrationAnchors(object.list = all_samples, dims = 1:30)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(seurat_integrated) <- "integrated"

seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 30)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:15)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.2)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:15)
DimPlot(seurat_integrated, split.by = "orig.ident", reduction = "umap")
DimPlot(seurat_integrated,reduction = "umap",label = TRUE,pt.size = 0.5)
#统计一下数量
table(seurat_integrated$orig.ident)



# 指定你的样本顺序

seurat_integrated_p <- subset(seurat_integrated, subset = orig.ident %in% c("GSM5599222_Norm3", "GSM5599223_Norm4","dPD_1_day14_1","dPD_1_day14_2","dPD1_DAY_7.1","dPD1_DAY_7.2","dWT_DAY_7.1","dWT_DAY_7.2", "dWT_day14_1","GSM5599226_cancer2","GSM5599227_cancer3"))

DimPlot(seurat_integrated_p, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(seurat_integrated_p, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "orig.ident")
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.2)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:15)
seurat_integrated_p <- subset(seurat_integrated, subset = orig.ident %in% c("GSM5599222_Norm3", "GSM5599223_Norm4","dPD_1_day14_1","dPD_1_day14_2","dPD1_DAY_7.1","dPD1_DAY_7.2","dWT_DAY_7.1","dWT_DAY_7.2", "dWT_day14_1","GSM5599226_cancer2","GSM5599227_cancer3"))

DimPlot(seurat_integrated_p, split.by = "orig.ident", reduction = "umap")
DimPlot(seurat_integrated_p,reduction = "umap",label = TRUE,pt.size = 0.5)
table(seurat_integrated$seurat_clusters,seurat_integrated$orig.ident)

seurat_integrated_p <- readRDS("combined_integrated.rds")




