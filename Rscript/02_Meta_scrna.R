#########Required packages###################
#############################################
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(DoubletFinder)
library(clustree)
library(glmGamPoi)
library(ggsci)
library(paletteer)
library(reshape2)
library(readxl)
library(pheatmap)
library(bbknnR)
############################################
############################################

meta_info<-readxl::read_xlsx('/home/ps/work/project/ICCOR/scRNA/Meta/Meta_info.xlsx')
###matrix
dataset1<-list.dirs('scRNA_Data',full.names = TRUE,recursive = F)
samplename1<-basename(dataset1)
seuList1 <- lapply(dataset1,function(x){ 
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
names(seuList1)<-samplename1

#ht
dataset2<-list.files('E-MTAB-15983',full.names = T)
seuList2 <- lapply(dataset2,function(x){
  seu = CreateSeuratObject(counts = Read10X_h5(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samplename2<-gsub('_sample_feature_bc_matrix.h5','',basename(dataset2))
names(seuList2)<-samplename2

seuList<-c(seuList1,seuList2)
samplename<-c(samplename1,samplename2)

for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(Sample[i],ncol(sce))
  sce[["orig.ident"]] = orig
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = '^RPL|^RPS|^MRPL|^MRPS')
  
  sce <- subset(sce, subset = nFeature_RNA >= 200 & nFeature_RNA <= 8000 &
                  nCount_RNA >= 500 &
                  percent.mt <= 20 & percent.ribo < 40)
  seuList[[i]]<-sce
}

#DecontX
library(decontX)
for(i in names(seuList)){
  sce <- seuList[[i]]
  counts <- GetAssayData(object = sce, layer='counts')
  decontX_results <- decontX(counts)
  sce$Contamination <- decontX_results$contamination
  ncell<-c(ncell,sum(decontX_results$contamination<0.2))
  seuList[[i]]<-sce
  print(paste0(i," Finished!"))
}

#doublet removal
for(i in names(seuList)){
  sce <- seuList[[i]]
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE,num.cores=48)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])    
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  DoubletRate = ncol(sce)*8*1e-6 
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, sct = FALSE)
  seuList[[i]]<-sce
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

meta_seu <- merge(seuList[[1]],
                  y = seuList[-1],
                  add.cell.ids = names(seuList))
c <- grep("pANN_0.25",colnames(meta_seu@meta.data))
meta_seu@meta.data <- meta_seu@meta.data[,-c]
c <- grep("DF.classifications_0.25",colnames(meta_seu@meta.data))
meta_seu@meta.data <- meta_seu@meta.data[,-c]
meta_seu@meta.data[["cell_name"]]<-colnames(meta_seu)
meta_seu <- subset(meta_seu, subset = DF.classifications=="Singlet" & Contamination<0.2)

meta_seu <- NormalizeData(meta_seu)
meta_seu <- FindVariableFeatures(meta_seu,nfeatures=3000)
meta_seu <- ScaleData(meta_seu,vars.to.regress='percent.mt')
meta_seu <- RunPCA(meta_seu, features = VariableFeatures(object = meta_seu))
meta_seu@reductions[["scVI"]]<-sce_for_scvi@reductions[["scVI"]]
meta_seu <- FindNeighbors(meta_seu, dims = 1:30,reduction="scVI")
meta_seu <- FindClusters(meta_seu,resolution = 0.3)
meta_seu <- RunUMAP(meta_seu,reduction="scVI", dims = 1:30)
meta_seu <- RunTSNE(meta_seu,reduction="scVI", dims = 1:30)

new.cluster.ids <- c('T cell','Epithelial cell','T cell','Myeloid cell','Myeloid cell','Epithelial cell','Fibroblasts','Endothelial cell',
                     'Epithelial cell','Epithelial cell','T cell','Fibroblasts','Smooth muscle cell','Plasma cell','Epithelial cell',
                     'Epithelial cell','Fibroblasts','B cell','Epithelial cell','Epithelial cell','Epithelial cell','T cell','Epithelial cell',
                     'Epithelial cell','Mast cell','Myeloid cell','Epithelial cell','Endothelial cell','Plasma cell','Endothelial cell',
                     'Myeloid cell','pDC','Myeloid cell','Epithelial cell','Epithelial cell','Epithelial cell','Epithelial cell')
names(new.cluster.ids) <- levels(meta_seu)
meta_seu <- RenameIdents(meta_seu, new.cluster.ids)
meta_seu@meta.data[["cell_type"]]<-meta_seu@active.ident

ncell_layer<-unlist(lapply(meta_seu@assays[["RNA"]]@layers, ncol))[159:316]
sketch_ncells <- pmin(
  ncell_layer,
  pmax(round(ncell_layer * 0.1), 100)
)
meta_seu <- SketchData(meta_seu,
                       ncells = sketch_ncells,method = "Uniform",
                       sketched.assay = "sketch",features = VariableFeatures(meta_seu))

#Epithelial
Epithelial<-subset(meta_seu,subset=cell_type=='Epithelial cell')
Epithelial <- NormalizeData(Epithelial)
Epithelial <- FindVariableFeatures(Epithelial,nfeatures=3000)
Epithelial <- ScaleData(Epithelial,vars.to.regress='percent.mt')
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunBBKNN(Epithelial, batch_key = "orig.ident",run_TSNE=TRUE,run_UMAP=TRUE)
Epithelial <- FindClusters(Epithelial,resolution = 0.5,graph.name='RNA_bbknn')
new.cluster.ids <- c('Squamous epithelial cell','Progenitor-like epithelial cell','Cycling epithelial cell','Progenitor-like epithelial cell','Squamous epithelial cell',
                     'Glandular epithelial cell','Squamous epithelial cell','Basal-like epithelial cell','Ciliated epithelial cell')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
DimPlot(Epithelial, reduction = "umap", label = F,raster=FALSE,cols = color,pt.size = 0.05,alpha = 0.3)
ncell_layer<-unlist(lapply(Epithelial@assays[["RNA"]]@layers, ncol))[157:312]
sketch_ncells <- pmin(
  ncell_layer,
  pmax(round(ncell_layer * 0.1), 100)
)
Epithelial <- SketchData(Epithelial,
                         ncells = sketch_ncells,
                         sketched.assay = "sketch",features = VariableFeatures(Epithelial))

####T cell
Tcell<-subset(meta_seu,subset=cell_type=='T cell')
Tcell <- NormalizeData(Tcell)
Tcell <- FindVariableFeatures(Tcell,nfeatures=3000)
Tcell <- ScaleData(Tcell,vars.to.regress='percent.mt')
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunBBKNN(Tcell, batch_key = "orig.ident",run_TSNE=TRUE,run_UMAP=TRUE)
Tcell <- FindClusters(Tcell,resolution = 0.7,graph.name='RNA_bbknn')
new.cluster.ids <- c('CD4 Tcm','CD8 CTL','CD4 naive T','CD4 Treg','Transitional T cell','CD8 Tex','CD8 Tem',
                     'CD8 Tem','CD56+ NK','Proliferating T cell','CD4 Tph','CD16+ NK','Unknown','Unknown')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
ncell_layer<-unlist(lapply(Tcell@assays[["RNA"]]@layers, ncol))[147:292]
sketch_ncells <- pmin(
  ncell_layer,
  pmax(round(ncell_layer * 0.1), 100)
)
Tcell <- SketchData(Tcell,
                    ncells = sketch_ncells,method = "Uniform",
                    sketched.assay = "sketch",features = VariableFeatures(Tcell))

####Myeloid
Myeloid<-subset(meta_seu,subset=cell_type %in% c('Myeloid cell','pDC','Mast cell'))
Myeloid <- NormalizeData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid,nfeatures=3000)
Myeloid <- ScaleData(Myeloid,vars.to.regress='percent.mt')
Myeloid <- RunPCA(Myeloid, features = VariableFeatures(object = Myeloid))
Myeloid <- RunBBKNN(Myeloid, batch_key = "orig.ident",run_TSNE=TRUE,run_UMAP=TRUE)
Myeloid <- FindClusters(Myeloid,resolution = 0.9,graph.name='RNA_bbknn')
new.cluster.ids<-c('TREM2+ Macro','Neutrophil','Neutrophil','Neutrophil','Monocytes','MERTK+ Macro','Mast cell','TAM',
                   'cDC2','MKI67+ Macro','Neutrophil','TREM2+ Macro','pDC','mDC','cDC1')
names(new.cluster.ids) <- levels(Myeloid)
Myeloid <- RenameIdents(Myeloid, new.cluster.ids)
Myeloid@meta.data[["cell_subtype"]]<-Myeloid@active.ident

####Fibroblasts
Fib<-subset(meta_seu,subset=cell_type=='Fibroblasts')
Fib <- NormalizeData(Fib)
Fib <- FindVariableFeatures(Fib,nfeatures=3000)
Fib <- ScaleData(Fib,vars.to.regress='percent.mt')
Fib <- RunPCA(Fib, features = VariableFeatures(object = Fib))
Fib <- RunBBKNN(Fib, batch_key = "orig.ident",run_TSNE=TRUE,run_UMAP=TRUE)
Fib <- FindClusters(Fib,resolution = 1,graph.name='RNA_bbknn')

new.cluster.ids <- c('THSD4+ mCAFs','CTHRC1+ myCAFs','CXCL2+ iCAFs','COL14A1+ mCAFs','COL14A1+ mCAFs','COL15A1+ mCAFs',
                     'ACTA2+ myCAFs','SULF1+ myCAFs','CD74+ apCAFs','COL15A1+ mCAFs','MCAM+ vCAFs','MKI67+ ProlifCAFs','Schwann cells')
names(new.cluster.ids) <- levels(Fib)
Fib <- RenameIdents(Fib, new.cluster.ids)
Fib@meta.data[["cell_subtype"]]<-Fib@active.ident

####Endothelial
Endo<-subset(meta_seu,subset=cell_type=='Endothelial cell')
Endo <- NormalizeData(Endo)
Endo <- FindVariableFeatures(Endo,nfeatures=3000)
Endo <- ScaleData(Endo,vars.to.regress='percent.mt')
Endo <- RunPCA(Endo, features = VariableFeatures(object = Endo))
Endo <- RunBBKNN(Endo, batch_key = "orig.ident",run_TSNE=TRUE,run_UMAP=TRUE)
Endo <- FindClusters(Endo,resolution = 0.6,graph.name='RNA_bbknn')
new.cluster.ids<-c('VenEC','VenEC','CapEC','Tip EC','VenEC','ArtEC','CapEC','Unknown','Mural cell','LymEC')
names(new.cluster.ids) <- levels(Endo)
Endo <- RenameIdents(Endo, new.cluster.ids)
Endo@meta.data[["cell_subtype"]]<-Endo@active.ident
