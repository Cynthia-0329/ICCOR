# Purpose: Process and annotate cervical cancer single-cell RNA-seq datasets.
# Input: Dataset-specific count matrices, metadata, and intermediate Seurat objects listed in each section.
# Output: Dataset-specific annotated Seurat objects
#
# Each dataset section preserves its original preprocessing, QC thresholds,
# integration method, clustering parameters, manual annotations, and output calls.

project_dir <- Sys.getenv("ICCOR_PROJECT_DIR", unset = "ICCOR")

# Dataset index
# 01  Cao et al., EMBO Journal (2023)
# 02  Guo et al., Clinical and Translational Medicine (2023)
# 03  Li et al., Cancer Cell International (2025)
# 04  Li et al., Communications Biology (2022), E-MTAB-11948
# 05  Li et al., Frontiers in Immunology (2022), E-MTAB-12305
# 06  Li et al., Journal of Medical Virology (2023), S-BSST1035
# 07  Lin et al., EBioMedicine (2023)
# 08  Qu et al., Cancer Communications (2023), GSE197461
# 09  Zhang et al., EBioMedicine (2023)
# 10  Sandoval et al., Cancer Research (2026), GSE297041
# 11  Yuan et al., Frontiers in Immunology (2026), GSE308792
# 12  Wu et al., Communications Biology (2025), E-MTAB-15983
# 13  Hyeon et al., Molecular Cancer (2025), GSE279998
# 14  Liu et al., Journal of Experimental and Clinical Cancer Research (2023), SCP1950
# 15  Dai et al., Cell Reports Medicine (2024), GSE236738
# 16  Peng et al., eLife (2025), SRP567748
# 17  Cao et al., Journal for ImmunoTherapy of Cancer (2025), treatment cohort
# =============================================================================
# Dataset 01: Cao et al., EMBO Journal (2023)
# =============================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(DoubletFinder)
library(clustree)
library(glmGamPoi)
# Sample integration, clustering, and annotation
mat <- list.dirs(file.path(project_dir, "Dataset01"), recursive = FALSE)
sceList <- lapply(mat,function(x){
  sce <- CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA > 1000 & nCount_RNA >500 & percent.mt < 15)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 0.5)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  nExp_poi <- round(0.075 *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:15, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  sceList[[i]]<-sce
}
for (i in seq_along(CC.sce)) {
  CC.sce[[i]][["pANN"]]<-CC.sce[[i]]@meta.data %>% select(contains('pANN'))
  CC.sce[[i]][["DF.classifications"]]<-CC.sce[[i]]@meta.data %>% select(contains('DF.classifications'))
}
CC.sce <- merge(CC.sce[[1]],
                y = CC.sce[-1],
                add.cell.ids = paste0("Sample", seq_along(sceList)))
CC.sce <- subset(CC.sce, subset = DF.classifications== "Singlet")
CC.sce <- NormalizeData(CC.sce, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce <- FindVariableFeatures(CC.sce, selection.method = "vst", nfeatures = 3000)
CC.sce <- ScaleData(CC.sce, vars.to.regress = "percent.mt")
CC.sce <- RunPCA(CC.sce, features = VariableFeatures(object = CC.sce))
CC.sce <- IntegrateLayers(
  object = CC.sce, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce)
)
CC.sce <- FindNeighbors(CC.sce, dims = 1:30,reduction="cca")
CC.sce <- FindClusters(CC.sce, resolution = 0.5)
CC.sce <- RunUMAP(CC.sce,reduction="cca", dims = 1:30)
CC.sce <- RunTSNE(CC.sce,reduction="cca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('T/NK','Fibroblasts','Epithelial','Fibroblasts','T/NK','Fibroblasts','Endothelial',
                     'Epithelial','Myeloid','T/NK',"B","T/NK",'T/NK','Plasma','Epithelial',
                     'Mast','Myeloid','Fibroblasts','T/NK','Epithelial','Fibroblasts')
names(new.cluster.ids) <- levels(CC.sce)
CC.sce <- RenameIdents(CC.sce, new.cluster.ids)
CC.sce@meta.data[["cell_type"]]<-CC.sce@active.ident

CC.sce.epi<-subset(CC.sce,cell_type=='Epithelial')
CC.sce.epi <- NormalizeData(CC.sce.epi, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce.epi <- FindVariableFeatures(CC.sce.epi, selection.method = "vst", nfeatures = 3000)
CC.sce.epi <- ScaleData(CC.sce.epi, vars.to.regress = "percent.mt")
CC.sce.epi <- RunPCA(CC.sce.epi, features = VariableFeatures(object = CC.sce.epi))
CC.sce.epi <- IntegrateLayers(
  object = CC.sce.epi, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce.epi)
)
CC.sce.epi <- FindNeighbors(CC.sce.epi, dims = 1:30,reduction="cca")
CC.sce.epi <- FindClusters(CC.sce.epi, resolution = 0.6)
CC.sce.epi <- RunUMAP(CC.sce.epi,reduction="cca", dims = 1:30)
CC.sce.epi <- RunTSNE(CC.sce.epi,reduction="cca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('EP0_MUC5B','EP1_KRT6A','EP1_KRT6A','EP4_EPCAM','EP1_KRT6A','EP1_KRT6A',
                     'EP2_POSTN','EP3_MKI67','EP1_KRT6A',"EP1_KRT6A","EP0_MUC5B",'Unknown','EP3_MKI67',
                     'EP0_MUC5B','EP1_KRT6A','Unknown','EP4_EPCAM','EP3_MKI67','EP1_KRT6A')
names(new.cluster.ids) <- levels(CC.sce.epi)
CC.sce.epi <- RenameIdents(CC.sce.epi, new.cluster.ids)
CC.sce.epi@meta.data[["cell_subtype"]]<-CC.sce.epi@active.ident

CC.sce.T<-subset(CC.sce,cell_type=='T/NK')
CC.sce.T <- NormalizeData(CC.sce.T, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce.T <- FindVariableFeatures(CC.sce.T, selection.method = "vst", nfeatures = 3000)
CC.sce.T <- ScaleData(CC.sce.T, vars.to.regress = "percent.mt")
CC.sce.T <- RunPCA(CC.sce.T, features = VariableFeatures(object = CC.sce.T))
CC.sce.T <- IntegrateLayers(
  object = CC.sce.T, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce.T)
)
CC.sce.T <- FindNeighbors(CC.sce.T, dims = 1:30,reduction="cca")
CC.sce.T <- FindClusters(CC.sce.T, resolution = 1)
CC.sce.T <- RunUMAP(CC.sce.T,reduction="cca", dims = 1:30)
CC.sce.T <- RunTSNE(CC.sce.T,reduction="cca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('Tex_HAVCR2','Tex_HAVCR2','Treg_FOXP3','Naive_CD4_CCR7','Tcyto_CD8_GZMK','Tm_CD8_IL7R','Tcyto_CD8_GZMK',
                     'Tprol_MKI67','NK_FCGR3A','Th17_CD4_IL17A',"Tcyto_CD8_GZMK","Tcyto_CD8_GZMK",'NK_KLRC1','Treg_FOXP3','Tfh_CD4_CXCL13',
                     'Unknown','Unknown','Tex_HAVCR2','Unknown','Unknown','Tprol_MKI67')
names(new.cluster.ids) <- levels(CC.sce.T)
CC.sce.T <- RenameIdents(CC.sce.T, new.cluster.ids)
CC.sce.T@meta.data[["cell_subtype"]]<-CC.sce.T@active.ident

CC.sce.B<-subset(CC.sce,cell_type=='B')
CC.sce.B <- NormalizeData(CC.sce.B, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce.B <- FindVariableFeatures(CC.sce.B, selection.method = "vst", nfeatures = 3000)
CC.sce.B <- ScaleData(CC.sce.B, vars.to.regress = "percent.mt")
CC.sce.B <- RunPCA(CC.sce.B, features = VariableFeatures(object = CC.sce.B))
CC.sce.B <- IntegrateLayers(
  object = CC.sce.B, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce.B)
)
CC.sce.B <- FindNeighbors(CC.sce.B, dims = 1:30,reduction="cca")
CC.sce.B <- FindClusters(CC.sce.B, resolution = 0.6)
CC.sce.B <- RunUMAP(CC.sce.B,reduction="cca", dims = 1:30)
CC.sce.B <- RunTSNE(CC.sce.B,reduction="ca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('B0_ABC_TNFRSF13B','PC_IGHA1','PC_IGHG4','B0_ABC_TNFRSF13B','B2_TC_MKI67',
                     'B1_MBC_S1PR1','B0_ABC_TNFRSF13B','B3_GCB_NEIL1','PC_IGHA1','PC_IGHA1')
names(new.cluster.ids) <- levels(CC.sce.B)
CC.sce.B <- RenameIdents(CC.sce.B, new.cluster.ids)
CC.sce.B@meta.data[["cell_subtype"]]<-CC.sce.B@active.ident


CC.sce.myeloid<-subset(CC.sce,cell_type=='Myeloid')
CC.sce.myeloid <- NormalizeData(CC.sce.myeloid, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce.myeloid <- FindVariableFeatures(CC.sce.myeloid, selection.method = "vst", nfeatures = 3000)
CC.sce.myeloid <- ScaleData(CC.sce.myeloid, vars.to.regress = "percent.mt")
CC.sce.myeloid <- RunPCA(CC.sce.myeloid, features = VariableFeatures(object = CC.sce.myeloid))
CC.sce.myeloid <- IntegrateLayers(
  object = CC.sce.myeloid, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce.myeloid)
)
CC.sce.myeloid <- FindNeighbors(CC.sce.myeloid, dims = 1:30,reduction="cca")
CC.sce.myeloid <- FindClusters(CC.sce.myeloid, resolution = 1)
CC.sce.myeloid <- RunUMAP(CC.sce.myeloid,reduction="cca", dims = 1:30)
CC.sce.myeloid <- RunTSNE(CC.sce.myeloid,reduction="cca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('Macro_C1QC','Mono_FCN1','DC_CD1C','Macro_TRAC','Mono_FCN1','Macro_C1QC','Macro_C1QC','DC_GPR183','Macro_C1QC',
                     'DC_CD1C','DC_LAMP3','Macro_C1QC','Mono_FCN1','DC_CLEC9A','Mono_FCN1','DC_LAMP3','Macro_C1QC','DC_GPR183')
names(new.cluster.ids) <- levels(CC.sce.myeloid)
CC.sce.myeloid <- RenameIdents(CC.sce.myeloid, new.cluster.ids)
CC.sce.myeloid@meta.data[["cell_subtype"]]<-CC.sce.myeloid@active.ident


CC.sce.fib<-subset(CC.sce,cell_type=='Fibroblasts')
CC.sce.fib <- NormalizeData(CC.sce.fib, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce.fib <- FindVariableFeatures(CC.sce.fib, selection.method = "vst", nfeatures = 3000)
CC.sce.fib <- ScaleData(CC.sce.fib, vars.to.regress = "percent.mt")
CC.sce.fib <- RunPCA(CC.sce.fib, features = VariableFeatures(object = CC.sce.fib))
CC.sce.fib <- IntegrateLayers(
  object = CC.sce.fib, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce.fib)
)
CC.sce.fib <- FindNeighbors(CC.sce.fib, dims = 1:30,reduction="cca")
CC.sce.fib <- FindClusters(CC.sce.fib, resolution = 0.8)
CC.sce.fib <- RunUMAP(CC.sce.fib,reduction="cca", dims = 1:30)
CC.sce.fib <- RunTSNE(CC.sce.fib,reduction="cca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('PVC0_MCAM','antiCAF_ID2','eCAF_DCN','iCAF_CXCL14','eCAF_DCN','mCAF_ACTG2','PVC1_ACTA2',
                     'PVC1_ACTA2','iCAF_CHI3L1','eCAF_DCN','antiCAF_ID2','Unknown','iCAF_CXCL14')
names(new.cluster.ids) <- levels(CC.sce.fib)
CC.sce.fib <- RenameIdents(CC.sce.fib, new.cluster.ids)
CC.sce.fib@meta.data[["cell_subtype"]]<-CC.sce.fib@active.ident

# =============================================================================
# Dataset 02: Guo et al., Clinical and Translational Medicine (2023)
# =============================================================================
# Sample integration, clustering, and annotation
mat <- list.dirs(file.path(project_dir, "Dataset02"), recursive = FALSE)
sceList <- lapply(mat,function(x){
  sce <- CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA <= mean(sce@meta.data[["nFeature_RNA"]])+2*sd(sce@meta.data[["nFeature_RNA"]]) & nFeature_RNA >= mean(sce@meta.data[["nFeature_RNA"]])-2*sd(sce@meta.data[["nFeature_RNA"]]) & nCount_RNA >= mean(sce@meta.data[["nCount_RNA"]])-2*sd(sce@meta.data[["nCount_RNA"]]) & nCount_RNA <= mean(sce@meta.data[["nCount_RNA"]])+2*sd(sce@meta.data[["nCount_RNA"]]) & percent.mt <= 10)
  sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  nExp_poi <- round(0.075 *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:15, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
  sceList[[i]]<-sce
}
seu_merged<-sceList
seu_merged <- merge(seu_merged[[1]],
                y = seu_merged[-1],
                add.cell.ids = paste0("Sample", seq_along(sceList)))
pann_columns <- grep("pANN_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-pann_columns]
classification_columns <- grep("DF.classifications_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-classification_columns]
seu_merged <- subset(seu_merged, DF.classifications== "Singlet")
seu_merged <- SCTransform(seu_merged, vars.to.regress = "percent.mt", verbose = FALSE)
seu_merged <- RunPCA(seu_merged, features = VariableFeatures(object = seu_merged))
seu_merged <- RunHarmony(seu_merged,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merged <- FindNeighbors(seu_merged, dims = 1:20,reduction="harmony")
seu_merged <- FindClusters(seu_merged, resolution = 0.3)
seu_merged <- RunUMAP(seu_merged,reduction="harmony", dims = 1:20)
seu_merged <- RunTSNE(seu_merged,reduction="harmony", dims = 1:20)
new.cluster.ids <- c('NK_T cells','Neutrophils','NK_T cells','NK_T cells','epithelial cells','NK_T cells','Fibroblasts','epithelial cells',
                     'myeloid cells','epithelial cells','NK_T cells','Mast cells','Plasma cells','myeloid cells','endothelial cells','Plasma cells','epithelial cells',
                     'Plasma cells','B cells','Plasma cells','epithelial cells','smooth muscle cells','epithelial cells',
                     'epithelial cells','dendritic cells','NK_T cells','epithelial cells')
names(new.cluster.ids) <- levels(seu_merged)
seu_merged <- RenameIdents(seu_merged, new.cluster.ids)
seu_merged@meta.data[["Cell_type"]]<-as.character(seu_merged@active.ident)

epithelial<-subset(seu_merged,Cell_type =='epithelial cells')
epithelial <- SCTransform(epithelial, vars.to.regress = "percent.mt", verbose = FALSE)
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
epithelial <- RunHarmony(epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
epithelial <- FindNeighbors(epithelial, dims = 1:30,reduction="harmony")
epithelial <- FindClusters(epithelial, resolution = 0.3)
epithelial <- RunUMAP(epithelial,reduction="harmony", dims = 1:30)
epithelial <- RunTSNE(epithelial,reduction="harmony", dims = 1:30)

# T-cell subclustering
t_nk<-subset(seu_merged,Cell_type =='NK_T cells')
t_nk <- SCTransform(t_nk, vars.to.regress = "percent.mt", verbose = FALSE)
t_nk <- RunPCA(t_nk, features = VariableFeatures(object = t_nk))
t_nk <- RunHarmony(t_nk,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
t_nk <- FindNeighbors(t_nk, dims = 1:30,reduction="harmony")
t_nk <- FindClusters(t_nk, resolution = 0.7)
t_nk <- RunUMAP(t_nk,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('CD4','CD8','CD8','CD8','CD4','CD4','CD8','CD4','CD8',
                     'CD4','CD4','CD8','CD8','CD8','CD8','CD8','CD4','CD4')
names(new.cluster.ids) <- levels(t_nk)
t_nk <- RenameIdents(t_nk, new.cluster.ids)
t_nk@meta.data[["Cell_subtype"]]<-as.character(t_nk@active.ident)

#CD4
cd4_t<-subset(t_nk,Cell_subtype =='CD4')
cd4_t <- RunPCA(cd4_t, features = VariableFeatures(object = cd4_t))
cd4_t <- RunHarmony(cd4_t,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
cd4_t <- FindNeighbors(cd4_t, dims = 1:30,reduction="harmony")
cd4_t <- FindClusters(cd4_t, resolution = 0.3)
cd4_t <- RunUMAP(cd4_t,reduction="harmony", dims = 1:30)
cd4_t <- RunTSNE(cd4_t,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('Tem','Treg','Tem','Naive','Tem','Th1-like','Tem','Th17','Treg')
names(new.cluster.ids) <- levels(cd4_t)
cd4_t <- RenameIdents(cd4_t, new.cluster.ids)
cd4_t@meta.data[["Cell_subtype"]]<-as.character(cd4_t@active.ident)

#CD8
cd8_t<-subset(t_nk,Cell_subtype =='CD8')
cd8_t <- RunPCA(cd8_t, features = VariableFeatures(object = cd8_t))
cd8_t <- RunHarmony(cd8_t,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
cd8_t <- FindNeighbors(cd8_t, dims = 1:30,reduction="harmony")
cd8_t <- FindClusters(cd8_t, resolution = 0.7)
cd8_t <- RunUMAP(cd8_t,reduction="harmony", dims = 1:30)
cd8_t <- RunTSNE(cd8_t,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('Tem','Tem','Trm','Tem','TemRA','Tem','Trm','Tem','Trm','IEL','Naive','MAIT','IEL','Tex')
names(new.cluster.ids) <- levels(cd8_t)
cd8_t <- RenameIdents(cd8_t, new.cluster.ids)
cd8_t@meta.data[["Cell_subtype"]]<-as.character(cd8_t@active.ident)
# Myeloid-cell subclustering
myeloid<-subset(seu_merged,Cell_type =='myeloid cells')
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid))
myeloid <- RunHarmony(myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
myeloid <- FindNeighbors(myeloid, dims = 1:30,reduction="harmony")
myeloid <- FindClusters(myeloid, resolution = 0.8)
myeloid <- RunUMAP(myeloid,reduction="harmony", dims = 1:30)
myeloid <- RunTSNE(myeloid,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('Monocytes','Macrophages','cDC2','cDC2','Macrophages','cDC2','Macrophages','Monocytes','Macrophages',
                     'cDC1','Macrophages','Macrophages','cDC2','pDC')
names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)
myeloid@meta.data[["Cell_subtype"]]<-as.character(myeloid@active.ident)

# =============================================================================
# Dataset 03: Li et al., Cancer Cell International (2025)
# Original script: single-cell/Li_Cancer.Cell.Int_2025/Zhengzhou_for_DB.R
# =============================================================================

# Sample integration, clustering, and annotation
mat <- list.dirs(file.path(project_dir, "Dataset03"), recursive = FALSE)
sceList <- lapply(mat,function(x){
  sce <- CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
CC.sce<-list()
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(sce@assays$RNA)) # Match hemoglobin genes to expression-matrix rows
  HB.genes <- rownames(sce@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]  # Keep matched genes and remove missing values
  sce[["percent.HB"]] <- PercentageFeatureSet(sce, features=HB.genes)
  if (i %in% 1:2){
    sce <- subset(sce, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & nCount_RNA >1000 & nCount_RNA < 25000 & percent.mt < 15 & percent.HB < 5)
  } else if (i == 3){
    sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & nCount_RNA >1000 & nCount_RNA < 50000 & percent.mt < 15 & percent.HB < 5)
  } else if (i == 4){
    sce <- subset(sce, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & nCount_RNA >500 & nCount_RNA < 40000 & percent.mt < 15 & percent.HB < 5)
  } else if (i == 5){
    sce <- subset(sce, subset = nFeature_RNA > 600 & nFeature_RNA < 5500 & nCount_RNA >1000 & nCount_RNA < 40000 & percent.mt < 15 & percent.HB < 5)
  } else if (i == 6){
    sce <- subset(sce, subset = nFeature_RNA > 600 & nFeature_RNA < 8000 & nCount_RNA >1000 & nCount_RNA < 50000 & percent.mt < 25 & percent.HB < 5)
  } else if (i == 7){
    sce <- subset(sce, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & nCount_RNA >1000 & nCount_RNA < 40000 & percent.mt < 15 & percent.HB < 5)
  } else if (i == 8){
    sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA >1000 & nCount_RNA < 30000 & percent.mt < 15 & percent.HB < 5)
  }
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  CC.sce[[i]]<-sce
}

CC.sce <- merge(CC.sce[[1]],
                y = CC.sce[-1],
                add.cell.ids = paste0("Sample", seq_along(sceList)))
CC.sce <- NormalizeData(CC.sce, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce <- FindVariableFeatures(CC.sce, selection.method = "vst", nfeatures = 2000)
CC.sce <- ScaleData(CC.sce, vars.to.regress = "percent.mt")
CC.sce <- RunPCA(CC.sce, features = VariableFeatures(object = CC.sce))
CC.sce <- IntegrateLayers(
  object = CC.sce, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce)
)
CC.sce <- FindNeighbors(CC.sce, dims = 1:30,reduction="cca")
CC.sce <- FindClusters(CC.sce, resolution = 0.7)
CC.sce <- RunUMAP(CC.sce,reduction="cca", dims = 1:30)
CC.sce <- RunTSNE(CC.sce,reduction="cca", dims = 1:30)


# Manual cluster annotation
new.cluster.ids <- c('Epithelial cells','Fibroblasts','Epithelial cells','T cells','T cells','Epithelial cells','Endothelial cells',
                     'Fibroblasts','T cells','Myeloid cells','Epithelial cells','Fibroblasts','Fibroblasts','Fibroblasts','B cells',
                     'T cells','T cells','Epithelial cells','Plasma cells','Myeloid cells','Mast cells','Myeloid cells','Endothelial cells','Epithelial cells')
names(new.cluster.ids) <- levels(CC.sce)
CC.sce <- RenameIdents(CC.sce, new.cluster.ids)
CC.sce@meta.data[["cell_type"]]<-CC.sce@active.ident

# Epithelial-cell subclustering
CC.epi <- subset(CC.sce,cell_type=='Epithelial cells')
CC.epi <- NormalizeData(CC.epi, normalization.method = "LogNormalize", scale.factor = 10000)
CC.epi <- FindVariableFeatures(CC.epi, selection.method = "vst", nfeatures = 2000)
CC.epi <- ScaleData(CC.epi, vars.to.regress = "percent.mt")
CC.epi <- RunPCA(CC.epi, features = VariableFeatures(object = CC.epi))
CC.epi <- IntegrateLayers(
  object = CC.epi, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.epi)
)
CC.epi <- FindNeighbors(CC.epi, dims = 1:30,reduction="cca")
CC.epi <- FindClusters(CC.epi, resolution = 0.6)
CC.epi <- RunUMAP(CC.epi,reduction="cca", dims = 1:30)
CC.epi <- RunTSNE(CC.epi,reduction="cca", dims = 1:30)

# Manual cluster annotation
new.cluster.ids <- c('C2-DST-NEpis','C6-TFF3-IAEpis','C1-IGFL1-Epis','C2-DST-NEpis','C2-DST-NEpis','C5-PCLAF-TAEpis','C4-KRTDAP-IAEpis',
                     'C3-MUC20-Tu','C5-PCLAF-TAEpis','C2-DST-NEpis','C8-NEURL1B-TAEpis','C2-DST-NEpis','C7-CENPF-TAEpis','C7-CENPF-TAEpis','C6-TFF3-IAEpis')
names(new.cluster.ids) <- levels(CC.epi)
CC.epi <- RenameIdents(CC.epi, new.cluster.ids)
CC.epi@meta.data[["cell_subtype"]]<-CC.epi@active.ident

# Fibroblast subclustering
CC.fib <- subset(CC.sce,cell_type=='Fibroblasts')
CC.fib <- NormalizeData(CC.fib, normalization.method = "LogNormalize", scale.factor = 10000)
CC.fib <- FindVariableFeatures(CC.fib, selection.method = "vst", nfeatures = 2000)
CC.fib <- ScaleData(CC.fib, vars.to.regress = "percent.mt")
CC.fib <- RunPCA(CC.fib, features = VariableFeatures(object = CC.fib))
CC.fib <- IntegrateLayers(
  object = CC.fib, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.fib)
)
CC.fib <- FindNeighbors(CC.fib, dims = 1:30,reduction="cca")
CC.fib <- FindClusters(CC.fib, resolution = 1)
CC.fib <- RunUMAP(CC.fib,reduction="cca", dims = 1:30)
CC.fib <- RunTSNE(CC.fib,reduction="cca", dims = 1:30)


# Manual cluster annotation
new.cluster.ids <- c('C2-MMP11-CAFs','C1-SPER4-IAFs','C4-RGS5-pericytes','C3-MUSTN1-myofibro','C1-SPER4-IAFs','C1-SPER4-IAFs',
                     'C2-MMP11-CAFs','C2-MMP11-CAFs','C5-DES-myofibro','C1-SPER4-IAFs','C1-SPER4-IAFs','C1-SPER4-IAFs',
                     'C3-MUSTN1-myofibro','C1-SPER4-IAFs','C3-MUSTN1-myofibro','C2-MMP11-CAFs','C1-SPER4-IAFs','C2-MMP11-CAFs')
names(new.cluster.ids) <- levels(CC.fib)
CC.fib <- RenameIdents(CC.fib, new.cluster.ids)
CC.fib@meta.data[["cell_subtype"]]<-CC.fib@active.ident

# Endothelial-cell subclustering
CC.ec <- subset(CC.sce,cell_type=='Endothelial cells')
CC.ec <- NormalizeData(CC.ec, normalization.method = "LogNormalize", scale.factor = 10000)
CC.ec <- FindVariableFeatures(CC.ec, selection.method = "vst", nfeatures = 2000)
CC.ec <- ScaleData(CC.ec, vars.to.regress = "percent.mt")
CC.ec <- RunPCA(CC.ec, features = VariableFeatures(object = CC.ec))
CC.ec <- IntegrateLayers(
  object = CC.ec, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.ec)
)
CC.ec <- FindNeighbors(CC.ec, dims = 1:30,reduction="cca")
CC.ec <- FindClusters(CC.ec, resolution = 0.9)
CC.ec <- RunUMAP(CC.ec,reduction="cca", dims = 1:30)
CC.ec <- RunTSNE(CC.ec,reduction="cca", dims = 1:30)


# Manual cluster annotation
new.cluster.ids <- c('C2-EDNRB-capillary ECs','C1-ACKR1-Venous ECs','C1-ACKR1-Venous ECs','C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs',
                     'C4-FBLBN5-arterial ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C2-EDNRB-capillary ECs','C3-CCL21-lymphatic ECs')
new.cluster.ids <- c('C1-ACKR1-Venous ECs','C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs','C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs',
                     'C2-EDNRB-capillary ECs','C4-FBLBN5-arterial ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C2-EDNRB-capillary ECs','C3-CCL21-lymphatic ECs')
names(new.cluster.ids) <- levels(CC.ec)
CC.ec <- RenameIdents(CC.ec, new.cluster.ids)
CC.ec@meta.data[["cell_subtype"]]<-CC.ec@active.ident

# Myeloid-cell subclustering
CC.mye <- subset(CC.sce,cell_type=='Myeloid cells')
CC.mye <- NormalizeData(CC.mye, normalization.method = "LogNormalize", scale.factor = 10000)
CC.mye <- FindVariableFeatures(CC.mye, selection.method = "vst", nfeatures = 2000)
CC.mye <- ScaleData(CC.mye, vars.to.regress = "percent.mt")
CC.mye <- RunPCA(CC.mye, features = VariableFeatures(object = CC.mye))
CC.mye <- IntegrateLayers(
  object = CC.mye, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.mye)
)
CC.mye <- FindNeighbors(CC.mye, dims = 1:30,reduction="cca")
CC.mye <- FindClusters(CC.mye, resolution = 0.6)
CC.mye <- RunUMAP(CC.mye,reduction="cca", dims = 1:30)
CC.mye <- RunTSNE(CC.mye,reduction="cca", dims = 1:30)


# Manual cluster annotation
new.cluster.ids <- c('C7-ISG15-Neus','C4-CD163-TAMs','C2-S100A8-IANs','C5-CD1C-cDC2','C1-C1QA-Macro','C6-CXCL8-TANs',
                     'C5-CD1C-cDC2','C5-CD1C-cDC2','C3-CXCR4-TANs','C8-LAMP3-cDC3','C1-C1QA-Macro','C1-C1QA-Macro','C8-LAMP3-cDC3','C5-CD1C-cDC2')
names(new.cluster.ids) <- levels(CC.mye)
CC.mye <- RenameIdents(CC.mye, new.cluster.ids)
CC.mye@meta.data[["cell_subtype"]]<-CC.mye@active.ident

# T/NK-cell subclustering
CC.NKT <- subset(CC.sce,cell_type=='T cells')
CC.NKT <- NormalizeData(CC.NKT, normalization.method = "LogNormalize", scale.factor = 10000)
CC.NKT <- FindVariableFeatures(CC.NKT, selection.method = "vst", nfeatures = 2000)
CC.NKT <- ScaleData(CC.NKT, vars.to.regress = "percent.mt")
CC.NKT <- RunPCA(CC.NKT, features = VariableFeatures(object = CC.NKT))
CC.NKT <- IntegrateLayers(
  object = CC.NKT, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.NKT)
)
CC.NKT <- FindNeighbors(CC.NKT, dims = 1:30,reduction="cca")
CC.NKT <- FindClusters(CC.NKT, resolution = 0.8)
CC.NKT <- RunUMAP(CC.NKT,reduction="cca", dims = 1:30)
CC.NKT <- RunTSNE(CC.NKT,reduction="cca", dims = 1:30)


# Manual cluster annotation
new.cluster.ids <- c('CD8-C1-ZNF683-Trm','CD4-C2-IL7R-Tcm','CD8-C3-GZMK-Tem','CD4-C4-FOXP3-Treg','CD8-C1-ZNF683-Trm','CD4-C2-IL7R-Tcm',
                     'CD8-C6-CXCL13-Tex','NK-C5-NKG7','CD8-C8-GPR183-Tcm','DP-C11-HIST1H1B','CD4-C10-CXCL13-Th1','NK-C5-NKG7','DP-C11-HIST1H1B')
names(new.cluster.ids) <- levels(CC.NKT)
CC.NKT <- RenameIdents(CC.NKT, new.cluster.ids)
CC.NKT@meta.data[["cell_subtype"]]<-CC.NKT@active.ident

# =============================================================================
# Dataset 04: Li et al., Communications Biology (2022), E-MTAB-11948
# Original script: single-cell/Li_Commun.Biol_2022/E-MTAB-11948.R
# =============================================================================

# Sample integration, clustering, and annotation
mat <- list.files(file.path(project_dir, "Dataset04"), full.names = TRUE)
sceList <- lapply(mat,function(x){
  counts = read.csv(x)
  counts = subset(counts,X!='0')
  colnames(counts)<-counts[1,]
  counts<-counts[-1,]
  rownames(counts)<-counts[,1]
  counts<-counts[,-1]
  sce <- CreateSeuratObject(counts = counts,
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})

sce.doubletFinder.list<-list()
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA >=200 & nCount_RNA >= 200 & percent.mt <= 10)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  nExp_poi <- round(0.075 *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:15, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  sce.doubletFinder.list[[i]]<-sce
}

seu_merged<-sce.doubletFinder.list
for (i in seq_along(seu_merged)) {
  seu_merged[[i]][["pANN"]]<-seu_merged[[i]]@meta.data %>% select(contains('pANN'))
  seu_merged[[i]][["DF.classifications"]]<-seu_merged[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

seu_merged <- merge(seu_merged[[1]],
                 y = seu_merged[-1],
                 add.cell.ids = paste0("Sample", seq_along(sceList)))
pann_columns <- grep("pANN_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-pann_columns]
classification_columns <- grep("DF.classifications_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-classification_columns]
seu_merged <- subset(seu_merged, DF.classifications== "Singlet")

seu_merged <- NormalizeData(seu_merged, normalization.method = "LogNormalize")
seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 2000)
seu_merged <- ScaleData(seu_merged, vars.to.regress = "percent.mt",features = VariableFeatures(object = seu_merged))
seu_merged <- RunPCA(seu_merged, features = VariableFeatures(object = seu_merged))
seu_merged <- RunHarmony(seu_merged,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merged <- FindNeighbors(seu_merged, dims = 1:30,reduction="harmony")
seu_merged <- FindClusters(seu_merged, resolution = 0.7)
seu_merged <- RunUMAP(seu_merged,reduction="harmony", dims = 1:30)
seu_merged <- RunTSNE(seu_merged,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('Epithelial cells','Fibroblasts','Endothelial cells','Epithelial cells','Fibroblasts','T cells',
                     'Epithelial cells','Smooth muscle cells','Fibroblasts','Smooth muscle cells','Fibroblasts',
                     'Endothelial cells','Neutrophils','Endothelial cells','Macrophages','Epithelial cells','T cells',
                     'Endothelial cells','Epithelial cells','Mast cells','Endothelial cells','B cells','Fibroblasts')
names(new.cluster.ids) <- levels(seu_merged)
seu_merged <- RenameIdents(seu_merged, new.cluster.ids)
seu_merged@meta.data[["Cell_type"]]<-as.character(seu_merged@active.ident)
seu_merged@meta.data[["Cell_type"]]<-factor(seu_merged@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblasts','Smooth muscle cells','Endothelial cells',
                                                  'T cells','Neutrophils','Macrophages','B cells','Mast cells'))
Idents(seu_merged)<-seu_merged@meta.data[["Cell_type"]]


# Epithelial-cell subclustering
epithelial<-subset(seu_merged,Cell_type =='Epithelial cells')
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
epithelial <- FindNeighbors(epithelial, dims = 1:30,reduction="pca")
epithelial <- FindClusters(epithelial, resolution = 0.5)
epithelial <- RunUMAP(epithelial,reduction="pca", dims = 1:30)
epithelial <- RunTSNE(epithelial,reduction="pca", dims = 1:30)


new.cluster.ids <- c('C2','C1','C1','C1','C5','C4','C2',
                     'C6','C2','C7','C3','C6','C4','C5')
names(new.cluster.ids) <- levels(epithelial)
epithelial <- RenameIdents(epithelial, new.cluster.ids)
epithelial@meta.data[["Cell_subtype"]]<-as.character(epithelial@active.ident)
epithelial@meta.data[["Cell_subtype"]]<-factor(epithelial@meta.data[["Cell_subtype"]],
                                         levels=c('C1','C2','C3','C4','C5','C6','C7'))
Idents(epithelial)<-epithelial@meta.data[["Cell_subtype"]]


# T-cell subclustering
t_cell<-subset(seu_merged,Cell_type =='T cells')
t_cell <- RunPCA(t_cell, features = VariableFeatures(object = t_cell))
t_cell <- RunHarmony(t_cell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
t_cell <- FindNeighbors(t_cell, dims = 1:30,reduction="harmony")
t_cell <- FindClusters(t_cell, resolution = 1)
t_cell <- RunUMAP(t_cell,reduction="harmony", dims = 1:30)
t_cell <- RunTSNE(t_cell,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('CXCR4+CD8','PDCD1+CD8','IL7R+Tm','IL7R+Tm','NK','IL7R+Tm','MKI67+CD8','Unknown','TNFRSF4+Treg','NK')
names(new.cluster.ids) <- levels(t_cell)
t_cell <- RenameIdents(t_cell, new.cluster.ids)
t_cell@meta.data[["Cell_subtype"]]<-as.character(t_cell@active.ident)
t_cell@meta.data[["Cell_subtype"]]<-factor(t_cell@meta.data[["Cell_subtype"]],
                                            levels=c('CXCR4+CD8','PDCD1+CD8','MKI67+CD8','NK','IL7R+Tm','TNFRSF4+Treg','Unknown'))
Idents(t_cell)<-t_cell@meta.data[["Cell_subtype"]]


# Fibroblast and smooth-muscle-cell subclustering
fibroblast<-subset(seu_merged,Cell_type %in% c('Fibroblasts','Smooth muscle cells'))
fibroblast <- RunPCA(fibroblast, features = VariableFeatures(object = fibroblast))
fibroblast <- RunHarmony(fibroblast,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
fibroblast <- FindNeighbors(fibroblast, dims = 1:30,reduction="harmony")
fibroblast <- FindClusters(object = fibroblast,resolution = 0.5)
fibroblast <- RunUMAP(fibroblast,reduction="harmony", dims = 1:30)
fibroblast <- RunTSNE(fibroblast,reduction="harmony", dims = 1:30)


# =============================================================================
# Dataset 05: Li et al., Frontiers in Immunology (2022), E-MTAB-12305
# Original script: single-cell/Li_Front.Immunol_2022/E-MTAB-12305.R
# =============================================================================

# Sample integration, clustering, and annotation
mat <- list.files(file.path(project_dir, "Dataset05"), full.names = TRUE)
sceList <- lapply(mat,function(x){
  counts = read.csv(x)
  counts = subset(counts,X!='0')
  colnames(counts)<-counts[1,]
  counts<-counts[-1,]
  rownames(counts)<-counts[,1]
  counts<-counts[,-1]
  sce <- CreateSeuratObject(counts = counts,
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})


sce.doubletFinder.list<-list()
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nCount_RNA >= 200 & percent.mt <= 10)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:30,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:30)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  nExp_poi <- round(0.075 *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:15, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  sce.doubletFinder.list[[i]]<-sce
}

seu_merged<-sce.doubletFinder.list
for (i in seq_along(seu_merged)) {
  seu_merged[[i]][["pANN"]]<-seu_merged[[i]]@meta.data %>% select(contains('pANN'))
  seu_merged[[i]][["DF.classifications"]]<-seu_merged[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

seu_merged <- merge(seu_merged[[1]],
                 y = seu_merged[-1],
                 add.cell.ids = paste0("Sample", seq_along(sceList)))
pann_columns <- grep("pANN_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-pann_columns]
classification_columns <- grep("DF.classifications_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-classification_columns]
seu_merged <- subset(seu_merged, DF.classifications== "Singlet")


seu_merged <- NormalizeData(seu_merged, normalization.method = "LogNormalize")
seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 2000)
seu_merged <- ScaleData(seu_merged, vars.to.regress = "percent.mt",features = VariableFeatures(object = seu_merged))
seu_merged <- RunPCA(seu_merged, features = VariableFeatures(object = seu_merged))
seu_merged <- RunHarmony(seu_merged,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merged <- FindNeighbors(seu_merged, dims = 1:30,reduction="harmony")
seu_merged <- FindClusters(seu_merged, resolution = 0.7)
seu_merged <- RunUMAP(seu_merged,reduction="harmony", dims = 1:30)
seu_merged <- RunTSNE(seu_merged,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('Fibroblasts','Epithelial cells','T/NK cells','T/NK cells','Endothelial cells','Fibroblasts',
                     'Epithelial cells','Smooth muscle cells','Epithelial cells','Macrophages','Neutrophils','Epithelial cells',
                     'Endothelial cells','Fibroblasts','Epithelial cells','Fibroblasts','Mast cells','Epithelial cells',
                     'Plasma cells','T/NK cells','B cells','Epithelial cells','Epithelial cells','Endothelial cells',
                     'Macrophages','Epithelial cells','Macrophages','Fibroblasts','Endothelial cells')
names(new.cluster.ids) <- levels(seu_merged)
seu_merged <- RenameIdents(seu_merged, new.cluster.ids)
seu_merged@meta.data[["Cell_type"]]<-as.character(seu_merged@active.ident)
seu_merged@meta.data[["Cell_type"]]<-factor(seu_merged@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblasts','Smooth muscle cells','Endothelial cells',
                                                  'T/NK cells','Neutrophils','Macrophages','B cells','Plasma cells','Mast cells'))
Idents(seu_merged)<-seu_merged@meta.data[["Cell_type"]]


# Epithelial-cell subclustering
epithelial<-subset(seu_merged,Cell_type =='Epithelial cells')
epithelial<-epithelial[rownames(seu_merged),]
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
epithelial <- FindNeighbors(epithelial, dims = 1:30,reduction="pca")
epithelial <- FindClusters(epithelial, resolution = 0.5)
epithelial <- RunUMAP(epithelial,reduction="pca", dims = 1:30)
epithelial <- RunTSNE(epithelial,reduction="pca", dims = 1:30)


# T-cell subclustering
t_nk<-subset(seu_merged,Cell_type =='T/NK cells')
t_nk <- RunPCA(t_nk, features = VariableFeatures(object = t_nk))
t_nk <- RunHarmony(t_nk,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
t_nk <- FindNeighbors(t_nk, dims = 1:30,reduction="harmony")
t_nk <- FindClusters(t_nk, resolution = 1)
t_nk <- RunUMAP(t_nk,reduction="harmony", dims = 1:30)
t_nk <- RunTSNE(t_nk,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('Cytotoxic CD8 T cells','Naive CD4 T cells','Effector NK cells','Reg T cells','Naive CD4 T cells',
                     'Effector memory CD8 T cells','Central memory T cells','Effector memory CD8 T cells','Exhausted CD8 T cells','Effector NK cells','Central memory T cells',
                     'Cytotoxic CD8 T cells','Central memory T cells','Cycling T cells','Naive CD4 T cells','Circulating NK cells')
names(new.cluster.ids) <- levels(t_nk)
t_nk <- RenameIdents(t_nk, new.cluster.ids)
t_nk@meta.data[["Cell_subtype"]]<-as.character(t_nk@active.ident)
Idents(t_nk)<-t_nk@meta.data[["Cell_subtype"]]


# Myeloid-cell subclustering
myeloid<-subset(seu_merged,Cell_type =='Macrophages')
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid))
myeloid <- RunHarmony(myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
myeloid <- FindNeighbors(myeloid, dims = 1:30,reduction="harmony")
myeloid <- FindClusters(myeloid, resolution = 0.8)
myeloid <- RunUMAP(myeloid,reduction="harmony", dims = 1:30)
myeloid <- RunTSNE(myeloid,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('C1-Ma-C1QA','C4-Ma-GZMA','C2-Ma-THBS1','C3-DC2','C2-Ma-THBS1','C1-Ma-C1QA','C8-Ma-SPP1',
                     'C6-Ma-MAFB','C7-Ma-MKI67','C3-DC2','C4-Ma-GZMA','C2-Ma-THBS1','C5-Ma-S100A8','C9-DC1','C6-Ma-MAFB')
names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)
myeloid@meta.data[["Cell_subtype"]]<-as.character(myeloid@active.ident)
myeloid@meta.data[["Cell_subtype"]]<-factor(myeloid@meta.data[["Cell_subtype"]],
                                                levels=c('C1-Ma-C1QA','C2-Ma-THBS1','C4-Ma-GZMA','C5-Ma-S100A8','C6-Ma-MAFB',
                                                         'C7-Ma-MKI67','C8-Ma-SPP1','C3-DC2','C9-DC1'))
Idents(myeloid)<-myeloid@meta.data[["Cell_subtype"]]


# Fibroblast subclustering
fibroblast<-subset(seu_merged,Cell_type == 'Fibroblasts')
fibroblast <- RunPCA(fibroblast, features = VariableFeatures(object = fibroblast))
fibroblast <- RunHarmony(fibroblast,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
fibroblast <- FindNeighbors(fibroblast, dims = 1:30,reduction="harmony")
fibroblast <- RunUMAP(fibroblast,reduction="harmony", dims = 1:30)
fibroblast <- RunTSNE(fibroblast,reduction="harmony", dims = 1:30)


# =============================================================================
# Dataset 06: Li et al., Journal of Medical Virology (2023), S-BSST1035
# Original script: single-cell/Li_J.Med.Virol_2023/S-BSST1035.R
# =============================================================================

# Sample integration, clustering, and annotation
mat <- list.files(file.path(project_dir, "Dataset06"), full.names = TRUE)
sceList <- lapply(mat,function(x){
  counts = read.csv(x)
  counts = subset(counts,X!='0')
  colnames(counts)<-counts[1,]
  counts<-counts[-1,]
  rownames(counts)<-counts[,1]
  counts<-counts[,-1]
  sce <- CreateSeuratObject(counts = counts,
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})


sce.doubletFinder.list<-list()
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA >=200 & nCount_RNA >= 800 & percent.mt <= 20)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  nExp_poi <- round(0.075 *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:15, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  sce.doubletFinder.list[[i]]<-sce
}

seu_merged<-sce.doubletFinder.list
for (i in seq_along(seu_merged)) {
  seu_merged[[i]][["pANN"]]<-seu_merged[[i]]@meta.data %>% select(contains('pANN'))
  seu_merged[[i]][["DF.classifications"]]<-seu_merged[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

seu_merged <- merge(seu_merged[[1]],
                 y = seu_merged[-1],
                 add.cell.ids = paste0("Sample", seq_along(sceList)))
pann_columns <- grep("pANN_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-pann_columns]
classification_columns <- grep("DF.classifications_0.25",colnames(seu_merged@meta.data))
seu_merged@meta.data <- seu_merged@meta.data[,-classification_columns]
seu_merged <- subset(seu_merged, DF.classifications== "Singlet")

seu_merged <- NormalizeData(seu_merged, normalization.method = "LogNormalize")
seu_merged <- FindVariableFeatures(seu_merged, selection.method = "vst", nfeatures = 2000)
seu_merged <- ScaleData(seu_merged, vars.to.regress = "percent.mt",features = VariableFeatures(object = seu_merged))
seu_merged <- RunPCA(seu_merged, features = VariableFeatures(object = seu_merged))
seu_merged <- RunHarmony(seu_merged,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merged <- FindNeighbors(seu_merged, dims = 1:20,reduction="harmony")
seu_merged <- FindClusters(seu_merged, resolution = 0.8)
seu_merged <- RunUMAP(seu_merged,reduction="harmony", dims = 1:20)
seu_merged <- RunTSNE(seu_merged,reduction="harmony", dims = 1:20)


new.cluster.ids <- c('Epithelial cells','Epithelial cells','Epithelial cells','Epithelial cells','NK/T cells','Epithelial cells',
                     'NK/T cells','Epithelial cells','Myeloid cells','Epithelial cells','Epithelial cells','Fibroblasts',
                     'Fibroblasts','B cells','NK/T cells','Neutrophils','Myeloid cells','NK/T cells','Plasma cells',
                     'Epithelial cells','Mast cells','Epithelial cells','Endothelial cells','Epithelial cells')
names(new.cluster.ids) <- levels(seu_merged)
seu_merged <- RenameIdents(seu_merged, new.cluster.ids)
seu_merged@meta.data[["Cell_type"]]<-as.character(seu_merged@active.ident)
seu_merged@meta.data[["Cell_type"]]<-factor(seu_merged@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblasts','Endothelial cells','NK/T cells',
                                                  'Neutrophils','Myeloid cells','B cells','Plasma cells','Mast cells'))
Idents(seu_merged)<-seu_merged@meta.data[["Cell_type"]]


# Epithelial-cell subclustering
epithelial<-subset(seu_merged,Cell_type =='Epithelial cells')
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
epithelial <- FindNeighbors(epithelial, dims = 1:30,reduction="pca")
epithelial <- FindClusters(epithelial, resolution = 0.3)
epithelial <- RunUMAP(epithelial,reduction="pca", dims = 1:30)
epithelial <- RunTSNE(epithelial,reduction="pca", dims = 1:30)


new.cluster.ids <- c('C0','C3','C5','C1','C4','C6','C2',
                     'C9','C8','C3','C5','C0','C1','C7')
names(new.cluster.ids) <- levels(epithelial)
epithelial <- RenameIdents(epithelial, new.cluster.ids)
epithelial@meta.data[["Cell_subtype"]]<-as.character(epithelial@active.ident)
epithelial@meta.data[["Cell_subtype"]]<-factor(epithelial@meta.data[["Cell_subtype"]],
                                            levels=c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'))
Idents(epithelial)<-epithelial@meta.data[["Cell_subtype"]]


# T-cell subclustering
t_cell<-subset(seu_merged,Cell_type =='NK/T cells')
t_cell <- RunPCA(t_cell, features = VariableFeatures(object = t_cell))
t_cell <- RunHarmony(t_cell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
t_cell <- FindNeighbors(t_cell, dims = 1:30,reduction="harmony")
t_cell <- FindClusters(t_cell, resolution = 1)
t_cell <- RunUMAP(t_cell,reduction="harmony", dims = 1:30)
t_cell <- RunTSNE(t_cell,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('Tissue resident CD8 T cells','Central memory CD4 T cells','CD160+ NK cells','Naive CD4 T cells','Tissue resident CD8 T cells','Treg CD4 T cells','Central memory CD4 T cells',
                     'Effector memory CD8 T cells','Proliferative cells','Exhausted CD8 T cells','Central memory CD4 T cells','Central memory CD4 T cells','Tissue resident NK cells')
names(new.cluster.ids) <- levels(t_cell)
t_cell <- RenameIdents(t_cell, new.cluster.ids)
t_cell@meta.data[["Cell_subtype"]]<-as.character(t_cell@active.ident)
Idents(t_cell)<-t_cell@meta.data[["Cell_subtype"]]


# Fibroblast subclustering
fibroblast<-subset(seu_merged,Cell_type == 'Fibroblasts')
fibroblast <- RunPCA(fibroblast, features = VariableFeatures(object = fibroblast))
fibroblast <- RunHarmony(fibroblast,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
fibroblast <- FindNeighbors(fibroblast, dims = 1:30,reduction="harmony")
fibroblast <- FindClusters(object = fibroblast,resolution = 1)
fibroblast <- RunUMAP(fibroblast,reduction="harmony", dims = 1:30)
fibroblast <- RunTSNE(fibroblast,reduction="harmony", dims = 1:30)


# Myeloid-cell subclustering
myeloid<-subset(seu_merged,Cell_type =='Myeloid cells')
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid))
myeloid <- RunHarmony(myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
myeloid <- FindNeighbors(myeloid, dims = 1:30,reduction="harmony")
myeloid <- FindClusters(myeloid, resolution = 1)
myeloid <- RunUMAP(myeloid,reduction="harmony", dims = 1:30)
myeloid <- RunTSNE(myeloid,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('C1-Ma-FCGR3A','C1-Ma-FCGR3A','C1-Ma-FCGR3A','C7-Ma6','C7-Ma6','C3-CD1C-cDC2','C4-Ma4',
                     'C0-Ma-CD74','C8-Ma7','Cycling cells','C1-Ma-FCGR3A','C2-Ma-CD68','C9-LAMP3 DCs')
names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)
myeloid@meta.data[["Cell_subtype"]]<-as.character(myeloid@active.ident)
#myeloid@meta.data[["Cell_subtype"]]<-factor(myeloid@meta.data[["Cell_subtype"]],
#                                                levels=c('MAIT','Trm','Naive','Tex','IEL','Tem','TemRA'))
# =============================================================================
# Dataset 07: Lin et al., EBioMedicine (2023)
# Original script: single-cell/Lin_EBioMedicine_2023/forDB.R
# =============================================================================


Lin_endo <- readRDS(file.path(project_dir, "Dataset07_input.rds"))

anno <- c("B_cell"="B Cell", "Cancer"="Epithelial Cell", "DC"="Myeloid Cell",
          "Epith"="Epithelial Cell",
          "Fibroblast"="Mesenchymal Cell",
          "Lymphatic endo_cell"="Lymphatic Endothelial Cell",
          "Macrophage"="Myeloid Cell", "Mast_cell"="Mast Cell",
          "Plasma_B_cell"="Plasma Cell",
          "Smooth_muscle_cell"="Mesenchymal Cell",
          "T_cell"="T Cell",
          "Vascular endo_cell"="Vascular Endothelial Cell")
Lin_endo[['cell_type']] = unname(anno[Lin_endo@meta.data$custom_anno])

Lin_endo <- subset(Lin_endo, cell_type == "Lymphatic Endothelial Cell" |
                 cell_type == "Vascular Endothelial Cell")
DefaultAssay(Lin_endo) <- "RNA"
Lin_endo@assays$integrated <- NULL
Lin_endo <- FindVariableFeatures(Lin_endo, selection.method = "vst",
                             nfeatures = 2000)
Lin_endo <- ScaleData(Lin_endo)
Lin_endo <- RunPCA(Lin_endo, features = VariableFeatures(object = Lin_endo),
               verbose = F)
Lin_endo <- RunHarmony(Lin_endo, group.by.vars = "orig.ident", verbose = F)
Lin_endo <- FindNeighbors(Lin_endo, reduction = "harmony", dims = 1:30)
Lin_endo <- RunUMAP(Lin_endo, reduction = "harmony", dims = 1:30)
Lin_endo <- RunTSNE(Lin_endo, reduction = "harmony", dims = 1:30)

Lin_endo <- FindClusters(Lin_endo, resolution = 0.15)

anno <- c("0"="vEC-c1", "1"="lEC", "2"="vEC-c2", "3"="vEC-c2", "4"="vEC-c3",
          "5"="Other", "6"="Other", "7"="vEC-c3")
Lin_endo[['cell_subtype']] = unname(anno[Lin_endo@meta.data$seurat_clusters])

Lin_endo <- subset(Lin_endo, cell_subtype != "Other")
Idents(Lin_endo) <- "cell_subtype"
Lin_endo$cell_subtype <- factor(Lin_endo$cell_subtype,
                         levels = c("vEC-c1", "vEC-c2", "vEC-c3",
                                    "lEC"))


# =============================================================================
# Dataset 08: Qu et al., Cancer Communications (2023), GSE197461
# Original script: single-cell/Qu_Cancer.Commun_2023/GSE197461.R
# =============================================================================

# Sample integration, clustering, and annotation
mat <- list.dirs(file.path(project_dir, "Dataset08"), recursive = FALSE)
sceList <- lapply(mat,function(x){
  sce <- CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
for (i in seq_along(sceList)) {
  sce<-sceList[[i]]
  sceList[[i]][["orig.ident"]] <- paste0("Sample", i)
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
}

sce <- merge(sceList[[1]],
             y = sceList[-1],
             add.cell.ids = paste0("Sample", seq_along(sceList)))
sce <- subset(sce, subset = nFeature_RNA >= 200 & percent.mt <= 20)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst")
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- RunHarmony(sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
sce <- FindNeighbors(sce, dims = 1:10,reduction="harmony")
sce <- FindClusters(sce, resolution = 0.8)
sce <- RunUMAP(sce,reduction="harmony", dims = 1:10)
sce <- RunTSNE(sce,reduction="harmony", dims = 1:10)


new.cluster.ids <- c('T cells','T cells','T cells','Fibroblasts','Epithelial cells','Fibroblasts','B/plasma cells','Epithelial cells',
                     'Myeloid cells','T cells','T cells','Epithelial cells','Smooth muscle cells','B/plasma cells','Myeloid cells','Epithelial cells',
                     'T cells','Mast cells','Endothelial cells','Epithelial cells','Myeloid cells','Epithelial cells')
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce@meta.data[["Cell_type"]]<-as.character(sce@active.ident)
sce@meta.data[["Cell_type"]]<-factor(sce@meta.data[["Cell_type"]],
                                     levels=c('Epithelial cells','Fibroblasts','Smooth muscle cells','Endothelial cells',
                                              'T cells','Myeloid cells','B/plasma cells','Mast cells'))
Idents(sce)<-sce@meta.data[["Cell_type"]]


# Myeloid-cell subclustering
myeloid<-subset(sce,Cell_type =='Myeloid cells')
myeloid <- RunPCA(myeloid, features = VariableFeatures(object = myeloid))
myeloid <- RunHarmony(myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
myeloid <- FindNeighbors(myeloid, dims = 1:30,reduction="harmony")
myeloid <- FindClusters(myeloid, resolution = 1)
myeloid <- RunUMAP(myeloid,reduction="harmony", dims = 1:30)
myeloid <- RunTSNE(myeloid,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('Macrophages','Monocytes','Dendritic cells','Macrophages','Macrophages','Monocytes','Macrophages',
                     'Proliferation','Macrophages','Neutrophils','Dendritic cells','Dendritic cells','Dendritic cells')
names(new.cluster.ids) <- levels(myeloid)
myeloid <- RenameIdents(myeloid, new.cluster.ids)
myeloid@meta.data[["Cell_subtype"]]<-as.character(myeloid@active.ident)
Idents(myeloid)<-myeloid@meta.data[["Cell_subtype"]]


# T-cell subclustering
t_nk<-subset(sce,Cell_type =='T cells')
t_nk <- RunPCA(t_nk, features = VariableFeatures(object = t_nk))
t_nk <- RunHarmony(t_nk,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
t_nk <- FindNeighbors(t_nk, dims = 1:20,reduction="harmony")
t_nk <- FindClusters(t_nk, resolution = 0.6)
t_nk <- RunUMAP(t_nk,reduction="harmony", dims = 1:20)


new.cluster.ids <- c('CD4','CD8','CD4','CD8','CD8','CD4','NK','CD8','CD4','CD8','NK','CD8','CD4')
names(new.cluster.ids) <- levels(t_nk)
t_nk <- RenameIdents(t_nk, new.cluster.ids)
t_nk@meta.data[["Cell_subtype"]]<-as.character(t_nk@active.ident)
Idents(t_nk)<-t_nk@meta.data[["Cell_subtype"]]

#CD4
cd4_t<-subset(t_nk,Cell_subtype =='CD4')
cd4_t <- RunPCA(cd4_t, features = VariableFeatures(object = cd4_t))
cd4_t <- RunHarmony(cd4_t,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
cd4_t <- FindNeighbors(cd4_t, dims = 1:20,reduction="harmony")
cd4_t <- FindClusters(cd4_t, resolution = 0.7)
cd4_t <- RunUMAP(cd4_t,reduction="harmony", dims = 1:20)
cd4_t <- RunTSNE(cd4_t,reduction="harmony", dims = 1:20)


new.cluster.ids <- c('CD4_naive','CD4_Treg','CD4_Treg','CD4_Tem','CD4_Tem','CD4_Th1-like','CD4_naive','CD4_naive','CD4_Treg','CD4_naive','CD4_naive','CD4_Tem')
names(new.cluster.ids) <- levels(cd4_t)
cd4_t <- RenameIdents(cd4_t, new.cluster.ids)
cd4_t@meta.data[["Cell_subtype"]]<-as.character(cd4_t@active.ident)
cd4_t@meta.data[["Cell_subtype"]]<-factor(cd4_t@meta.data[["Cell_subtype"]],
                                            levels=c('Tem','Th1-like','Th17','Naive','Treg'))
Idents(cd4_t)<-cd4_t@meta.data[["Cell_subtype"]]


#CD8
cd8_t<-subset(t_nk,Cell_subtype =='CD8')
cd8_t <- RunPCA(cd8_t, features = VariableFeatures(object = cd8_t))
cd8_t <- RunHarmony(cd8_t,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
cd8_t <- FindNeighbors(cd8_t, dims = 1:20,reduction="harmony")
cd8_t <- FindClusters(cd8_t, resolution = 0.7)
cd8_t <- RunUMAP(cd8_t,reduction="harmony", dims = 1:20)
cd8_t <- RunTSNE(cd8_t,reduction="harmony", dims = 1:20)


new.cluster.ids <- c('CD8_Trm','CD8_Tem','CD8_Tex','CD8_Tex','CD8_Cycling','CD8_Trm','CD8_Trm','CD8_Trm','CD8_Trm','CD8_MAIT','CD8_Cycling','CD8_Tex')
names(new.cluster.ids) <- levels(cd8_t)
cd8_t <- RenameIdents(cd8_t, new.cluster.ids)
cd8_t@meta.data[["Cell_subtype"]]<-as.character(cd8_t@active.ident)
Idents(cd8_t)<-cd8_t@meta.data[["Cell_subtype"]]


# Epithelial-cell subclustering
epithelial<-subset(sce,Cell_type =='Epithelial cells')
epithelial <- RunPCA(epithelial, features = VariableFeatures(object = epithelial))
epithelial <- RunHarmony(epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
epithelial <- FindNeighbors(epithelial, dims = 1:30,reduction="harmony")
epithelial <- FindClusters(epithelial, resolution = 1)
epithelial <- RunUMAP(epithelial,reduction="harmony", dims = 1:30)
epithelial <- RunTSNE(epithelial,reduction="harmony", dims = 1:30)


# =============================================================================
# Dataset 09: Zhang et al., EBioMedicine (2023)
# Original script: single-cell/Zhang_EBioMedicine_2023/forDB.R
# =============================================================================

#########Huang Ebiomedicine
Huang_scRNA <- readRDS(file.path(project_dir, "Dataset09_all_cells.rds"))
Huang_scRNA <- RunTSNE(Huang_scRNA,reduction="pca", dims = 1:20)


Huang_epi <- readRDS(file.path(project_dir, "Dataset09_epithelial.rds"))
Huang_epi <- RunTSNE(Huang_epi,reduction="pca", dims = 1:20)

Huang_T <- readRDS(file.path(project_dir, "Dataset09_t_cells.rds"))

Huang_fib <- subset(Huang_scRNA,cell_type=='Fibroblasts')
Huang_fib <- NormalizeData(Huang_fib, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_fib <- FindVariableFeatures(Huang_fib, selection.method = "vst", nfeatures = 3000)
Huang_fib <- ScaleData(Huang_fib, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_fib))
Huang_fib <- RunPCA(Huang_fib, features = VariableFeatures(object = Huang_fib))
Huang_fib <- RunHarmony(Huang_fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Huang_fib <- FindNeighbors(Huang_fib, dims = 1:20,reduction="harmony")
Huang_fib <- RunUMAP(Huang_fib,reduction="harmony", dims = 1:20)
Huang_fib <- RunTSNE(Huang_fib,reduction="harmony", dims = 1:20)


Huang_endo <- subset(Huang_scRNA,cell_type=='Endothelial cell')
Huang_endo <- NormalizeData(Huang_endo, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_endo <- FindVariableFeatures(Huang_endo, selection.method = "vst", nfeatures = 3000)
Huang_endo <- ScaleData(Huang_endo, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_endo))
Huang_endo <- RunPCA(Huang_endo, features = VariableFeatures(object = Huang_endo))
Huang_endo <- RunHarmony(Huang_endo,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Huang_endo <- FindNeighbors(Huang_endo, dims = 1:20,reduction="harmony")
Huang_endo <- RunUMAP(Huang_endo,reduction="harmony", dims = 1:20)
Huang_endo <- RunTSNE(Huang_endo,reduction="harmony", dims = 1:20)


Huang_macro <- subset(Huang_scRNA,cell_type=='Macrophages')
Huang_macro <- NormalizeData(Huang_macro, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_macro <- FindVariableFeatures(Huang_macro, selection.method = "vst", nfeatures = 3000)
Huang_macro <- ScaleData(Huang_macro, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_macro))
Huang_macro <- RunPCA(Huang_macro, features = VariableFeatures(object = Huang_macro))
Huang_macro <- RunHarmony(Huang_macro,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Huang_macro <- FindNeighbors(Huang_macro, dims = 1:20,reduction="harmony")
Huang_macro <- RunUMAP(Huang_macro,reduction="harmony", dims = 1:20)
Huang_macro <- RunTSNE(Huang_macro,reduction="harmony", dims = 1:20)


Huang_B <- subset(Huang_scRNA,cell_type=='B cells')
Huang_B <- NormalizeData(Huang_B, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_B <- FindVariableFeatures(Huang_B, selection.method = "vst", nfeatures = 3000)
Huang_B <- ScaleData(Huang_B, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_B))
Huang_B <- RunPCA(Huang_B, features = VariableFeatures(object = Huang_B))
Huang_B <- RunHarmony(Huang_B,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Huang_B <- FindNeighbors(Huang_B, dims = 1:20,reduction="harmony")
Huang_B <- RunUMAP(Huang_B,reduction="harmony", dims = 1:20)
Huang_B <- RunTSNE(Huang_B,reduction="harmony", dims = 1:20)


# =============================================================================
# Dataset 10: Sandoval et al., Cancer Research (2026), GSE297041
# Original script: new-data/single-cell/Sandoval_Cancer.Res_2026(GSE297041)/GSE297041.R
# =============================================================================

seu <- readRDS(file.path(project_dir, "Dataset10_input.rds"))
seu <- RunTSNE(seu,dims = 1:30)
new.cluster.ids <- c('Immune','Immune','Immune','Immune','Tumor','Immune','Tumor','Tumor','Tumor','Tumor','Immune','Tumor','Stroma',
                     'Immune','Tumor','Stroma','Immune','Stroma','Tumor','Tumor','Tumor','Tumor','Immune','Immune','Immune')
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu@meta.data[["cell_type"]]<-seu@active.ident
seu@meta.data[["cell_type"]]<-factor(seu@meta.data[["cell_type"]],levels=unique(new.cluster.ids))

# Immune-cell subclustering
Immune<-subset(seu,subset=cell_type =='Immune')
Immune <- SCTransform(Immune, vars.to.regress = c('nCount_RNA', 'percent.mito'),variable.features.n = 2000)
Immune <- RunPCA(Immune, features = VariableFeatures(object = Immune))
# Immune-cell subclustering <- RunHarmony(Immune,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Immune <- FindNeighbors(Immune, dims = 1:30,reduction="pca")
Immune <- FindClusters(Immune, resolution = 0.5)
Immune <- RunUMAP(Immune,reduction="pca", dims = 1:30)
Immune <- RunTSNE(Immune,dims = 1:30)


new.cluster.ids <- c('CD8+ Tcells','HCAR2/3+ Neutrophils','CD8+ Tcells','Monocytes','Treg','CXCR2+ Neutrophils','SPP1+ Macrophages',
                     'CD4+ Tcells','DCs','TREM2+ Macrophages','CCL3/4+ Neutrophils','Plasma cells','CD8+ Tcells','B cells',
                     'Unknown','Proliferating Tcells','Plasma cells','Mast','SPP1+ Macrophages','Plasma cells','pDC','Monocytes')
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
Immune@meta.data[["cell_subtype"]]<-Immune@active.ident
Immune@meta.data[["cell_subtype"]]<-factor(Immune@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


# =============================================================================
# Dataset 11: Yuan et al., Frontiers in Immunology (2026), GSE308792
# Original script: new-data/single-cell/Yuan_Front.Immunol_2026(GSE308792)/GSE308792.R
# =============================================================================

files <- list.dirs(file.path(project_dir, "Dataset11"), recursive = FALSE)
seuList <- lapply(files,function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})

for (i in seq_along(seuList)) {
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 200 & nFeature_RNA <= 7500 & nCount_RNA >= 500 & nCount_RNA <= 40000 & percent.mt <= 20)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  seuList[[i]]<-sce
}

for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")
seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge, selection.method = "vst", nfeatures = 2000)
seu_merge <- ScaleData(seu_merge, vars.to.regress = "percent.mt")
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:30,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)

new.cluster.ids <- c('T cell','T cell','Epithelial cell','Epithelial cell','Epithelial cell','Epithelial cell',
                     'Epithelial cell','Epithelial cell','Epithelial cell','Fibroblast','Macrophages','Plasma cell',
                     'Endothelial cell','B cell','Epithelial cell','pDC','Epithelial cell','Epithelial cell',
                     'Epithelial cell','Mast cell')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))


#T
Tcell<-subset(seu_merge,subset=cell_type == 'T cell')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:30,reduction="harmony")
Tcell <- FindClusters(Tcell, resolution = 0.3)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:30)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('CD8 Tex','CD4 Treg','CD4 Naive T','CD8 Naive T','CD4 Naive T','CD8 Tex','CD16+ NK',
                     'CD4 Tph','Memory T cells','CD56+ NK','Unknown')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids))[c(2:4,6:8,1,5,9)])


# =============================================================================
# Dataset 12: Wu et al., Communications Biology (2025), E-MTAB-15983
# Original script: new-data/single-cell/Wu_Commun.Biol_2025/E-MTAB-15983.R
# =============================================================================

files <- list.files(file.path(project_dir, "Dataset12"), full.names = TRUE)
seuList <- lapply(files,function(x){
  seu = CreateSeuratObject(counts = Read10X_h5(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
# TODO: Sample 1 is skipped here, but downstream code expects DoubletFinder fields for every sample.
for(i in 2:length(seuList)){
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 500 & nFeature_RNA <= 10000 & nCount_RNA >= 500 & nCount_RNA <= 20000 & percent.mt <= 20)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  seuList[[i]]<-sce
}

for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")
seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge, selection.method = "vst", nfeatures = 2000)
seu_merge <- ScaleData(seu_merge, vars.to.regress = "percent.mt")
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:30,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)

new.cluster.ids <- c('NK/T','Epithelial','NK/T','Fibroblast','Myeloid','NK/T','Endothelial','Epithelial','B','CAF','Epithelial',
                     'NK/T','Plasma','Mast cell','Myeloid','Epithelial','Myeloid','Epithelial','Fibroblast','Fibroblast','Myeloid',
                     'Epithelial','Epithelial','Endothelial','Myeloid')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))

Epithelial<-subset(Epithelial,subset=cell_type=='Epithelial')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:30,reduction="harmony")
Epithelial <- FindClusters(Epithelial, resolution = 0.2)
Epithelial <- RunUMAP(Epithelial,reduction="harmony", dims = 1:30)
Epithelial <- RunTSNE(Epithelial,reduction="harmony", dims = 1:30)
Epithelial$Cell_subtype<-factor(paste0('Epi',Epithelial@meta.data$seurat_clusters),levels = paste0('Epi',levels(Epithelial@meta.data$seurat_clusters)))
Idents(Epithelial)<-Epithelial$Cell_subtype

Immune<-subset(seu_merge,subset=cell_type %in% c('NK/T','Myeloid','B','Plasma'))
Immune <- NormalizeData(Immune, normalization.method = "LogNormalize", scale.factor = 10000)
Immune <- FindVariableFeatures(Immune, selection.method = "vst", nfeatures = 2000)
Immune <- ScaleData(Immune, vars.to.regress = "percent.mt")
Immune <- RunPCA(Immune, features = VariableFeatures(object = Immune))
Immune <- RunHarmony(Immune,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Immune <- FindNeighbors(Immune, dims = 1:30,reduction="harmony")
Immune <- FindClusters(Immune, resolution = 0.6)
Immune <- RunUMAP(Immune,reduction="harmony", dims = 1:30)
Immune <- RunTSNE(Immune,reduction="harmony", dims = 1:30)


new.cluster.ids <- c('CD8 T','CD4 Naive','CD8 T','CD4 Treg','CD163 Mac','CD163 Mac','B','CD4 Tph','NK','Proliferating NK/CD8 T','Plasma','Neu',
                     'CD8 T','Mono/Macro','Mono/Macro','Proliferating Mono/Mac','Mono/Macro','Mono/Macro','Mono/Macro','CD4 Treg','pDC')
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
Immune@meta.data[["cell_subtype"]]<-Immune@active.ident
Immune@meta.data[["cell_subtype"]]<-factor(Immune@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


# =============================================================================
# Dataset 13: Hyeon et al., Molecular Cancer (2025), GSE279998
# Original script: new-data/single-cell/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.R
# =============================================================================

sample_dirs <- list.dirs(file.path(project_dir, "Dataset13"), recursive = FALSE)
seuList <- lapply(sample_dirs,function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})

for (i in seq_along(seuList)) {
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 500 & nFeature_RNA <= 5000 & nCount_RNA >= 2000 & percent.mt <= 10)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  seuList[[i]]<-sce
}

for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")
seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge, selection.method = "vst", nfeatures = 2000)
seu_merge <- ScaleData(seu_merge, vars.to.regress = "percent.mt")
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:24,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:24)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:24)

new.cluster.ids <- c('Epithelial','CD8 T cell','Treg','CD4 T cell','NK cell','CD8 T cell','Plasma cell','CAF','CD8 T cell',
                     'Macrophage','B cell','CAF/SMC','Endothelial','Plasma cell','Neutrophil','Proliferating T cell',
                     'Mast','Epithelial','pDC','Plasma cell','Plasma cell')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))

# Epithelial-cell subclustering
Epithelial<-subset(seu_merge,subset=cell_type=='Epithelial')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:50,reduction="harmony")
Epithelial <- FindClusters(Epithelial, resolution = 0.3)
Epithelial <- RunUMAP(Epithelial,reduction="harmony", dims = 1:50)
Epithelial <- RunTSNE(Epithelial,reduction="harmony", dims = 1:50)
new.cluster.ids <- c('AKR1C2','GADD45B','CEACAM5','KRT1/3','Proliferating','GADD45B','Unknown','MUC5AC','PGK1','Unknown','Unknown','CXCL10')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


#CAF
Fib<-subset(seu_merge,subset=cell_type=='CAF')
Fib <- NormalizeData(Fib, normalization.method = "LogNormalize", scale.factor = 10000)
Fib <- FindVariableFeatures(Fib, selection.method = "vst", nfeatures = 2000)
Fib <- ScaleData(Fib, vars.to.regress = "percent.mt")
Fib <- RunPCA(Fib, features = VariableFeatures(object = Fib))
Fib <- RunHarmony(Fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Fib <- FindNeighbors(Fib, dims = 1:50,reduction="harmony")
Fib <- FindClusters(Fib, resolution = 0.7)
Fib <- RunUMAP(Fib,reduction="harmony", dims = 1:50)
Fib <- RunTSNE(Fib,reduction="harmony", dims = 1:50)
new.cluster.ids <- c('PDGFD','CCN5','CXCL1','PDGFD','PLN','PDGFD','CCL21','Unknown','DES','PDGFD','RGS5','CXCL10','Unknown','Unknown')
names(new.cluster.ids) <- levels(Fib)
Fib <- RenameIdents(Fib, new.cluster.ids)
Fib@meta.data[["cell_subtype"]]<-Fib@active.ident
Fib@meta.data[["cell_subtype"]]<-factor(Fib@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type %in% c('Macrophage','pDC','Mast'))
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:27,reduction="harmony")
Mye <- FindClusters(Mye, resolution = 0.2)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:27)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:27)
new.cluster.ids <- c('C1QB+ Mac','FCN1+ Mac','Mast','CD1C+ DC','C1QB+ Mac','pDC','LAMP3+ DC','CLEC9A+ Mac','Proliferating','C1QB+ Mac')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))

#T
Tcell<-subset(seu_merge,subset=cell_type %in% c('CD8 T cell','Treg','CD4 T cell','NK cell'))
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:21,reduction="harmony")
Tcell <- FindClusters(Tcell, resolution = 0.3)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:21)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:21)
new.cluster.ids <- c('IL7R+ T','GZMK+ T','LAG3+ T','Treg','FCER1G+ NK','IL7R+ T','FCGR3A+ NK','CD8+ T','IFIT3+ T')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


# =============================================================================
# Dataset 14: Liu et al., Journal of Experimental and Clinical Cancer Research (2023), SCP1950
# Original script: new-data/single-cell/Liu_J.Exp.Clin.Cancer.Res_2023/SCP1950.R
# =============================================================================

sample_dirs <- list.dirs(file.path(project_dir, "Dataset14"), recursive = FALSE)
seuList <- lapply(sample_dirs,function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})

for (i in seq_along(seuList)) {
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 500 & nFeature_RNA <= 4000 & nCount_RNA < 8000 & percent.mt <= 10)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  seuList[[i]]<-sce
}

for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")
seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge, selection.method = "vst", nfeatures = 2000)
seu_merge <- ScaleData(seu_merge, vars.to.regress = "percent.mt")
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:30,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)

new.cluster.ids <- c('Basal cells','Mesenchymal stem cells','Basal cells','Cancer stem cells','Cancer stem cells',
                     'CD8 T cells','Proliferating Epithelial','Macrophage','Epithelial cells','T cells','naive B cells',
                     'Treg','myofibroblast','Cancer stem cells','mature B cells','Epithelial cells','Endothelial cells',
                     'Pericytes','T cells','CD141+CLEC9A+ DC','γδT','Lymphatic endothelial cells')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))

# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type %in% c('Macrophage','CD141+CLEC9A+ DC'))
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:30,reduction="harmony")
Mye <- FindClusters(Mye, resolution = 0.6)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:30)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('C0-Res','','','C3-DC','C2-M2','C1-TAM','','C4-M1','','','','C3-DC','','','')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))

# =============================================================================
# Dataset 15: Dai et al., Cell Reports Medicine (2024), GSE236738
# Original script: new-data/single-cell/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.R
# =============================================================================

sample_dirs <- list.dirs(file.path(project_dir, "Dataset15"), recursive = FALSE)
seuList <- lapply(sample_dirs,function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})

for (i in seq_along(seuList)) {
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA <= 6000 & nCount_RNA > 2000 & percent.mt <= 20)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  seuList[[i]]<-sce
}

for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")
seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge, selection.method = "vst", nfeatures = 2000)
seu_merge <- ScaleData(seu_merge, vars.to.regress = "percent.mt")
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:50,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.4)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:50)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:50)

new.cluster.ids <- c('T cell','T cell','Monocyte','T cell','T cell','Plasma cell','Fibroblast','Tumor cell','Tumor cell',
                     'Tumor cell','T cell','Tumor cell','Neutrophil','Mast cell','B cell','Fibroblast','T cell','Tumor cell',
                     'Blood cell','Endothelial cell','pDC')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))

# Epithelial-cell subclustering
Epithelial<-subset(seu_merge,subset=cell_type=='Tumor cell')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:50,reduction="harmony")
Epithelial <- FindClusters(Epithelial, resolution = 0.3)
Epithelial <- RunUMAP(Epithelial,reduction="harmony", dims = 1:50)
Epithelial <- RunTSNE(Epithelial,reduction="harmony", dims = 1:50)
Epithelial@meta.data[["cell_subtype"]]<-paste0('Epi',Epithelial$seurat_clusters)
# TODO: These levels come from broad cell-type labels and may convert the Epi labels to NA.
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
Idents(Epithelial)<-Epithelial$cell_subtype


# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type %in% c('Monocyte','pDC','Neutrophil'))
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:50,reduction="harmony")
Mye <- FindClusters(Mye, resolution = 0.4)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:50)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:50)

new.cluster.ids <- c('Macrophage','Monocyte','Neutrophil','Macrophage','Monocyte','Macrophage','Unknown','Unknown','pDC')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))

#T
Tcell<-subset(seu_merge,subset=cell_type == 'T cell')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:50,reduction="harmony")
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:50)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:50)

new.cluster.ids <- c('CD4 Tna','CD8 T','CD4 Treg','CD56+ NK','CD8 T','CD16+ NK','CD8 T','CD4 Th1-like','Proliferating T','CD8 Tex',
                     'Unknown','Unknown','CD56+ NK','Unknown')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


# =============================================================================
# Dataset 16: Peng et al., eLife (2025), SRP567748
# Original script: new-data/single-cell/Peng_eLife_2025/SRP567748.R
# =============================================================================

sample_dirs <- list.dirs(file.path(project_dir, "Dataset16"), recursive = FALSE)
seuList <- lapply(paste0(sample_dirs,'/filtered_feature_bc_matrix'),function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})

for (i in seq_along(seuList)) {
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 200 & percent.mt <= 25)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = FALSE)
  seuList[[i]]<-sce
}

for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")

seu_merge <- readRDS(file.path(project_dir, "Dataset16_merged.rds"))
seu_merge <- subset(seu_merge, subset = nFeature_RNA <= 8000 & nCount_RNA<=40000)

seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge)
seu_merge <- ScaleData(seu_merge)
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:30,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.4)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)

new.cluster.ids <- c('T cell','Neutrophil','T cell','Epithelial cell','Epithelial cell','Macrophage','Epithelial cell',
                     'Epithelial cell','Plasma cell','Fibroblast','Epithelial cell','Epithelial cell','T cell',
                     'Epithelial cell','Endothelial cell','B cell','Fibroblast','Mast cell','Epithelial cell','T cell')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))

# Epithelial-cell subclustering
Epithelial<-subset(seu_merge,subset=cell_type=='Epithelial cell')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial)
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:30,reduction="pca")
Epithelial <- FindClusters(Epithelial, resolution = 0.3)
Epithelial <- RunUMAP(Epithelial,reduction="pca", dims = 1:30)
Epithelial <- RunTSNE(Epithelial,reduction="pca", dims = 1:30)
new.cluster.ids <- c('Epi_02_IGLC2','Epi_01_NTS','Epi_03_TMPRSS11E','Epi_05_CCL5','Epi_03_TMPRSS11E','Epi_05_CCL5',
                     'Epi_02_IGLC2','Epi_10_CYSTM1','Epi_01_NTS','Epi_04_TFF2','Epi_06_TMC5','Epi_06_TMC5','Epi_09_SST',
                     'Epi_02_IGLC2','Epi_07_CAPS','Epi_12_RRAD','Epi_06_TMC5','Epi_02_IGLC2','Epi_08_SCGB3A1','Epi_07_CAPS')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids)))


# Neutrophil subclustering
Neu<-subset(seu_merge,subset=cell_type =='Neutrophil')
Neu@assays[["RNA"]]@meta.data<-data.frame()
Neu <- NormalizeData(Neu, normalization.method = "LogNormalize", scale.factor = 10000)
Neu <- FindVariableFeatures(Neu, selection.method = "vst", nfeatures = 2000)
Neu <- ScaleData(Neu)
Neu <- RunPCA(Neu, features = VariableFeatures(object = Neu))
Neu <- RunHarmony(Neu,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Neu <- FindNeighbors(Neu, dims = 1:30,reduction="harmony")
Neu <- FindClusters(Neu, resolution = 0.5)
Neu <- RunUMAP(Neu,reduction="harmony", dims = 1:30)
Neu <- RunTSNE(Neu,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('Neu_01','Neu_01','Neu_01','Neu_01','Neu_01','Neu_01','Neu_01','Neu_02','Neu_03')
names(new.cluster.ids) <- levels(Neu)
Neu <- RenameIdents(Neu, new.cluster.ids)
Neu@meta.data[["cell_subtype"]]<-Neu@active.ident
Neu@meta.data[["cell_subtype"]]<-factor(Neu@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))

#T
Tcell<-subset(seu_merge,subset=cell_type == 'T cell')
Tcell@assays[["RNA"]]@meta.data<-data.frame()
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell)
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:30,reduction="harmony")
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:30)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('Exhausted T','Exhausted T','Treg','Naive T','Cytotoxic T','Tph','Proliferating T','Cytotoxic T',
                     'Exhausted T','Cytotoxic T','Exhausted T','CD16+ NK','Activated T','CD56+ NK')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))


# B-cell and plasma-cell subclustering
Plasma<-subset(seu_merge,subset=cell_type %in% c('B cell','Plasma cell'))
Plasma@assays[["RNA"]]@meta.data<-data.frame()
Plasma <- NormalizeData(Plasma, normalization.method = "LogNormalize", scale.factor = 10000)
Plasma <- FindVariableFeatures(Plasma, selection.method = "vst", nfeatures = 2000)
Plasma <- ScaleData(Plasma)
Plasma <- RunPCA(Plasma, features = VariableFeatures(object = Plasma))
Plasma <- RunHarmony(Plasma,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Plasma <- FindNeighbors(Plasma, dims = 1:30,reduction="harmony")
Plasma <- FindClusters(Plasma, resolution = 0.6)
Plasma <- RunUMAP(Plasma,reduction="harmony", dims = 1:30)
Plasma <- RunTSNE(Plasma,reduction="harmony", dims = 1:30)
new.cluster.ids <- c('Plasma/B_02','Plasma/B_01','Plasma/B_01','Plasma/B_01','Plasma/B_01','Plasma/B_04','Plasma/B_02',
                     'Plasma/B_03','Plasma/B_05','Plasma/B_06','Plasma/B_01','Plasma/B_05')
names(new.cluster.ids) <- levels(Plasma)
Plasma <- RenameIdents(Plasma, new.cluster.ids)
Plasma@meta.data[["cell_subtype"]]<-Plasma@active.ident
Plasma@meta.data[["cell_subtype"]]<-factor(Plasma@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids)))


# =============================================================================
# Dataset 17: Cao et al., Journal for ImmunoTherapy of Cancer (2025), treatment cohort
# Original script: new-data/single-cell/Cao_J.Immunother.Cancer_2025/Treatment.R
# =============================================================================

filename <- list.dirs(file.path(project_dir, "Dataset17"), recursive = FALSE)
seuList <- lapply(filename,function(x){
  sce <- CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
for (i in seq_along(seuList)) {
  sce<-seuList[[i]]
  sce[["orig.ident"]] <- paste0("Sample", i)
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & nCount_RNA <= 45000 & percent.mt <= 25)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  # TODO: The parameter sweep uses PCs 1:15, while doubletFinder below uses PCs 1:20.
sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion
  homotypic.prop <- modelHomotypic(sce$seurat_clusters)
  # Estimate the expected doublet count
  #DoubletRate = ncol(sce)*8*1e-6
  DoubletRate = 0.08
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, sct = FALSE)
  seuList[[i]]<-sce
}
for (i in seq_along(seuList)) {
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                    y = seuList[-1],
                    add.cell.ids = paste0("Sample", seq_along(seuList)))
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")
seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge, selection.method = "vst", nfeatures = 2000)
seu_merge <- ScaleData(seu_merge, vars.to.regress = "percent.mt")
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:15,reduction="harmony")
seu_merge <- FindClusters(seu_merge, resolution = 0.1)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:15)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:15)

new.cluster.ids <- c('Epithelial','T/NK','Epithelial','Fibroblast','Myeloid','B/Plasma','B/Plasma','Endothelial','Fibroblast','Mast',
                     'T/NK','Epithelial','Epithelial','Epithelial')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))


# Epithelial-cell subclustering
Epithelial <- subset(seu_merge,subset=cell_type=='Epithelial')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
# NOTE: Harmony is calculated above, but the original downstream epithelial workflow uses PCA.
Epithelial <- FindNeighbors(Epithelial, dims = 1:15,reduction="pca")
Epithelial <- FindClusters(Epithelial, resolution = 0.2)
Epithelial <- RunUMAP(Epithelial,reduction="pca", dims = 1:15)
Epithelial <- RunTSNE(Epithelial,reduction="pca", dims = 1:15)
new.cluster.ids <- c('EP0_SPP1','EP8_SERPINA1','EP0_SPP1','EP3_CEACAM5','Unknown','EP1_EPCAM','EP0_SPP1','EP6_NEAT1',
                     'EP3_CEACAM5','EP2_MYC','EP3_CEACAM5','EP4_MUC5B','EP0_SPP1','EP0_SPP1','EP0_SPP1','EP6_NEAT1')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids)))


# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type =='Myeloid')
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:15,reduction="harmony")
Mye <- FindClusters(Mye, resolution = 0.4)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:15)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:15)
new.cluster.ids <- c('Macro_CCL18','Macro_IL1B','Macro_IL1B','Macro_SPP1','Unknown','Macro_C1QC','Unknown','Macro_SPP1','cDC_CD1C',
                     'pDC_GZMB','Macro_MKI67','Macro_IL1B','Macro_SPP1','mDC_LAMP3','Macro_SPP1')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))

#T
Tcell<-subset(seu_merge,subset=cell_type =='T/NK')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:15,reduction="harmony")
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:15)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:15)
new.cluster.ids <- c('CD4_Tcm','CD4_Treg','CD8_CTL','CD4_Tn','CD8_Tex_PDCD1','CD8_Trm_ZNF683','CD16+_NK','CD8_Tn',
                     'CD56+_NK','Proliferating_T','CD8_Tn','Proliferating_T',
                     'CD16+_NK','Proliferating_T','Unknown','CD4_Treg')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))

