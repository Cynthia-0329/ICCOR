# Purpose: Process and annotate 18 cervical cancer single-cell RNA-seq datasets.
# Input: Dataset-specific count matrices, metadata, and intermediate Seurat objects listed in each section.
# Output: Dataset-specific annotated Seurat objects and shared subtype marker collections.
#
# Each dataset section preserves its original preprocessing, QC thresholds,
# integration method, clustering parameters, manual annotations, and output calls.
# Run only the required section unless all source data and intermediate objects
# are available. Several sections update the same subtype marker RData file.

project_dir <- Sys.getenv("ICCOR_PROJECT_DIR", unset = "/Volumes/Cynthia/项目/宫颈癌数据库") # TODO: update for your system
local_project_dir <- Sys.getenv("ICCOR_LOCAL_PROJECT_DIR", unset = "/Users/cynthia/项目") # TODO: update for your system
user_project_dir <- Sys.getenv("ICCOR_USER_PROJECT_DIR", unset = path.expand("~/项目")) # TODO: update for your system

# Dataset index
# 01  Cao et al., EMBO Journal (2023)
# 02  Guo et al., Clinical and Translational Medicine (2023)
# 03  Li et al., Cancer Cell International (2025)
# 04  Li et al., Communications Biology (2022), E-MTAB-11948
# 05  Li et al., Frontiers in Immunology (2022), E-MTAB-12305
# 06  Li et al., Journal of Medical Virology (2023), S-BSST1035
# 07  Li et al., Molecular Therapy Nucleic Acids (2021), GSE168652
# 08  Lin et al., EBioMedicine (2023)
# 09  Qu et al., Cancer Communications (2023), GSE197461
# 10  Zhang et al., EBioMedicine (2023)
# 11  Sandoval et al., Cancer Research (2026), GSE297041
# 12  Yuan et al., Frontiers in Immunology (2026), GSE308792
# 13  Wu et al., Communications Biology (2025), E-MTAB-15983
# 14  Hyeon et al., Molecular Cancer (2025), GSE279998
# 15  Liu et al., Journal of Experimental and Clinical Cancer Research (2023), SCP1950
# 16  Dai et al., Cell Reports Medicine (2024), GSE236738
# 17  Peng et al., eLife (2025), SRP567748
# 18  Cao et al., Journal for ImmunoTherapy of Cancer (2025), treatment cohort
# =============================================================================
# Dataset 01: Cao et al., EMBO Journal (2023)
# Original script: single-cell/Cao_EMBO.J_2023/Lifang_for_DB.R
# =============================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(DoubletFinder)
library(clustree)
library(glmGamPoi)
library(ggsci)
# Sample integration, clustering, and annotation
mat <- c(file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P1/CA_N_5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P1/CA_5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P2/R2020003P_5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P2/R2020002CA_5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P3/R2020012P_5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P3/R2020012ca_5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P4/P4_Normal/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P4/P4_Tumor/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P5/P5_Tumor/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P6/P6_Normal/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "宫颈癌数据库/单细胞/matrix/P6/P6_Tumor/filtered_feature_bc_matrix"))
sceList <- lapply(mat,function(x){ 
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samples<-c('P1_N','P1_CA','P2_N','P2_CA','P3_N','P3_CA','P4_N','P4_CA','P5_CA','P6_N','P6_CA')
type<-c('Normal','CA','Normal','CA','Normal','CA','Normal','CA','CA','Normal','CA')
orig <- NULL
tissue <- NULL
CC.sce.doubletFinder.list<-list()
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  tissue = rep(type[i],ncol(sceList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["tissue"]] = tissue
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA > 1000 & nCount_RNA >500 & percent.mt < 15)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  #sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 0.5)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
  CC.sce.doubletFinder.list[[i]]<-sce
  print(paste0(samples[i]," Finished!"))
}

CC.sce<-CC.sce.doubletFinder.list
table(CC.sce.doubletFinder.list[[1]]@meta.data[["DF.classifications_0.25_11_465"]])
table(CC.sce.doubletFinder.list[[2]]@meta.data[["DF.classifications_0.25_2_1149"]])
table(CC.sce.doubletFinder.list[[3]]@meta.data[["DF.classifications_0.25_10_642"]])
table(CC.sce.doubletFinder.list[[4]]@meta.data[["DF.classifications_0.25_31_913"]])
table(CC.sce.doubletFinder.list[[5]]@meta.data[["DF.classifications_0.25_2_1149"]])
table(CC.sce.doubletFinder.list[[6]]@meta.data[["DF.classifications_0.25_28_992"]])
table(CC.sce.doubletFinder.list[[7]]@meta.data[["DF.classifications_0.25_30_2397"]])
table(CC.sce.doubletFinder.list[[8]]@meta.data[["DF.classifications_0.25_6_392"]])

for(i in 1:length(CC.sce)){
  CC.sce[[i]][["pANN"]]<-CC.sce[[i]]@meta.data %>% select(contains('pANN'))
  CC.sce[[i]][["DF.classifications"]]<-CC.sce[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}


CC.sce <- merge(CC.sce[[1]],
                y = CC.sce[-1],
                add.cell.ids = samples)
c <- grep("pANN_0.25",colnames(CC.sce@meta.data))
CC.sce@meta.data <- CC.sce@meta.data[,-c]
c <- grep("DF.classifications_0.25",colnames(CC.sce@meta.data))
CC.sce@meta.data <- CC.sce@meta.data[,-c]
DimPlot(CC.sce, group.by = 'DF.classifications', label = TRUE,raster=FALSE)


CC.sce <- subset(CC.sce, subset = DF.classifications== "Singlet")
saveRDS(CC.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.sce.rds"))
Li.paired.sce<-readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.sce.rds"))
#Li.paired.sce <- SCTransform(Li.paired.sce, vars.to.regress = "percent.mt", verbose = FALSE)
VlnPlot(Li.paired.sce,features = c('nFeature_RNA','nCount_RNA','percent.mt'),group.by = 'project.name')
Li.paired.sce <- NormalizeData(Li.paired.sce, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.sce <- FindVariableFeatures(Li.paired.sce, selection.method = "vst", nfeatures = 3000)
Li.paired.sce <- ScaleData(Li.paired.sce, vars.to.regress = "percent.mt")
Li.paired.sce <- RunPCA(Li.paired.sce, features = VariableFeatures(object = Li.paired.sce))
#Li.paired.sce <- RunHarmony(Li.paired.sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.sce <- IntegrateLayers(
  object = Li.paired.sce, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.sce)
)
ElbowPlot(Li.paired.sce, ndims = 50, reduction = "pca")
DimPlot(Li.paired.sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.sce, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
Li.paired.sce <- FindNeighbors(Li.paired.sce, dims = 1:30,reduction="cca")
Li.paired.sce <- FindClusters(object = Li.paired.sce,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.sce@meta.data, prefix = "RNA_snn_res.")
Li.paired.sce <- FindClusters(Li.paired.sce, resolution = 0.5)
Li.paired.sce <- RunUMAP(Li.paired.sce,reduction="cca", dims = 1:30)
Li.paired.sce <- RunTSNE(Li.paired.sce,reduction="cca", dims = 1:30)
DimPlot(Li.paired.sce, reduction = "umap", label = TRUE,raster=FALSE)
color<-colorRampPalette(brewer.pal(12, 'Paired'))(11)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2', 'KRT8','KRT18','KRT19','TSPAN8','MUC5B',   #epithelial cells
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CDH5', 'EMCN','PECAM1', 'PCDH17','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','MZB1','CD38',                      #plasma cells
          'MS4A1','CD79B','IGKC','CD79A',                     #B cells
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CPA3','TPSAB1','KIT')                     #mast cells
DotPlot(Li.paired.sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))


# Manual cluster annotation
new.cluster.ids <- c('T/NK','Fibroblasts','Epithelial','Fibroblasts','T/NK','Fibroblasts','Endothelial',
                     'Epithelial','Myeloid','T/NK',"B","T/NK",'T/NK','Plasma','Epithelial',
                     'Mast','Myeloid','Fibroblasts','T/NK','Epithelial','Fibroblasts')
names(new.cluster.ids) <- levels(Li.paired.sce)
Li.paired.sce <- RenameIdents(Li.paired.sce, new.cluster.ids)
Li.paired.sce@meta.data[["cell_type"]]<-Li.paired.sce@active.ident
Li.paired.sce@meta.data[["cell_type"]]<-factor(Li.paired.sce@meta.data[["cell_type"]],levels=c('Epithelial','Fibroblasts','T/NK',
                                                                                 'Endothelial','Myeloid','B','Plasma','Mast'))
Idents(Li.paired.sce)<-Li.paired.sce@meta.data[["cell_type"]]
DimPlot(Li.paired.sce, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.sce.anno.rds"))



Li.paired.epi<-subset(Li.paired.sce,cell_type=='Epithelial')
Li.paired.epi <- NormalizeData(Li.paired.epi, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.epi <- FindVariableFeatures(Li.paired.epi, selection.method = "vst", nfeatures = 3000)
Li.paired.epi <- ScaleData(Li.paired.epi, vars.to.regress = "percent.mt")
Li.paired.epi <- RunPCA(Li.paired.epi, features = VariableFeatures(object = Li.paired.epi))
Li.paired.epi <- RunHarmony(Li.paired.epi,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.epi <- IntegrateLayers(
  object = Li.paired.epi, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.epi)
)
ElbowPlot(Li.paired.epi, ndims = 50, reduction = "harmony")
DimPlot(Li.paired.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.epi, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
Li.paired.epi <- FindNeighbors(Li.paired.epi, dims = 1:30,reduction="harmony")
Li.paired.epi <- FindClusters(object = Li.paired.epi,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.epi@meta.data, prefix = "RNA_snn_res.")
Li.paired.epi <- FindClusters(Li.paired.epi, resolution = 0.6)
Li.paired.epi <- RunUMAP(Li.paired.epi,reduction="harmony", dims = 1:30)
Li.paired.epi <- RunTSNE(Li.paired.epi,reduction="harmony", dims = 1:30)
DimPlot(Li.paired.epi, reduction = "umap", label = TRUE,raster=FALSE)

Li.paired.epi<-JoinLayers(Li.paired.epi)
marker<-c('MUC5B','KRT6A','EPCAM','POSTN','KRT14','MKI67')
marker<-c('CDKN2A',
          'KRT6A','KRT5',
          'KRT7','KRT17',
          'S100A8','S100A9',
          'COL17A1','POSTN','SNAI2',
          'MKI67','UBE2C','CDC20',
          'MUC6','MUC5B',
          'HLA-A','HLA-B','HLA-C')
VlnPlot(Li.paired.epi,features = c('MUC5B','KRT6A','EPCAM','POSTN','KRT14','MKI67'))
DotPlot(Li.paired.epi, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('EP0_MUC5B','EP1_KRT6A','EP1_KRT6A','EP4_EPCAM','EP1_KRT6A','EP1_KRT6A',
                     'EP2_POSTN','EP3_MKI67','EP1_KRT6A',"EP1_KRT6A","EP0_MUC5B",'Unknown','EP3_MKI67',
                     'EP0_MUC5B','EP1_KRT6A','Unknown','EP4_EPCAM','EP3_MKI67','EP1_KRT6A')
names(new.cluster.ids) <- levels(Li.paired.epi)
Li.paired.epi <- RenameIdents(Li.paired.epi, new.cluster.ids)
Li.paired.epi@meta.data[["cell_subtype"]]<-Li.paired.epi@active.ident
Li.paired.epi@meta.data[["cell_subtype"]]<-factor(Li.paired.epi@meta.data[["cell_subtype"]],levels=c('EP0_MUC5B','EP1_KRT6A','EP2_POSTN',
                                                                                                     'EP3_MKI67','EP4_EPCAM','Unknown'))
Idents(Li.paired.epi)<-Li.paired.epi@meta.data[["cell_subtype"]]
DimPlot(Li.paired.epi, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.epi.anno.rds"))




Li.paired.T<-subset(Li.paired.sce,cell_type=='T/NK')
Li.paired.T <- NormalizeData(Li.paired.T, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.T <- FindVariableFeatures(Li.paired.T, selection.method = "vst", nfeatures = 3000)
Li.paired.T <- ScaleData(Li.paired.T, vars.to.regress = "percent.mt")
Li.paired.T <- RunPCA(Li.paired.T, features = VariableFeatures(object = Li.paired.T))
Li.paired.T <- RunHarmony(Li.paired.T,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.T <- IntegrateLayers(
  object = Li.paired.T, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.T)
)
ElbowPlot(Li.paired.T, ndims = 50, reduction = "harmony")
DimPlot(Li.paired.T, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.T, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
Li.paired.T <- FindNeighbors(Li.paired.T, dims = 1:30,reduction="harmony")
Li.paired.T <- FindClusters(object = Li.paired.T,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.T@meta.data, prefix = "RNA_snn_res.")
Li.paired.T <- FindClusters(Li.paired.T, resolution = 1)
Li.paired.T <- RunUMAP(Li.paired.T,reduction="harmony", dims = 1:30)
Li.paired.T <- RunTSNE(Li.paired.T,reduction="harmony", dims = 1:30)
DimPlot(Li.paired.T, reduction = "umap", label = TRUE,raster=FALSE)

Li.paired.T<-JoinLayers(Li.paired.T)
marker<-c('CD4','CCR7','TCF7','LEF1','SELL','FAS',	
          'RORC','IL17F','IL17A',
          'CXCL13','IL21','CD40LG','ICOS',
          'FOXP3','CTLA4','TIGIT','IL2RA','IKZF2',
          'STMN1','TYMS','HMGB2','UBE2C','MKI67','TOP2A','CDK1',
          'CD8A','CD8B',
          'HAVCR2','LAG3','GZMB','PDCD1',
          'GZMK','GZMA','GZMH',
          'IL7R','LTB','ITGAE',
          'KLRC1','KLRD1','FCGR3A','GNLY','PRF1','NKG7','CX3CR1',
          'KLRB1')
DotPlot(Li.paired.T, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('Tex_HAVCR2','Tex_HAVCR2','Treg_FOXP3','Naive_CD4_CCR7','Tcyto_CD8_GZMK','Tm_CD8_IL7R','Tcyto_CD8_GZMK',
                     'Tprol_MKI67','NK_FCGR3A','Th17_CD4_IL17A',"Tcyto_CD8_GZMK","Tcyto_CD8_GZMK",'NK_KLRC1','Treg_FOXP3','Tfh_CD4_CXCL13',
                     'Unknown','Unknown','Tex_HAVCR2','Unknown','Unknown','Tprol_MKI67')
names(new.cluster.ids) <- levels(Li.paired.T)
Li.paired.T <- RenameIdents(Li.paired.T, new.cluster.ids)
Li.paired.T@meta.data[["cell_subtype"]]<-Li.paired.T@active.ident
Li.paired.T@meta.data[["cell_subtype"]]<-factor(Li.paired.T@meta.data[["cell_subtype"]],levels=c('Tm_CD8_IL7R','Tcyto_CD8_GZMK','Tex_HAVCR2',
                                                                                                 'Naive_CD4_CCR7','Th17_CD4_IL17A','Tfh_CD4_CXCL13','Treg_FOXP3',
                                                                                                 'Tprol_MKI67','NK_FCGR3A','NK_KLRC1','Unknown'))
Idents(Li.paired.T)<-Li.paired.T@meta.data[["cell_subtype"]]
DimPlot(Li.paired.T, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.T, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.T.anno.rds"))



Li.paired.B<-subset(Li.paired.sce,cell_type=='B' & orig.ident!='P1_N')
Li.paired.B <- NormalizeData(Li.paired.B, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.B <- FindVariableFeatures(Li.paired.B, selection.method = "vst", nfeatures = 3000)
Li.paired.B <- ScaleData(Li.paired.B, vars.to.regress = "percent.mt")
Li.paired.B <- RunPCA(Li.paired.B, features = VariableFeatures(object = Li.paired.B))
Li.paired.B <- RunHarmony(Li.paired.B,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.B <- IntegrateLayers(
  object = Li.paired.B, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.B)
)
ElbowPlot(Li.paired.B, ndims = 50, reduction = "harmony")
DimPlot(Li.paired.B, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.B, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
Li.paired.B <- FindNeighbors(Li.paired.B, dims = 1:30,reduction="harmony")
Li.paired.B <- FindClusters(object = Li.paired.B,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.B@meta.data, prefix = "RNA_snn_res.")
Li.paired.B <- FindClusters(Li.paired.B, resolution = 0.6)
Li.paired.B <- RunUMAP(Li.paired.B,reduction="harmony", dims = 1:30)
Li.paired.B <- RunTSNE(Li.paired.B,reduction="harmony", dims = 1:30)
DimPlot(Li.paired.B, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('TNFRSF13B','CD83','BANK1','NR4A2','HLA-DRA',
          'S1PR1','BACH2','KLF4','FCER2',
          'MKI67','STMN1',
          'NEIL1','MME','AICDA','BCL6','TCL1A',
          'MZB1','XBP1','IGHA1',
          'SDC1','CD38','IGHG4')
DotPlot(Li.paired.B, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('B0_ABC_TNFRSF13B','PC_IGHA1','PC_IGHG4','B0_ABC_TNFRSF13B','B2_TC_MKI67',
                     'B1_MBC_S1PR1','B0_ABC_TNFRSF13B','B3_GCB_NEIL1','PC_IGHA1','PC_IGHA1')
names(new.cluster.ids) <- levels(Li.paired.B)
Li.paired.B <- RenameIdents(Li.paired.B, new.cluster.ids)
Li.paired.B@meta.data[["cell_subtype"]]<-Li.paired.B@active.ident
Li.paired.B@meta.data[["cell_subtype"]]<-factor(Li.paired.B@meta.data[["cell_subtype"]],levels=c('B0_ABC_TNFRSF13B','B1_MBC_S1PR1','B2_TC_MKI67',
                                                                                                 'B3_GCB_NEIL1','PC_IGHA1','PC_IGHG4'))
Idents(Li.paired.B)<-Li.paired.B@meta.data[["cell_subtype"]]
DimPlot(Li.paired.B, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.B, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.B.anno.rds"))



Li.paired.myeloid<-subset(Li.paired.sce,cell_type=='Myeloid')
Li.paired.myeloid <- NormalizeData(Li.paired.myeloid, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.myeloid <- FindVariableFeatures(Li.paired.myeloid, selection.method = "vst", nfeatures = 3000)
Li.paired.myeloid <- ScaleData(Li.paired.myeloid, vars.to.regress = "percent.mt")
Li.paired.myeloid <- RunPCA(Li.paired.myeloid, features = VariableFeatures(object = Li.paired.myeloid))
Li.paired.myeloid <- RunHarmony(Li.paired.myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.myeloid <- IntegrateLayers(
  object = Li.paired.myeloid, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.myeloid)
)
ElbowPlot(Li.paired.myeloid, ndims = 50, reduction = "harmony")
DimPlot(Li.paired.myeloid, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.myeloid, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
Li.paired.myeloid <- FindNeighbors(Li.paired.myeloid, dims = 1:30,reduction="harmony")
Li.paired.myeloid <- FindClusters(object = Li.paired.myeloid,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.myeloid@meta.data, prefix = "RNA_snn_res.")
Li.paired.myeloid <- FindClusters(Li.paired.myeloid, resolution = 1)
Li.paired.myeloid <- RunUMAP(Li.paired.myeloid,reduction="harmony", dims = 1:30)
Li.paired.myeloid <- RunTSNE(Li.paired.myeloid,reduction="harmony", dims = 1:30)
DimPlot(Li.paired.myeloid, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('FCN1','HLA-DRA',
          'C1QC','TREM2',
          'TRAC','CD247','TRBC2'
          )
marker<-c('S100A8','S100A9','VCAN',
          'HP','INHBA','IGFBP2','DEFB1','CAMP',
          'C1QA','C1QB','C1QC','TREM2','SLC40A1','MS4A4A',
          #'HLA-DPA1','HLA-DQA2',
          'TRBC2','TRAC','CCL5','CXCL13',
          'MKI67','TOP2A','TYMS',
          'CLEC9A','XCR1','CADM1','CLNK',
          'LAMP3','IDO1','CD274','CCL19',
          'CD1C','CD1E','FCER1A','CLEC10A',
          'HLA-DRA','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1')

DotPlot(Li.paired.myeloid, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('Macro_C1QC','Mono_FCN1','DC_CD1C','Macro_TRAC','Mono_FCN1','Macro_C1QC','Macro_C1QC','DC_GPR183','Macro_C1QC',
                     'DC_CD1C','DC_LAMP3','Macro_C1QC','Mono_FCN1','DC_CLEC9A','Mono_FCN1','DC_LAMP3','Macro_C1QC','DC_GPR183')
names(new.cluster.ids) <- levels(Li.paired.myeloid)
Li.paired.myeloid <- RenameIdents(Li.paired.myeloid, new.cluster.ids)
Li.paired.myeloid@meta.data[["cell_subtype"]]<-Li.paired.myeloid@active.ident
Li.paired.myeloid@meta.data[["cell_subtype"]]<-factor(Li.paired.myeloid@meta.data[["cell_subtype"]],levels=c('Mono_FCN1','Macro_C1QC','Macro_TRAC',
                                                                                                 'DC_CLEC9A','DC_LAMP3','DC_CD1C','DC_GPR183'))
Idents(Li.paired.myeloid)<-Li.paired.myeloid@meta.data[["cell_subtype"]]
DimPlot(Li.paired.myeloid, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.myeloid, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.myeloid.anno.rds"))



Li.paired.fib<-subset(Li.paired.sce,cell_type=='Fibroblasts')
Li.paired.fib <- NormalizeData(Li.paired.fib, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.fib <- FindVariableFeatures(Li.paired.fib, selection.method = "vst", nfeatures = 3000)
Li.paired.fib <- ScaleData(Li.paired.fib, vars.to.regress = "percent.mt")
Li.paired.fib <- RunPCA(Li.paired.fib, features = VariableFeatures(object = Li.paired.fib))
Li.paired.fib <- RunHarmony(Li.paired.fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.fib <- IntegrateLayers(
  object = Li.paired.fib, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.fib)
)
ElbowPlot(Li.paired.fib, ndims = 50, reduction = "harmony")
DimPlot(Li.paired.fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.fib, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
Li.paired.fib <- FindNeighbors(Li.paired.fib, dims = 1:30,reduction="harmony")
Li.paired.fib <- FindClusters(object = Li.paired.fib,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.fib@meta.data, prefix = "RNA_snn_res.")
Li.paired.fib <- FindClusters(Li.paired.fib, resolution = 0.8)
Li.paired.fib <- RunUMAP(Li.paired.fib,reduction="harmony", dims = 1:30)
Li.paired.fib <- RunTSNE(Li.paired.fib,reduction="harmony", dims = 1:30)
DimPlot(Li.paired.fib, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('MCAM','RGS5',
          'ACTA2','ACTG2','MYH11','TAGLN',
          'COL1A1','DCN','MMP2',
          'CXCL14','CXCL12',
          'IL11','IL24','IL32','CHI3L1','CXCL5',
          'SOD2','VEGFA','HGF')
DotPlot(Li.paired.fib, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('PVC0_MCAM','antiCAF_ID2','eCAF_DCN','iCAF_CXCL14','eCAF_DCN','mCAF_ACTG2','PVC1_ACTA2',
                     'PVC1_ACTA2','iCAF_CHI3L1','eCAF_DCN','antiCAF_ID2','Unknown','iCAF_CXCL14')
names(new.cluster.ids) <- levels(Li.paired.fib)
Li.paired.fib <- RenameIdents(Li.paired.fib, new.cluster.ids)
Li.paired.fib@meta.data[["cell_subtype"]]<-Li.paired.fib@active.ident
Li.paired.fib@meta.data[["cell_subtype"]]<-factor(Li.paired.fib@meta.data[["cell_subtype"]],levels=c('PVC0_MCAM','PVC1_ACTA2','mCAF_ACTG2',
                                                                                                 'eCAF_DCN','iCAF_CXCL14','iCAF_CHI3L1',
                                                                                                 'antiCAF_ID2','Unknown'))
Idents(Li.paired.fib)<-Li.paired.fib@meta.data[["cell_subtype"]]
DimPlot(Li.paired.fib, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.fib, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.fib.anno.rds"))




# Cluster without manual subtype annotation
Li.paired.endo<-subset(Li.paired.sce,cell_type=='Endothelial')
Li.paired.endo <- NormalizeData(Li.paired.endo, normalization.method = "LogNormalize", scale.factor = 10000)
Li.paired.endo <- FindVariableFeatures(Li.paired.endo, selection.method = "vst", nfeatures = 3000)
Li.paired.endo <- ScaleData(Li.paired.endo, vars.to.regress = "percent.mt")
Li.paired.endo <- RunPCA(Li.paired.endo, features = VariableFeatures(object = Li.paired.endo))
Li.paired.endo <- RunHarmony(Li.paired.endo,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Li.paired.endo <- IntegrateLayers(
  object = Li.paired.endo, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = Li.paired.endo)
)
ElbowPlot(Li.paired.endo, ndims = 50, reduction = "harmony")
DimPlot(Li.paired.endo, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Li.paired.endo, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
Li.paired.endo <- FindNeighbors(Li.paired.endo, dims = 1:30,reduction="harmony")
Li.paired.endo <- FindClusters(object = Li.paired.endo,resolution = seq(0.1,1,by=0.1))
clustree(Li.paired.endo@meta.data, prefix = "RNA_snn_res.")
Li.paired.endo <- FindClusters(Li.paired.endo, resolution = 0.8)
Li.paired.endo <- RunUMAP(Li.paired.endo,reduction="harmony", dims = 1:30)
Li.paired.endo <- RunTSNE(Li.paired.endo,reduction="harmony", dims = 1:30)
DimPlot(Li.paired.endo, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Li.paired.endo, file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.endo.anno.rds"))


####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Li.paired.epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.epi.anno.rds"))
Idents(Li.paired.epi)<-Li.paired.epi@meta.data[["cell_subtype"]]
Li.paired.epi_markers <- FindAllMarkers(Li.paired.epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Li.paired.epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_EMBO.J_2023_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Li.paired.T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.T.anno.rds"))
Idents(Li.paired.T)<-Li.paired.T@meta.data[["cell_subtype"]]
Li.paired.T<-JoinLayers(Li.paired.T)
Li.paired.T_markers <- FindAllMarkers(Li.paired.T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Li.paired.T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_EMBO.J_2023_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Li.paired.B <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.B.anno.rds"))
Idents(Li.paired.B)<-Li.paired.B@meta.data[["cell_subtype"]]
Li.paired.B<-JoinLayers(Li.paired.B)
Li.paired.B_markers <- FindAllMarkers(Li.paired.B, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Li.paired.B_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_EMBO.J_2023_B_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Li.paired.myeloid <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.myeloid.anno.rds"))
Idents(Li.paired.myeloid)<-Li.paired.myeloid@meta.data[["cell_subtype"]]
Li.paired.myeloid<-JoinLayers(Li.paired.myeloid)
Li.paired.myeloid_markers <- FindAllMarkers(Li.paired.myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Li.paired.myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_EMBO.J_2023_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Li.paired.fib <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Cao_EMBO.J_2023/Li.paired.fib.anno.rds"))
Idents(Li.paired.fib)<-Li.paired.fib@meta.data[["cell_subtype"]]
Li.paired.fib<-JoinLayers(Li.paired.fib)
Li.paired.fib_markers <- FindAllMarkers(Li.paired.fib, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Li.paired.fib_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_EMBO.J_2023_fib_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 02: Guo et al., Clinical and Translational Medicine (2023)
# Original script: single-cell/Guo_Clin.Transl.Med_2023/HuaKQ_scRNA.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- list.files(file.path(local_project_dir, "K14小鼠/Hua_scRNA_matrix"),full.names = T)
sceList <- lapply(mat,function(x){
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samples<-list.files(file.path(local_project_dir, "K14小鼠/Hua_scRNA_matrix"))
type<-c('CC','HSIL','HSIL','N_neg','N_neg','Normal','Normal','CC','CC')
orig <- NULL
location <- NULL
sce.doubletFinder.list<-list()
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  location = rep(type[i],ncol(sceList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["location"]] = location
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA <= mean(sce@meta.data[["nFeature_RNA"]])+2*sd(sce@meta.data[["nFeature_RNA"]]) & nFeature_RNA >= mean(sce@meta.data[["nFeature_RNA"]])-2*sd(sce@meta.data[["nFeature_RNA"]]) & nCount_RNA >= mean(sce@meta.data[["nCount_RNA"]])-2*sd(sce@meta.data[["nCount_RNA"]]) & nCount_RNA <= mean(sce@meta.data[["nCount_RNA"]])+2*sd(sce@meta.data[["nCount_RNA"]]) & percent.mt <= 10)
  #sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  #sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  #sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
  sce.doubletFinder.list[[i]]<-sce
  print(paste0(samples[i]," Finished!"))
}


Hua.sce<-sce.doubletFinder.list
for(i in 1:length(Hua.sce)){
  Hua.sce[[i]][["pANN"]]<-Hua.sce[[i]]@meta.data %>% select(contains('pANN'))
  Hua.sce[[i]][["DF.classifications"]]<-Hua.sce[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

Hua.sce <- merge(Hua.sce[[1]],
                y = Hua.sce[-1],
                add.cell.ids = samples)
c <- grep("pANN_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
c <- grep("DF.classifications_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
DimPlot(Hua.sce, reduction = 'pca',group.by = 'DF.classifications', label = TRUE,raster=FALSE)
Hua.sce <- subset(Hua.sce, DF.classifications== "Singlet")

saveRDS(Hua.sce, file.path(local_project_dir, "K14小鼠/数据分析/补测/Hua.sce.rds"))
#Hua.sce <- NormalizeData(Hua.sce, normalization.method = "LogNormalize")
#Hua.sce <- FindVariableFeatures(Hua.sce, selection.method = "vst", nfeatures = 3000)
#Hua.sce <- ScaleData(Hua.sce, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.sce))
Hua.sce <- SCTransform(Hua.sce, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.sce <- RunPCA(Hua.sce, features = VariableFeatures(object = Hua.sce))
Hua.sce <- RunHarmony(Hua.sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.sce, ndims = 50, reduction = "harmony")
Hua.sce <- FindNeighbors(Hua.sce, dims = 1:20,reduction="harmony")
Hua.sce <- FindClusters(object = Hua.sce,resolution = seq(0.1,1,by=0.1))
clustree(Hua.sce@meta.data, prefix = "RNA_snn_res.")
Hua.sce <- FindClusters(Hua.sce, resolution = 0.3)
Hua.sce <- RunUMAP(Hua.sce,reduction="harmony", dims = 1:20)
Hua.sce <- RunTSNE(Hua.sce,reduction="harmony", dims = 1:20)
marker<-c('EPCAM', 'CDH1', 'KRT5','TP63',             #epithelial cells
          'ITGAX', 'CSF1R', 'CD14','FCGR3A',          #myeloid cells
          'CLDN5', 'VWF','CDH5', 'KDR',               #endothelial cells
          'COL1A1', 'COL1A2', 'LUM',                  #fibroblasts
          'MZB1','CD79A', 'MS4A1',                    #B cells
          'CD3D', 'CD3E', 'CD2',                      #T cells 
          'KIT', 'IL1RL1','MS4A2')                    #mast cells
marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
marker<-c('COL1A1','DCN','COL3A1',          #fibroblasts
          'NKG7','CCL5','GZMA','CD3G', 'CD3D','CD3E',                           #T cells
          'CD68', 'CSF1R', 'CD163','LYZ',         #myeloid cells
          'NCF1','SORL1',                     #FCGR2A+monocytes
          #'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5',              #smooth muscle cells
          'PECAM1','CDH5', 'VWF',             #ECs
          'JCHAIN','MZB1', 'IGHGH1','IGHG3','IGKC','XBP1',                     #plasma cells
          'MS4A1','CD79B','CD19','CD79A','BANK1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
DotPlot(Hua.sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))
DimPlot(Hua.sce, reduction = "tsne", label = TRUE,raster=FALSE)

new.cluster.ids <- c('NK_T cells','Neutrophils','NK_T cells','NK_T cells','epithelial cells','NK_T cells','Fibroblasts','epithelial cells',
                     'myeloid cells','epithelial cells','NK_T cells','Mast cells','Plasma cells','myeloid cells','endothelial cells','Plasma cells','epithelial cells',
                     'Plasma cells','B cells','Plasma cells','epithelial cells','smooth muscle cells','epithelial cells',
                     'epithelial cells','dendritic cells','NK_T cells','epithelial cells')
names(new.cluster.ids) <- levels(Hua.sce)
Hua.sce <- RenameIdents(Hua.sce, new.cluster.ids)
Hua.sce@meta.data[["Cell_type"]]<-as.character(Hua.sce@active.ident)
Hua.sce@meta.data[["Cell_type"]]<-factor(Hua.sce@meta.data[["Cell_type"]],
                                                  levels=c('myeloid cells','NK_T cells','epithelial cells','Neutrophils',
                                                           'B cells','Mast cells','Plasma cells','smooth muscle cells',
                                                           'Fibroblasts','endothelial cells','dendritic cells'))
Idents(Hua.sce)<-Hua.sce@meta.data[["Cell_type"]]
DimPlot(Hua.sce, reduction = "umap",raster=FALSE)

saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.sce.anno.rds"))

library(AUCell)
Hua.epi<-subset(Hua.sce,Cell_type %in% c('Epithelial','Squamous_Epithelial'))
Hua.epi<-subset(Hua.sce,Cell_type =='epithelial cells')
Hua.epi <- NormalizeData(Hua.epi, normalization.method = "LogNormalize")
Hua.epi <- FindVariableFeatures(Hua.epi, selection.method = "vst", nfeatures = 3000)
Hua.epi <- ScaleData(Hua.epi, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.epi))
Hua.epi <- SCTransform(Hua.epi, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.epi <- RunPCA(Hua.epi, features = VariableFeatures(object = Hua.epi))
Hua.epi <- RunHarmony(Hua.epi,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.epi, ndims = 50, reduction = "harmony")
Hua.epi <- FindNeighbors(Hua.epi, dims = 1:30,reduction="pca")
Hua.epi <- FindClusters(object = Hua.epi,resolution = seq(0.1,1,by=0.1))
clustree(Hua.epi@meta.data, prefix = "SCT_snn_res.")
Hua.epi <- FindClusters(Hua.epi, resolution = 0.3)
Hua.epi <- RunUMAP(Hua.epi,reduction="pca", dims = 1:30)
Hua.epi <- RunTSNE(Hua.epi,reduction="pca", dims = 1:30)
DimPlot(Hua.epi, reduction = "umap",raster=FALSE,label = T)
marker<-c('SLC5A8','DERL3',
          'CDH16','CDH17','VSIG1','CTSE',
          'CASP14','PRSS27','CALML5')
DotPlot(Hua.epi, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

saveRDS(Hua.epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.epi.anno.rds"))

Idents(Hua.epi)<-factor(Hua.epi@meta.data[["location"]],levels=c("N_neg","Normal","HSIL",'CC'))
Hua_epi_markers <- FindAllMarkers(Hua.epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
genes<-c()
for(type in unique(Hua_epi_markers$cluster)){
  top_gene<-subset(Hua_epi_markers,cluster == type)[1:20,]
  #top100gene %>%  top_n(n = 100, wt = avg_log2FC) -> top100gene
  genes<-c(genes,top_gene$gene)
}
for(type in unique(Cervix_epi_markers$cluster)){
  top_gene<-subset(Cervix_epi_markers,cluster == type)[1:100,]
  #top100gene %>%  top_n(n = 100, wt = avg_log2FC) -> top100gene
  genes<-c(genes,top_gene$gene)
}
DotPlot(Cervix.epi, features = str_to_title(genes))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))


expr<-GetAssayData(JoinLayers(Hua.epi[['RNA']]), layer = 'data')
expr<-GetAssayData(Hua.epi@assays[["SCT"]], layer ="data")
expr<-AverageExpression(Hua.epi,group.by = 'location',layer='data',assays = 'SCT')
expr<-expr[['SCT']]
genes<-list()
for(type in unique(Cervix_epi_markers$cluster)){
  top100gene<-subset(Cervix_epi_markers,cluster == type)[1:100,]
  #top100gene %>%  top_n(n = 100, wt = avg_log2FC) -> top100gene
  genes[type] = list(top100=toupper(top100gene$gene))
}
AUC_result <- AUCell_run(expr, genes)
AUC_result <- getAUC(AUC_result)
AUC_result <- data.frame(t(AUC_result))
AUC_result$Group <- Hua.epi@meta.data[["seurat_clusters"]]
Hua_AUC_mean<-stats::aggregate(AUC_result[,1:5],by=list(type=AUC_result$Group),mean)
rownames(Hua_AUC_mean)<-Hua_AUC_mean$type
Hua_AUC_mean<-Hua_AUC_mean[,-1]
pheatmap(t(Hua_AUC_mean),show_colnames = T,show_rownames = T,cluster_rows=F,cluster_cols=F,scale = 'column',
         #color=colorRampPalette(rev(c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3","#F5F5F5",
         #                              "#C7EAE5","#35978F","#003C30")))(100),
         color=colorRampPalette(rev(paletteer_d("RColorBrewer::Spectral")))(200))

AUC_result<-melt(AUC_result)
colnames(AUC_result)<-c('Stage','Type','Score')
AUC_result$Stage<-factor(AUC_result$Stage,levels = c('N_neg','Normal','HSIL','CC'))
ggplot(data = AUC_result, aes(x = Stage, y = Score)) +
  geom_boxplot(aes(color = Stage), width = 0.8,size=0.7) +
  #geom_jitter(aes(color = Stage), width = 0.5, size = 1,alpha=0.8)+
  scale_color_nejm()+
  scale_fill_nejm()+
  labs(x='')+
  theme_classic()+
  labs(x="",y="Geneset Scores")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  stat_compare_means()+
  facet_wrap( ~ Type, scales = "free_y")

gsva_scores <- gsvaParam(expr,genes)
gsva_scores <- gsva(gsva_scores)
gsva_scores <- as.data.frame(t(gsva_scores))
gsva_scores <- gsva_scores[colnames(Hua.epi),]
gsva_scores$Group <- Hua.epi@meta.data[["orig.ident"]]
Hua_AUC_mean<-stats::aggregate(gsva_scores[,1:5],by=list(type=gsva_scores$Group),mean)
rownames(Hua_AUC_mean)<-Hua_AUC_mean$type
Hua_AUC_mean<-Hua_AUC_mean[,-1]
pheatmap(t(Hua_AUC_mean),show_colnames = T,show_rownames = T,cluster_rows=F,cluster_cols=F,#scale = 'column',
         #color=colorRampPalette(rev(c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3","#F5F5F5",
         #                              "#C7EAE5","#35978F","#003C30")))(100),
         color=colorRampPalette(rev(paletteer_d("RColorBrewer::Spectral")))(200))

gsva_scores<-melt(gsva_scores)
colnames(gsva_scores)<-c('Stage','Type','Score')
gsva_scores$Stage<-factor(gsva_scores$Stage,levels = c('N_neg','Normal','HSIL','CC'))
ggplot(data = gsva_scores, aes(x = Stage, y = Score)) +
  geom_boxplot(aes(color = Stage), width = 0.8,size=0.7) +
  #geom_jitter(aes(color = Stage), width = 0.5, size = 1,alpha=0.8)+
  scale_color_nejm()+
  scale_fill_nejm()+
  labs(x='')+
  theme_classic()+
  labs(x="",y="Geneset Scores")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.title=element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  stat_compare_means()+
  facet_wrap( ~ Type, scales = "free_y")


genes<-list()
gsva_result=list()
for (top_num in c(30,50,75,100,150,200,250,300)){
  for(type in unique(Cervix_epi_markers$cluster)){
    top_gene<-subset(Cervix_epi_markers,cluster == type)
    top_gene %>% group_by(cluster) %>% top_n(n = top_num, wt = avg_log2FC) -> top_gene
    #order<-order(top_gene$p_val_adj,-top_gene$avg_log2FC)[1:top_num]
    #top_gene<-top_gene[order,]$gene
    genes[type] = list(top=toupper(top_gene$gene))
  }
  #AUC_result <- AUCell_run(expr, genes)
  #AUC_result <- getAUC(AUC_result)
  #AUC_result <- data.frame(t(AUC_result))
  #AUC_result$Group <- Huang_epi@meta.data[["cell_type_disease"]]
  #AUC_result$Group <- Huang_epi@meta.data[["orig.ident"]]
  #Huang_AUC_mean<-stats::aggregate(AUC_result[,1:5],by=list(type=AUC_result$Group),mean)
  #rownames(Huang_AUC_mean)<-Huang_AUC_mean$type
  gsva_scores <- gsvaParam(expr,genes)
  gsva_scores <- gsva(gsva_scores)
  gsva_result[paste0('Top',top_num)] = list(score=gsva_scores)
  save(gsva_result,file='gsva_result.RData')
  gsva_scores <- data.frame(t(gsva_scores))
  gsva_scores <- gsva_scores[colnames(Hua.epi),]
  gsva_scores$Stage<-Hua.epi@meta.data[["location"]]
  gsva_scores_mean<-stats::aggregate(gsva_scores[,1:5],by=list(type=gsva_scores$Stage),mean)
  rownames(gsva_scores_mean)<-gsva_scores_mean$type
  gsva_scores_mean<-gsva_scores_mean[c('N_neg','Normal','HSIL','CC'),-1]
  pdf(file = paste0('Hua_gsva_heatmap_Top',top_num,'.pdf'),width = 10,height = 7)
  p<-pheatmap(t(gsva_scores_mean),show_colnames = T,show_rownames = T,cluster_rows=F,cluster_cols=F,scale = 'column',
           color=colorRampPalette(rev(paletteer_d("RColorBrewer::Spectral")))(200))
  print(p)
  dev.off()
  
  
  gsva_scores<-melt(gsva_scores,variable.name = 'Stage')
  colnames(gsva_scores)<-c('Stage','Type','Score')
  gsva_scores$Stage<-factor(gsva_scores$Stage,levels = c('N_neg','Normal','HSIL','CC'))
  pdf(file = paste0('Hua_gsva_boxplot_Top',top_num,'.pdf'),width = 10,height = 7)
  p<-ggplot(data = gsva_scores, aes(x = Stage, y = Score)) +
    geom_boxplot(aes(color = Stage), width = 0.8,size=0.7) +
    #geom_jitter(aes(color = Stage), width = 0.5, size = 1,alpha=0.8)+
    scale_color_nejm()+
    scale_fill_nejm()+
    labs(x='')+
    theme_classic()+
    labs(x="",y="Geneset Scores")+
    guides(fill=guide_legend(title=NULL))+
    theme(legend.title=element_blank())+
    theme(axis.line = element_line(colour = "black"))+
    stat_compare_means()+
    facet_wrap( ~ Type, scales = "free_y")
  print(p)
  dev.off()
}

DotPlot(Hua.epi, features = marker,group.by='location')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))



# T-cell subclustering
Hua.NKT<-subset(Hua.sce,Cell_type =='NK_T cells')
Hua.NKT <- SCTransform(Hua.NKT, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.NKT <- RunPCA(Hua.NKT, features = VariableFeatures(object = Hua.NKT))
Hua.NKT <- RunHarmony(Hua.NKT,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.NKT, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.NKT, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.NKT, ndims = 50, reduction = "harmony")
Hua.NKT <- FindNeighbors(Hua.NKT, dims = 1:30,reduction="harmony")
Hua.NKT <- FindClusters(object = Hua.NKT,resolution = seq(0.1,1,by=0.1))
clustree(Hua.NKT@meta.data, prefix = "SCT_snn_res.")
Hua.NKT <- FindClusters(Hua.NKT, resolution = 0.7)
Hua.NKT <- RunUMAP(Hua.NKT,reduction="harmony", dims = 1:30)
DimPlot(Hua.NKT, reduction = "umap",raster=FALSE,label = T)

marker<-c('GZMK','CXCR4','CST7',    #CD8 Tem
          'XCL1','CAPG','NR4A1',    #Trm
          'CD160','KIR2DL4','KLRC2', #IEL
          'HAVCR2','CXCL13','PDCD1', #Tex
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL',  #naive T
          'FOXP3','CTLA4','TNFRSF18', #CD4 Treg
          'CCL5',   #Tems
          'IL17A','CCR6','CTSH',  #Th17
          'IFNG', #Th1-like
          'CD4','CD8A','CD8B','NKG7'
          )
DotPlot(Hua.NKT, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))



new.cluster.ids <- c('CD4','CD8','CD8','CD8','CD4','CD4','CD8','CD4','CD8',
                     'CD4','CD4','CD8','CD8','CD8','CD8','CD8','CD4','CD4')
names(new.cluster.ids) <- levels(Hua.NKT)
Hua.NKT <- RenameIdents(Hua.NKT, new.cluster.ids)
Hua.NKT@meta.data[["Cell_subtype"]]<-as.character(Hua.NKT@active.ident)
Idents(Hua.NKT)<-Hua.NKT@meta.data[["Cell_subtype"]]
DimPlot(Hua.NKT, reduction = "umap",raster=FALSE)

FeaturePlot(Hua.NKT,features = c('CD4','CD8A','CD8B','NKG7'),raster=FALSE)


#CD4
Hua.CD4<-subset(Hua.NKT,Cell_subtype =='CD4')
Hua.CD4 <- RunPCA(Hua.CD4, features = VariableFeatures(object = Hua.CD4))
Hua.CD4 <- RunHarmony(Hua.CD4,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.CD4, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.CD4, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.CD4, ndims = 50, reduction = "harmony")
Hua.CD4 <- FindNeighbors(Hua.CD4, dims = 1:30,reduction="harmony")
Hua.CD4 <- FindClusters(object = Hua.CD4,resolution = seq(0.1,1,by=0.1))
clustree(Hua.CD4@meta.data, prefix = "SCT_snn_res.")
Hua.CD4 <- FindClusters(Hua.CD4, resolution = 0.3)
Hua.CD4 <- RunUMAP(Hua.CD4,reduction="harmony", dims = 1:30)
Hua.CD4 <- RunTSNE(Hua.CD4,reduction="harmony", dims = 1:30)
DimPlot(Hua.CD4, reduction = "umap",raster=FALSE,label = T)

marker<-c('FOXP3','CTLA4','TNFRSF18', #CD4 Treg
          'CCL5','GZMK','CXCR4',   #Tems
          'LEF1','CCR7','SELL',   #naive
          'IL17A','CCR6','CTSH',  #Th17
          'CXCL13','PDCD1','IFNG'  #Th1-like
)
DotPlot(Hua.CD4, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))


new.cluster.ids <- c('Tem','Treg','Tem','Naive','Tem','Th1-like','Tem','Th17','Treg')
names(new.cluster.ids) <- levels(Hua.CD4)
Hua.CD4 <- RenameIdents(Hua.CD4, new.cluster.ids)
Hua.CD4@meta.data[["Cell_subtype"]]<-as.character(Hua.CD4@active.ident)
Hua.CD4@meta.data[["Cell_subtype"]]<-factor(Hua.CD4@meta.data[["Cell_subtype"]],
                                         levels=c('Tem','Th1-like','Th17','Naive','Treg'))
Idents(Hua.CD4)<-Hua.CD4@meta.data[["Cell_subtype"]]
DimPlot(Hua.CD4, reduction = "umap",raster=FALSE)

saveRDS(Hua.CD4, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.CD4T.anno.rds"))

#CD8
Hua.CD8<-subset(Hua.NKT,Cell_subtype =='CD8')
Hua.CD8 <- RunPCA(Hua.CD8, features = VariableFeatures(object = Hua.CD8))
Hua.CD8 <- RunHarmony(Hua.CD8,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.CD8, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.CD8, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.CD8, ndims = 50, reduction = "harmony")
Hua.CD8 <- FindNeighbors(Hua.CD8, dims = 1:30,reduction="harmony")
Hua.CD8 <- FindClusters(object = Hua.CD8,resolution = seq(0.1,1,by=0.1))
clustree(Hua.CD8@meta.data, prefix = "SCT_snn_res.")
Hua.CD8 <- FindClusters(Hua.CD8, resolution = 0.7)
Hua.CD8 <- RunUMAP(Hua.CD8,reduction="harmony", dims = 1:30)
Hua.CD8 <- RunTSNE(Hua.CD8,reduction="harmony", dims = 1:30)
DimPlot(Hua.CD8, reduction = "umap",raster=FALSE,label = T)

marker<-c('GZMK','CXCR4','CST7',    #CD8 Tem
          'XCL1','CAPG','NR4A1',    #Trm
          'CD160','KIR2DL4','KLRC2', #IEL
          'HAVCR2','CXCL13','PDCD1', #Tex
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL',  #naive T
          'CD4','CD8A','CD8B','NKG7'
)
DotPlot(Hua.CD8, features = c('XCL1','CAPG','NR4A1'))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))


new.cluster.ids <- c('Tem','Tem','Trm','Tem','TemRA','Tem','Trm','Tem','Trm','IEL','Naive','MAIT','IEL','Tex')
names(new.cluster.ids) <- levels(Hua.CD8)
Hua.CD8 <- RenameIdents(Hua.CD8, new.cluster.ids)
Hua.CD8@meta.data[["Cell_subtype"]]<-as.character(Hua.CD8@active.ident)
Hua.CD8@meta.data[["Cell_subtype"]]<-factor(Hua.CD8@meta.data[["Cell_subtype"]],
                                            levels=c('MAIT','Trm','Naive','Tex','IEL','Tem','TemRA'))
Idents(Hua.CD8)<-Hua.CD8@meta.data[["Cell_subtype"]]
DimPlot(Hua.CD8, reduction = "umap",raster=FALSE)

saveRDS(Hua.CD8, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.CD8T.anno.rds"))


# Myeloid-cell subclustering
Hua.myeloid<-subset(Hua.sce,Cell_type =='myeloid cells')
Hua.myeloid <- RunPCA(Hua.myeloid, features = VariableFeatures(object = Hua.myeloid))
Hua.myeloid <- RunHarmony(Hua.myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.myeloid, ndims = 50, reduction = "harmony")
Hua.myeloid <- FindNeighbors(Hua.myeloid, dims = 1:30,reduction="harmony")
Hua.myeloid <- FindClusters(object = Hua.myeloid,resolution = seq(0.1,1,by=0.1))
clustree(Hua.myeloid@meta.data, prefix = "SCT_snn_res.")
Hua.myeloid <- FindClusters(Hua.myeloid, resolution = 0.8)
Hua.myeloid <- RunUMAP(Hua.myeloid,reduction="harmony", dims = 1:30)
Hua.myeloid <- RunTSNE(Hua.myeloid,reduction="harmony", dims = 1:30)
DimPlot(Hua.myeloid, reduction = "tsne",raster=FALSE,label = T)

marker<-c('XCR1','CLEC9A',
          'CD207','CD1C','FCER1A',
          'CLEC4C','LILRA4')
DotPlot(Hua.myeloid, features = c('TOP2A','MKI67'))+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('Monocytes','Macrophages','cDC2','cDC2','Macrophages','cDC2','Macrophages','Monocytes','Macrophages',
                     'cDC1','Macrophages','Macrophages','cDC2','pDC')
names(new.cluster.ids) <- levels(Hua.myeloid)
Hua.myeloid <- RenameIdents(Hua.myeloid, new.cluster.ids)
Hua.myeloid@meta.data[["Cell_subtype"]]<-as.character(Hua.myeloid@active.ident)
Hua.myeloid@meta.data[["Cell_subtype"]]<-factor(Hua.myeloid@meta.data[["Cell_subtype"]],
                                            levels=c('MAIT','Trm','Naive','Tex','IEL','Tem','TemRA'))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
DimPlot(Hua.myeloid, reduction = "umap",raster=FALSE)
saveRDS(Hua.myeloid, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.myeloid.anno.rds"))


# Cluster without manual subtype annotation
Hua.neu<-subset(Hua.sce,Cell_type =='Neutrophils')
Hua.neu <- RunPCA(Hua.neu, features = VariableFeatures(object = Hua.neu))
Hua.neu <- RunHarmony(Hua.neu,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.neu, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.neu, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.neu, ndims = 50, reduction = "harmony")
Hua.neu <- FindNeighbors(Hua.neu, dims = 1:30,reduction="harmony")
Hua.neu <- FindClusters(object = Hua.neu,resolution = seq(0.1,1,by=0.1))
clustree(Hua.neu@meta.data, prefix = "SCT_snn_res.")
Hua.neu <- RunUMAP(Hua.neu,reduction="harmony", dims = 1:30)
Hua.neu <- RunTSNE(Hua.neu,reduction="harmony", dims = 1:30)
DimPlot(Hua.neu, reduction = "umap",raster=FALSE,label = T)
saveRDS(Hua.neu, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.neutrophils.rds"))


Hua.fib<-subset(Hua.sce,Cell_type =='Fibroblasts')
Hua.fib <- RunPCA(Hua.fib, features = VariableFeatures(object = Hua.fib))
Hua.fib <- RunHarmony(Hua.fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.fib, ndims = 50, reduction = "harmony")
Hua.fib <- FindNeighbors(Hua.fib, dims = 1:30,reduction="harmony")
Hua.fib <- FindClusters(object = Hua.fib,resolution = seq(0.1,1,by=0.1))
clustree(Hua.fib@meta.data, prefix = "SCT_snn_res.")
Hua.fib <- RunUMAP(Hua.fib,reduction="harmony", dims = 1:30)
Hua.fib <- RunTSNE(Hua.fib,reduction="harmony", dims = 1:30)
DimPlot(Hua.fib, reduction = "umap",raster=FALSE,label = T)
saveRDS(Hua.fib, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.fibroblasts.rds"))


Hua.endo<-subset(Hua.sce,Cell_type =='endothelial cells')
Hua.endo <- RunPCA(Hua.endo, features = VariableFeatures(object = Hua.endo))
Hua.endo <- RunHarmony(Hua.endo,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.endo, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.endo, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.endo, ndims = 50, reduction = "harmony")
Hua.endo <- FindNeighbors(Hua.endo, dims = 1:30,reduction="harmony")
Hua.endo <- FindClusters(object = Hua.endo,resolution = seq(0.1,1,by=0.1))
clustree(Hua.endo@meta.data, prefix = "SCT_snn_res.")
Hua.endo <- RunUMAP(Hua.endo,reduction="harmony", dims = 1:30)
Hua.endo <- RunTSNE(Hua.endo,reduction="harmony", dims = 1:30)
DimPlot(Hua.endo, reduction = "umap",raster=FALSE,label = T)
saveRDS(Hua.endo, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.endothelial.rds"))


Hua.B<-subset(Hua.sce,Cell_type %in% c('B cells','Plasma cells'))
Hua.B <- RunPCA(Hua.B, features = VariableFeatures(object = Hua.B))
Hua.B <- RunHarmony(Hua.B,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.B, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.B, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.B, ndims = 50, reduction = "harmony")
Hua.B <- FindNeighbors(Hua.B, dims = 1:30,reduction="harmony")
Hua.B <- FindClusters(object = Hua.B,resolution = seq(0.1,1,by=0.1))
clustree(Hua.B@meta.data, prefix = "SCT_snn_res.")
Hua.B <- RunUMAP(Hua.B,reduction="harmony", dims = 1:30)
Hua.B <- RunTSNE(Hua.B,reduction="harmony", dims = 1:30)
DimPlot(Hua.B, reduction = "umap",raster=FALSE,label = T)
saveRDS(Hua.B, file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.B.Plasma.rds"))


####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Hua.epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.epi.anno.rds"))
Idents(Hua.epi)<-Hua.epi@meta.data[["SCT_snn_res.0.3"]]
Hua.epi_markers <- FindAllMarkers(Hua.epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Guo_Clin.Transl.Med_2023_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.CD4T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.CD4T.anno.rds"))
Idents(Hua.CD4T)<-Hua.CD4T@meta.data[["Cell_subtype"]]
Hua.CD4T_markers <- FindAllMarkers(Hua.CD4T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.CD4T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Guo_Clin.Transl.Med_2023_CD4T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.CD8T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.CD8T.anno.rds"))
Idents(Hua.CD8T)<-Hua.CD8T@meta.data[["Cell_subtype"]]
Hua.CD8T_markers <- FindAllMarkers(Hua.CD8T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.CD8T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Guo_Clin.Transl.Med_2023_CD8T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.myeloid <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Guo_Clin.Transl.Med_2023/Hua.myeloid.anno.rds"))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
Hua.myeloid_markers <- FindAllMarkers(Hua.myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Guo_Clin.Transl.Med_2023_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 03: Li et al., Cancer Cell International (2025)
# Original script: single-cell/Li_Cancer.Cell.Int_2025/Zhengzhou_for_DB.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- c(file.path(local_project_dir, "郑大单细胞数据/CESC/PN1/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PN5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PT1/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PT2/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PT3/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PT4/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PT5/filtered_feature_bc_matrix"),
         file.path(local_project_dir, "郑大单细胞数据/CESC/PT6/filtered_feature_bc_matrix"))
sceList <- lapply(mat,function(x){ 
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samples<-c('PN1','PN5','PT1','PT2','PT3','PT4','PT5','PT6')
type<-c('Normal','Normal','Tumor','Tumor','Tumor','Tumor','Tumor','Tumor')
orig <- NULL
tissue <- NULL
CC.sce<-list()
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  tissue = rep(type[i],ncol(sceList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["tissue"]] = tissue
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(sce@assays$RNA)) # Match hemoglobin genes to expression-matrix rows
  HB.genes <- rownames(sce@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]  # Keep matched genes and remove missing values
  sce[["percent.HB"]] <- PercentageFeatureSet(sce, features=HB.genes)
  if (samples[i]  %in% c('PN1','PN5')){
    sce <- subset(sce, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & nCount_RNA >1000 & nCount_RNA < 25000 & percent.mt < 15 & percent.HB < 5)
  } else if (samples[i] == 'PT1'){
    sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & nCount_RNA >1000 & nCount_RNA < 50000 & percent.mt < 15 & percent.HB < 5)
  } else if (samples[i] == 'PT2'){
    sce <- subset(sce, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & nCount_RNA >500 & nCount_RNA < 40000 & percent.mt < 15 & percent.HB < 5)
  } else if (samples[i] == 'PT3'){
    sce <- subset(sce, subset = nFeature_RNA > 600 & nFeature_RNA < 5500 & nCount_RNA >1000 & nCount_RNA < 40000 & percent.mt < 15 & percent.HB < 5)
  } else if (samples[i] == 'PT4'){
    sce <- subset(sce, subset = nFeature_RNA > 600 & nFeature_RNA < 8000 & nCount_RNA >1000 & nCount_RNA < 50000 & percent.mt < 25 & percent.HB < 5)
  } else if (samples[i] == 'PT5'){
    sce <- subset(sce, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & nCount_RNA >1000 & nCount_RNA < 40000 & percent.mt < 15 & percent.HB < 5)
  } else if (samples[i] == 'PT6'){
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
                add.cell.ids = samples)
CC.sce <- NormalizeData(CC.sce, normalization.method = "LogNormalize", scale.factor = 10000)
CC.sce <- FindVariableFeatures(CC.sce, selection.method = "vst", nfeatures = 2000)
CC.sce <- ScaleData(CC.sce, vars.to.regress = "percent.mt")
CC.sce <- RunPCA(CC.sce, features = VariableFeatures(object = CC.sce))
CC.sce <- IntegrateLayers(
  object = CC.sce, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",features = VariableFeatures(object = CC.sce)
)
ElbowPlot(CC.sce, ndims = 50, reduction = "pca")
DimPlot(CC.sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(CC.sce, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
CC.sce <- FindNeighbors(CC.sce, dims = 1:30,reduction="cca")
CC.sce <- FindClusters(object = CC.sce,resolution = seq(0.1,1,by=0.1))
clustree(CC.sce@meta.data, prefix = "RNA_snn_res.")
CC.sce <- FindClusters(CC.sce, resolution = 0.7)
CC.sce <- RunUMAP(CC.sce,reduction="cca", dims = 1:30)
CC.sce <- RunTSNE(CC.sce,reduction="cca", dims = 1:30)
DimPlot(CC.sce, reduction = "tsne", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2', 'KRT8','KRT18','KRT19','TSPAN8','MUC5B',   #epithelial cells
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CDH5', 'EMCN','PECAM1', 'PCDH17','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','MZB1','CD38',                      #plasma cells
          'MS4A1','CD79B','IGKC','CD79A',                     #B cells
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CPA3','TPSAB1','KIT')                     #mast cells
DotPlot(CC.sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('Epithelial cells','Fibroblasts','Epithelial cells','T cells','T cells','Epithelial cells','Endothelial cells',
                     'Fibroblasts','T cells','Myeloid cells','Epithelial cells','Fibroblasts','Fibroblasts','Fibroblasts','B cells',
                     'T cells','T cells','Epithelial cells','Plasma cells','Myeloid cells','Mast cells','Myeloid cells','Endothelial cells','Epithelial cells')
names(new.cluster.ids) <- levels(CC.sce)
CC.sce <- RenameIdents(CC.sce, new.cluster.ids)
CC.sce@meta.data[["cell_type"]]<-CC.sce@active.ident
CC.sce@meta.data[["cell_type"]]<-factor(CC.sce@meta.data[["cell_type"]],levels=c('Epithelial cells','Fibroblasts','Endothelial cells',
                                                                                 'T cells','Myeloid cells','B cells','Plasma cells','Mast cells'))
Idents(CC.sce)<-CC.sce@meta.data[["cell_type"]]
DimPlot(CC.sce, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(CC.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.sce.anno.rds"))


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
ElbowPlot(CC.epi, ndims = 50, reduction = "cca")
DimPlot(CC.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(CC.epi, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
CC.epi <- FindNeighbors(CC.epi, dims = 1:30,reduction="cca")
CC.epi <- FindClusters(object = CC.epi,resolution = seq(0.1,1,by=0.1))
clustree(CC.epi@meta.data, prefix = "RNA_snn_res.")
CC.epi <- FindClusters(CC.epi, resolution = 0.6)
CC.epi <- RunUMAP(CC.epi,reduction="cca", dims = 1:30)
CC.epi <- RunTSNE(CC.epi,reduction="cca", dims = 1:30)
DimPlot(CC.epi, reduction = "umap", label = TRUE,raster=FALSE,split.by = 'tissue')

marker<-c('HIST1H1D','HIST1H1B','RRM2','HIST1H4C','UBE2C','ASPM','TOP2A','CCNB1','CENPF','BPIFB1','SCGB2A1','PIGR','MUC5B','TFF3',
          'CLSPN','PCNA','TK1','TYMS','PCLAF','S100A7','KETDAP','RHCG','SBSN','CEACAM6','MUC20','CAV1','CCN1','COL17A1','DST','CXCL14',
          'SPRR3','SPRR1B','SLPI','LCN2','IGFL1')
DotPlot(CC.epi, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('C2-DST-NEpis','C6-TFF3-IAEpis','C1-IGFL1-Epis','C2-DST-NEpis','C2-DST-NEpis','C5-PCLAF-TAEpis','C4-KRTDAP-IAEpis',
                     'C3-MUC20-Tu','C5-PCLAF-TAEpis','C2-DST-NEpis','C8-NEURL1B-TAEpis','C2-DST-NEpis','C7-CENPF-TAEpis','C7-CENPF-TAEpis','C6-TFF3-IAEpis')
names(new.cluster.ids) <- levels(CC.epi)
CC.epi <- RenameIdents(CC.epi, new.cluster.ids)
CC.epi@meta.data[["cell_subtype"]]<-CC.epi@active.ident
CC.epi@meta.data[["cell_subtype"]]<-factor(CC.epi@meta.data[["cell_subtype"]],levels=c('C1-IGFL1-Epis','C2-DST-NEpis','C3-MUC20-Tu','C4-KRTDAP-IAEpis',
                                                                                    'C5-PCLAF-TAEpis','C6-TFF3-IAEpis','C7-CENPF-TAEpis','C8-NEURL1B-TAEpis'))
Idents(CC.epi)<-CC.epi@meta.data[["cell_subtype"]]
DimPlot(CC.epi, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(CC.epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.epi.anno.rds"))



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
ElbowPlot(CC.fib, ndims = 50, reduction = "pca")
DimPlot(CC.fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(CC.fib, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
CC.fib <- FindNeighbors(CC.fib, dims = 1:30,reduction="cca")
CC.fib <- FindClusters(object = CC.fib,resolution = seq(0.1,1,by=0.1))
clustree(CC.fib@meta.data, prefix = "RNA_snn_res.")
CC.fib <- FindClusters(CC.fib, resolution = 1)
CC.fib <- RunUMAP(CC.fib,reduction="cca", dims = 1:30)
CC.fib <- RunTSNE(CC.fib,reduction="cca", dims = 1:30)
DimPlot(CC.fib, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('CNN1','PCP4','ACTG2','DES',
          'NOTCH3','HIGD1B','STEAP4','NDUFA4L2','RGS5',
          'ADIRF','ACTA2','ECRG4','MYH11','MUSTN1',
          'CTHRC1','COL5A2','WNT5A','SFRP2','MMP11',
          'CFD','PTGDS','SFRP1','SFRP4','CCN5')
DotPlot(CC.fib, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('C2-MMP11-CAFs','C1-SPER4-IAFs','C4-RGS5-pericytes','C3-MUSTN1-myofibro','C1-SPER4-IAFs','C1-SPER4-IAFs',
                     'C2-MMP11-CAFs','C2-MMP11-CAFs','C5-DES-myofibro','C1-SPER4-IAFs','C1-SPER4-IAFs','C1-SPER4-IAFs',
                     'C3-MUSTN1-myofibro','C1-SPER4-IAFs','C3-MUSTN1-myofibro','C2-MMP11-CAFs','C1-SPER4-IAFs','C2-MMP11-CAFs')
names(new.cluster.ids) <- levels(CC.fib)
CC.fib <- RenameIdents(CC.fib, new.cluster.ids)
CC.fib@meta.data[["cell_subtype"]]<-CC.fib@active.ident
CC.fib@meta.data[["cell_subtype"]]<-factor(CC.fib@meta.data[["cell_subtype"]],levels=c('C1-SPER4-IAFs','C2-MMP11-CAFs','C3-MUSTN1-myofibro',
                                                                                       'C4-RGS5-pericytes','C5-DES-myofibro'))
Idents(CC.fib)<-CC.fib@meta.data[["cell_subtype"]]
DimPlot(CC.fib, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(CC.fib, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.fib.anno.rds"))




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
ElbowPlot(CC.ec, ndims = 50, reduction = "pca")
DimPlot(CC.ec, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(CC.ec, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
CC.ec <- FindNeighbors(CC.ec, dims = 1:30,reduction="cca")
CC.ec <- FindClusters(object = CC.ec,resolution = seq(0.1,1,by=0.1))
clustree(CC.ec@meta.data, prefix = "RNA_snn_res.")
CC.ec <- FindClusters(CC.ec, resolution = 0.9)
CC.ec <- RunUMAP(CC.ec,reduction="cca", dims = 1:30)
CC.ec <- RunTSNE(CC.ec,reduction="cca", dims = 1:30)
DimPlot(CC.ec, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('EFNB2','ARL15','IGFBP3','CXCL12','SEMA3G',
          'TCIM','APOLD1','INSR','COL4A2','RGCC','COL4A1',
          'C7','MMRN1','CCL14','CLU','ACKR1','CCL21')
DotPlot(CC.ec, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('C2-EDNRB-capillary ECs','C1-ACKR1-Venous ECs','C1-ACKR1-Venous ECs','C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs',
                     'C4-FBLBN5-arterial ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C2-EDNRB-capillary ECs','C3-CCL21-lymphatic ECs')
new.cluster.ids <- c('C1-ACKR1-Venous ECs','C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs','C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs',
                     'C2-EDNRB-capillary ECs','C4-FBLBN5-arterial ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C3-CCL21-lymphatic ECs','C2-EDNRB-capillary ECs','C3-CCL21-lymphatic ECs')
names(new.cluster.ids) <- levels(CC.ec)
CC.ec <- RenameIdents(CC.ec, new.cluster.ids)
CC.ec@meta.data[["cell_subtype"]]<-CC.ec@active.ident
CC.ec@meta.data[["cell_subtype"]]<-factor(CC.ec@meta.data[["cell_subtype"]],levels=c('C1-ACKR1-Venous ECs','C2-EDNRB-capillary ECs',
                                                                                     'C3-CCL21-lymphatic ECs','C4-FBLBN5-arterial ECs'))
Idents(CC.ec)<-CC.ec@meta.data[["cell_subtype"]]
DimPlot(CC.ec, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(CC.ec, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.ECs.anno.rds"))



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
ElbowPlot(CC.mye, ndims = 50, reduction = "pca")
DimPlot(CC.mye, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(CC.mye, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
CC.mye <- FindNeighbors(CC.mye, dims = 1:30,reduction="cca")
CC.mye <- FindClusters(object = CC.mye,resolution = seq(0.1,1,by=0.1))
clustree(CC.mye@meta.data, prefix = "RNA_snn_res.")
CC.mye <- FindClusters(CC.mye, resolution = 0.6)
CC.mye <- RunUMAP(CC.mye,reduction="cca", dims = 1:30)
CC.mye <- RunTSNE(CC.mye,reduction="cca", dims = 1:30)
DimPlot(CC.mye, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('FSCN1','CST7','CRIP1','LAMP3','IDO1',
          'RSAD2','IFIT1','IFIT2','IFIT3','CSF3R','HCAR3',
          'FCGBP','PKIB','FCER1A','HLA-DQB2','S100B',
          'SLC40A1','FOLR2','RNASE1','SELENOP',
          'IFITM2','SMIM25','FCGR3B','CXCR4','CXCL8',
          'VCAN','S100A12',
          'S100A9','S100A8','FCN1',
          'SPP1','APOE','C1QA','APOC1','C1QB',
          'CD163','CD1C','ISG15')
DotPlot(CC.mye, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('C7-ISG15-Neus','C4-CD163-TAMs','C2-S100A8-IANs','C5-CD1C-cDC2','C1-C1QA-Macro','C6-CXCL8-TANs',
                     'C5-CD1C-cDC2','C5-CD1C-cDC2','C3-CXCR4-TANs','C8-LAMP3-cDC3','C1-C1QA-Macro','C1-C1QA-Macro','C8-LAMP3-cDC3','C5-CD1C-cDC2')
names(new.cluster.ids) <- levels(CC.mye)
CC.mye <- RenameIdents(CC.mye, new.cluster.ids)
CC.mye@meta.data[["cell_subtype"]]<-CC.mye@active.ident
CC.mye@meta.data[["cell_subtype"]]<-factor(CC.mye@meta.data[["cell_subtype"]],levels=c('C1-C1QA-Macro','C2-S100A8-IANs',
                                                                                     'C3-CXCR4-TANs','C4-CD163-TAMs','C5-CD1C-cDC2',
                                                                                     'C6-CXCL8-TANs','C7-ISG15-Neus','C8-LAMP3-cDC3'))
Idents(CC.mye)<-CC.mye@meta.data[["cell_subtype"]]
DimPlot(CC.mye, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(CC.mye, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.myeloid.anno.rds"))


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
ElbowPlot(CC.NKT, ndims = 50, reduction = "pca")
DimPlot(CC.NKT, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(CC.NKT, group.by = 'orig.ident', reduction = "cca",raster=FALSE)
CC.NKT <- FindNeighbors(CC.NKT, dims = 1:30,reduction="cca")
CC.NKT <- FindClusters(object = CC.NKT,resolution = seq(0.1,1,by=0.1))
clustree(CC.NKT@meta.data, prefix = "RNA_snn_res.")
CC.NKT <- FindClusters(CC.NKT, resolution = 0.8)
CC.NKT <- RunUMAP(CC.NKT,reduction="cca", dims = 1:30)
CC.NKT <- RunTSNE(CC.NKT,reduction="cca", dims = 1:30)
DimPlot(CC.NKT, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('CD4','CD8A','CD8B','NKG7')
marker<-c('CD8A','CCL5','CD8B','GZMB','XCL1','IL7R','CTSW','HLA-DRB1','KLRB1','IFNG','GZMK','CST7','GZMA','CCL4','CCL4L2','TNFRSF4',
          'BATF','TNFRSF18','FOXP3','CTLA4','GNLY','TYROBP','TRDC','FCER1G','KLRD1','CXCL13','KRT86','CCL3','LINC00892','CD40LG',
          'TIGIT','PPP1R2C','LINC02195','IL2RA','PDCD1','RBPJ','DUSP4',
          'ZNF683','GPR183','HIST1H1B')
DotPlot(CC.NKT, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = rev(paletteer_d("RColorBrewer::Spectral")))

# Manual cluster annotation
new.cluster.ids <- c('CD8-C1-ZNF683-Trm','CD4-C2-IL7R-Tcm','CD8-C3-GZMK-Tem','CD4-C4-FOXP3-Treg','CD8-C1-ZNF683-Trm','CD4-C2-IL7R-Tcm',
                     'CD8-C6-CXCL13-Tex','NK-C5-NKG7','CD8-C8-GPR183-Tcm','DP-C11-HIST1H1B','CD4-C10-CXCL13-Th1','NK-C5-NKG7','DP-C11-HIST1H1B')
names(new.cluster.ids) <- levels(CC.NKT)
CC.NKT <- RenameIdents(CC.NKT, new.cluster.ids)
CC.NKT@meta.data[["cell_subtype"]]<-CC.NKT@active.ident
CC.NKT@meta.data[["cell_subtype"]]<-factor(CC.NKT@meta.data[["cell_subtype"]],levels=c('CD8-C1-ZNF683-Trm','CD4-C2-IL7R-Tcm','CD8-C3-GZMK-Tem','CD4-C4-FOXP3-Treg',
                                                                                       'NK-C5-NKG7','CD8-C6-CXCL13-Tex','CD8-C8-GPR183-Tcm','CD4-C10-CXCL13-Th1','DP-C11-HIST1H1B'))
Idents(CC.NKT)<-CC.NKT@meta.data[["cell_subtype"]]
DimPlot(CC.NKT, reduction = "tsne", label = TRUE,raster=FALSE)
saveRDS(CC.NKT, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.NKT.anno.rds"))

####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

CC.epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.epi.anno.rds"))
Idents(CC.epi)<-CC.epi@meta.data[["cell_subtype"]]
CC.epi<-JoinLayers(CC.epi)
CC.epi_markers <- FindAllMarkers(CC.epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
CC.epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Cancer.Cell.Int_2025_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

CC.NKT <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.NKT.anno.rds"))
Idents(CC.NKT)<-CC.NKT@meta.data[["cell_subtype"]]
CC.NKT<-JoinLayers(CC.NKT)
CC.NKT_markers <- FindAllMarkers(CC.NKT, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
CC.NKT_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Cancer.Cell.Int_2025_NKT_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

CC.myeloid <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.myeloid.anno.rds"))
Idents(CC.myeloid)<-CC.myeloid@meta.data[["cell_subtype"]]
CC.myeloid<-JoinLayers(CC.myeloid)
CC.myeloid_markers <- FindAllMarkers(CC.myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
CC.myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Cancer.Cell.Int_2025_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

CC.fib <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.fib.anno.rds"))
Idents(CC.fib)<-CC.fib@meta.data[["cell_subtype"]]
CC.fib<-JoinLayers(CC.fib)
CC.fib_markers <- FindAllMarkers(CC.fib, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
CC.fib_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Cancer.Cell.Int_2025_fib_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

CC.ECs <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Cancer.Cell.Int_2025/Li.ECs.anno.rds"))
Idents(CC.ECs)<-CC.ECs@meta.data[["cell_subtype"]]
CC.ECs<-JoinLayers(CC.ECs)
CC.ECs_markers <- FindAllMarkers(CC.ECs, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
CC.ECs_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Cancer.Cell.Int_2025_ECs_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 04: Li et al., Communications Biology (2022), E-MTAB-11948
# Original script: single-cell/Li_Commun.Biol_2022/E-MTAB-11948.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/matrix"),full.names = T)
sceList <- lapply(mat,function(x){
  counts = read.csv(x)
  counts = subset(counts,X!='0')
  colnames(counts)<-counts[1,]
  counts<-counts[-1,]
  rownames(counts)<-counts[,1]
  counts<-counts[,-1]
  sce = CreateSeuratObject(counts = counts,
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})


samples<-c('Sample1T','Sample2T','Sample3T','Sample1N','Sample2N','Sample3N')
type<-c('Tumor','Tumor','Tumor','Normal','Normal','Normal')
orig <- NULL
group <- NULL
sce.doubletFinder.list<-list()
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  group = rep(type[i],ncol(sceList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["group"]] = group
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA >=200 & nCount_RNA >= 200 & percent.mt <= 10)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
  sce.doubletFinder.list[[i]]<-sce
  print(paste0(samples[i]," Finished!"))
}

Hua.sce<-sce.doubletFinder.list
for(i in 1:length(Hua.sce)){
  Hua.sce[[i]][["pANN"]]<-Hua.sce[[i]]@meta.data %>% select(contains('pANN'))
  Hua.sce[[i]][["DF.classifications"]]<-Hua.sce[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

Hua.sce <- merge(Hua.sce[[1]],
                 y = Hua.sce[-1],
                 add.cell.ids = samples)
c <- grep("pANN_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
c <- grep("DF.classifications_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
DimPlot(Hua.sce, reduction = 'pca',group.by = 'DF.classifications', label = TRUE,raster=FALSE)
Hua.sce <- subset(Hua.sce, DF.classifications== "Singlet")

saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.sce.rds"))
Hua.sce <- NormalizeData(Hua.sce, normalization.method = "LogNormalize")
Hua.sce <- FindVariableFeatures(Hua.sce, selection.method = "vst", nfeatures = 2000)
Hua.sce <- ScaleData(Hua.sce, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.sce))
#Hua.sce <- SCTransform(Hua.sce, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.sce <- RunPCA(Hua.sce, features = VariableFeatures(object = Hua.sce))
Hua.sce <- RunHarmony(Hua.sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.sce, ndims = 50, reduction = "harmony")
Hua.sce <- FindNeighbors(Hua.sce, dims = 1:30,reduction="harmony")
Hua.sce <- FindClusters(object = Hua.sce,resolution = seq(0.1,1,by=0.1))
clustree(Hua.sce@meta.data, prefix = "RNA_snn_res.")
Hua.sce <- FindClusters(Hua.sce, resolution = 0.7)
Hua.sce <- RunUMAP(Hua.sce,reduction="harmony", dims = 1:30)
Hua.sce <- RunTSNE(Hua.sce,reduction="harmony", dims = 1:30)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
DotPlot(Hua.sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))
DimPlot(Hua.sce, reduction = "tsne", label = TRUE,raster=FALSE)

new.cluster.ids <- c('Epithelial cells','Fibroblasts','Endothelial cells','Epithelial cells','Fibroblasts','T cells',
                     'Epithelial cells','Smooth muscle cells','Fibroblasts','Smooth muscle cells','Fibroblasts',
                     'Endothelial cells','Neutrophils','Endothelial cells','Macrophages','Epithelial cells','T cells',
                     'Endothelial cells','Epithelial cells','Mast cells','Endothelial cells','B cells','Fibroblasts')
names(new.cluster.ids) <- levels(Hua.sce)
Hua.sce <- RenameIdents(Hua.sce, new.cluster.ids)
Hua.sce@meta.data[["Cell_type"]]<-as.character(Hua.sce@active.ident)
Hua.sce@meta.data[["Cell_type"]]<-factor(Hua.sce@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblasts','Smooth muscle cells','Endothelial cells',
                                                  'T cells','Neutrophils','Macrophages','B cells','Mast cells'))
Idents(Hua.sce)<-Hua.sce@meta.data[["Cell_type"]]
DimPlot(Hua.sce, reduction = "umap",raster=FALSE)

saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.sce.anno.rds"))


# Epithelial-cell subclustering
Hua.epi<-subset(Hua.sce,Cell_type =='Epithelial cells')
#Hua.epi <- NormalizeData(Hua.epi, normalization.method = "LogNormalize")
#Hua.epi <- FindVariableFeatures(Hua.epi, selection.method = "vst", nfeatures = 3000)
#Hua.epi <- ScaleData(Hua.epi, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.epi))
#Hua.epi <- SCTransform(Hua.epi, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.epi <- RunPCA(Hua.epi, features = VariableFeatures(object = Hua.epi))
#Hua.epi <- RunHarmony(Hua.epi,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
#DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.epi, ndims = 50, reduction = "pca")
Hua.epi <- FindNeighbors(Hua.epi, dims = 1:30,reduction="pca")
Hua.epi <- FindClusters(object = Hua.epi,resolution = seq(0.1,1,by=0.1))
clustree(Hua.epi@meta.data, prefix = "RNA_snn_res.")
Hua.epi <- FindClusters(Hua.epi, resolution = 0.5)
Hua.epi <- RunUMAP(Hua.epi,reduction="pca", dims = 1:30)
Hua.epi <- RunTSNE(Hua.epi,reduction="pca", dims = 1:30)
DimPlot(Hua.epi, reduction = "umap",raster=FALSE,label = T)

marker<-c('MMP1','SPRR1B','KRT16','CSTA','S100A9',
          'CD74','IL32',
          'CCDC80','IER5','MAFB',
          'UBE2C','TOP2A','ANLN',
          'CLU','SCGB3A1','MUC5B')
DotPlot(Hua.epi, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('C2','C1','C1','C1','C5','C4','C2',
                     'C6','C2','C7','C3','C6','C4','C5')
names(new.cluster.ids) <- levels(Hua.epi)
Hua.epi <- RenameIdents(Hua.epi, new.cluster.ids)
Hua.epi@meta.data[["Cell_subtype"]]<-as.character(Hua.epi@active.ident)
Hua.epi@meta.data[["Cell_subtype"]]<-factor(Hua.epi@meta.data[["Cell_subtype"]],
                                         levels=c('C1','C2','C3','C4','C5','C6','C7'))
Idents(Hua.epi)<-Hua.epi@meta.data[["Cell_subtype"]]
DimPlot(Hua.epi, reduction = "tsne",raster=FALSE)

saveRDS(Hua.epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.epi.anno.rds"))


# T-cell subclustering
Hua.T<-subset(Hua.sce,Cell_type =='T cells')
#Hua.T <- NormalizeData(Hua.T, normalization.method = "LogNormalize")
#Hua.T <- FindVariableFeatures(Hua.T, selection.method = "vst", nfeatures = 3000)
#Hua.T <- ScaleData(Hua.T, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.T))
#Hua.T <- SCTransform(Hua.T, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.T <- RunPCA(Hua.T, features = VariableFeatures(object = Hua.T))
Hua.T <- RunHarmony(Hua.T,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.T, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.T, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.T, ndims = 50, reduction = "harmony")
Hua.T <- FindNeighbors(Hua.T, dims = 1:30,reduction="harmony")
Hua.T <- FindClusters(object = Hua.T,resolution = seq(0.1,1,by=0.1))
clustree(Hua.T@meta.data, prefix = "RNA_snn_res.")
Hua.T <- FindClusters(Hua.T, resolution = 1)
Hua.T <- RunUMAP(Hua.T,reduction="harmony", dims = 1:30)
Hua.T <- RunTSNE(Hua.T,reduction="harmony", dims = 1:30)
DimPlot(Hua.T, reduction = "umap",raster=FALSE,label = T)

marker<-c('GZMK','CXCR4','CST7',    #CD8 Tem
          'XCL1','CAPG','NR4A1',    #Trm
          'CD160','KIR2DL4','KLRC2', #IEL
          'HAVCR2','CXCL13','PDCD1', #Tex
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL',  #naive T
          'FOXP3','CTLA4','TNFRSF18', #CD4 Treg
          'CCL5',   #Tems
          'IL17A','CCR6','CTSH',  #Th17
          'IFNG', #Th1-like
          'CD4','CD8A','CD8B','NKG7'
)
marker<-c('CD8A','NKG7','IL7R','TNFRSF4')
marker<-c('GZMK','CXCR4','CXCR2','CX3CR1',
          'PDCD1','TIGIT','CTLA4','HAVCR2','LAG3','CD274','HLA-DPA1','HLA-DRA','HLA-DRB1',
          'MKI67','TOP2A','CCNB1')
DotPlot(Hua.T, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('CXCR4+CD8','PDCD1+CD8','IL7R+Tm','IL7R+Tm','NK','IL7R+Tm','MKI67+CD8','Unknown','TNFRSF4+Treg','NK')
names(new.cluster.ids) <- levels(Hua.T)
Hua.T <- RenameIdents(Hua.T, new.cluster.ids)
Hua.T@meta.data[["Cell_subtype"]]<-as.character(Hua.T@active.ident)
Hua.T@meta.data[["Cell_subtype"]]<-factor(Hua.T@meta.data[["Cell_subtype"]],
                                            levels=c('CXCR4+CD8','PDCD1+CD8','MKI67+CD8','NK','IL7R+Tm','TNFRSF4+Treg','Unknown'))
Idents(Hua.T)<-Hua.T@meta.data[["Cell_subtype"]]
DimPlot(Hua.T, reduction = "tsne",raster=FALSE)

saveRDS(Hua.T, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.T.anno.rds"))


# Fibroblast and smooth-muscle-cell subclustering
Hua.fib<-subset(Hua.sce,Cell_type %in% c('Fibroblasts','Smooth muscle cells'))
Hua.fib <- RunPCA(Hua.fib, features = VariableFeatures(object = Hua.fib))
Hua.fib <- RunHarmony(Hua.fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.fib, ndims = 50, reduction = "harmony")
Hua.fib <- FindNeighbors(Hua.fib, dims = 1:30,reduction="harmony")
Hua.fib <- FindClusters(object = Hua.fib,resolution = seq(0.1,1,by=0.1))
clustree(Hua.fib@meta.data, prefix = "RNA_snn_res.")
Hua.fib <- FindClusters(object = Hua.fib,resolution = 0.5)
Hua.fib <- RunUMAP(Hua.fib,reduction="harmony", dims = 1:30)
Hua.fib <- RunTSNE(Hua.fib,reduction="harmony", dims = 1:30)
DimPlot(Hua.fib, reduction = "tsne",raster=FALSE,label = T)

VlnPlot(Hua.fib,features = c('DCN','COL1A2','ACTA2'))
marker<-c('DCN','COL1A2','ACTA2',
          'IL6','IL8','CXCL1','CXCL2','CCL2','CXCL12')
DotPlot(Hua.fib, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

          
saveRDS(Hua.fib, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.fib.rds"))

####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Hua.epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.epi.anno.rds"))
Idents(Hua.epi)<-Hua.epi@meta.data[["Cell_subtype"]]
Hua.epi<-JoinLayers(Hua.epi)
Hua.epi_markers <- FindAllMarkers(Hua.epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Commun.Biol_2022_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.T.anno.rds"))
Idents(Hua.T)<-Hua.T@meta.data[["Cell_subtype"]]
Hua.T<-JoinLayers(Hua.T)
Hua.T_markers <- FindAllMarkers(Hua.T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Commun.Biol_2022_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.fib <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.fib.rds"))
Idents(Hua.fib)<-Hua.fib@meta.data[["seurat_clusters"]]
Hua.fib<-JoinLayers(Hua.fib)
Hua.fib_markers <- FindAllMarkers(Hua.fib, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.fib_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Commun.Biol_2022_fib_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 05: Li et al., Frontiers in Immunology (2022), E-MTAB-12305
# Original script: single-cell/Li_Front.Immunol_2022/E-MTAB-12305.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/matrix"),full.names = T)
sceList <- lapply(mat,function(x){
  counts = read.csv(x)
  counts = subset(counts,X!='0')
  colnames(counts)<-counts[1,]
  counts<-counts[-1,]
  rownames(counts)<-counts[,1]
  counts<-counts[,-1]
  sce = CreateSeuratObject(counts = counts,
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})


samples<-c('H1','H2','L1','N1','N2','N3','T1','T2','T3','T4')
type<-c('HSIL','HSIL','Metastatic','Normal','Normal','Normal','CC','CC','CC','CC')
orig <- NULL
group <- NULL
sce.doubletFinder.list<-list()
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  group = rep(type[i],ncol(sceList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["group"]] = group
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nCount_RNA >= 200 & percent.mt <= 10)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:30,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:30)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
  sce.doubletFinder.list[[i]]<-sce
  print(paste0(samples[i]," Finished!"))
}

Hua.sce<-sce.doubletFinder.list
for(i in 1:length(Hua.sce)){
  Hua.sce[[i]][["pANN"]]<-Hua.sce[[i]]@meta.data %>% select(contains('pANN'))
  Hua.sce[[i]][["DF.classifications"]]<-Hua.sce[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

Hua.sce <- merge(Hua.sce[[1]],
                 y = Hua.sce[-1],
                 add.cell.ids = samples)
c <- grep("pANN_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
c <- grep("DF.classifications_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
DimPlot(Hua.sce, reduction = 'pca',group.by = 'DF.classifications', label = TRUE,raster=FALSE)
Hua.sce <- subset(Hua.sce, DF.classifications== "Singlet")

# TODO: This original path writes into the Li_Commun.Biol_2022 dataset directory; verify before running.
saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.sce.rds"))

Hua.sce <- NormalizeData(Hua.sce, normalization.method = "LogNormalize")
Hua.sce <- FindVariableFeatures(Hua.sce, selection.method = "vst", nfeatures = 2000)
Hua.sce <- ScaleData(Hua.sce, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.sce))
#Hua.sce <- SCTransform(Hua.sce, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.sce <- RunPCA(Hua.sce, features = VariableFeatures(object = Hua.sce))
Hua.sce <- RunHarmony(Hua.sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.sce, ndims = 50, reduction = "harmony")
Hua.sce <- FindNeighbors(Hua.sce, dims = 1:30,reduction="harmony")
Hua.sce <- FindClusters(object = Hua.sce,resolution = seq(0.1,1,by=0.1))
clustree(Hua.sce@meta.data, prefix = "RNA_snn_res.")
Hua.sce <- FindClusters(Hua.sce, resolution = 0.7)
Hua.sce <- RunUMAP(Hua.sce,reduction="harmony", dims = 1:30)
Hua.sce <- RunTSNE(Hua.sce,reduction="harmony", dims = 1:30)
DimPlot(Hua.sce, reduction = "tsne", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
DotPlot(Hua.sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))


new.cluster.ids <- c('Fibroblasts','Epithelial cells','T/NK cells','T/NK cells','Endothelial cells','Fibroblasts',
                     'Epithelial cells','Smooth muscle cells','Epithelial cells','Macrophages','Neutrophils','Epithelial cells',
                     'Endothelial cells','Fibroblasts','Epithelial cells','Fibroblasts','Mast cells','Epithelial cells',
                     'Plasma cells','T/NK cells','B cells','Epithelial cells','Epithelial cells','Endothelial cells',
                     'Macrophages','Epithelial cells','Macrophages','Fibroblasts','Endothelial cells')
names(new.cluster.ids) <- levels(Hua.sce)
Hua.sce <- RenameIdents(Hua.sce, new.cluster.ids)
Hua.sce@meta.data[["Cell_type"]]<-as.character(Hua.sce@active.ident)
Hua.sce@meta.data[["Cell_type"]]<-factor(Hua.sce@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblasts','Smooth muscle cells','Endothelial cells',
                                                  'T/NK cells','Neutrophils','Macrophages','B cells','Plasma cells','Mast cells'))
Idents(Hua.sce)<-Hua.sce@meta.data[["Cell_type"]]
DimPlot(Hua.sce, reduction = "umap",raster=FALSE)

saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/E-MTAB-12305.sce.anno.rds"))

# Epithelial-cell subclustering
Hua.epi<-subset(Hua.sce,Cell_type =='Epithelial cells')
Hua.epi<-Hua.epi[rownames(Hua.sce),]
#Hua.epi <- NormalizeData(Hua.epi, normalization.method = "LogNormalize")
#Hua.epi <- FindVariableFeatures(Hua.epi, selection.method = "vst", nfeatures = 3000)
#Hua.epi <- ScaleData(Hua.epi, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.epi))
#Hua.epi <- SCTransform(Hua.epi, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.epi <- RunPCA(Hua.epi, features = VariableFeatures(object = Hua.epi))
#Hua.epi <- RunHarmony(Hua.epi,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
#DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.epi, ndims = 50, reduction = "pca")
Hua.epi <- FindNeighbors(Hua.epi, dims = 1:30,reduction="pca")
Hua.epi <- FindClusters(object = Hua.epi,resolution = seq(0.1,1,by=0.1))
clustree(Hua.epi@meta.data, prefix = "RNA_snn_res.")
Hua.epi <- FindClusters(Hua.epi, resolution = 0.5)
Hua.epi <- RunUMAP(Hua.epi,reduction="pca", dims = 1:30)
Hua.epi <- RunTSNE(Hua.epi,reduction="pca", dims = 1:30)
DimPlot(Hua.epi, reduction = "umap",raster=FALSE,label = T)


# T-cell subclustering
Hua.NKT<-subset(Hua.sce,Cell_type =='T/NK cells')
Hua.NKT <- RunPCA(Hua.NKT, features = VariableFeatures(object = Hua.NKT))
Hua.NKT <- RunHarmony(Hua.NKT,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.NKT, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.NKT, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.NKT, ndims = 50, reduction = "harmony")
Hua.NKT <- FindNeighbors(Hua.NKT, dims = 1:30,reduction="harmony")
Hua.NKT <- FindClusters(object = Hua.NKT,resolution = seq(0.1,1,by=0.1))
clustree(Hua.NKT@meta.data, prefix = "RNA_snn_res.")
Hua.NKT <- FindClusters(Hua.NKT, resolution = 1)
Hua.NKT <- RunUMAP(Hua.NKT,reduction="harmony", dims = 1:30)
Hua.NKT <- RunTSNE(Hua.NKT,reduction="harmony", dims = 1:30)
DimPlot(Hua.NKT, reduction = "umap",raster=FALSE,label = T)

marker<-c('GZMK','CXCR4','CST7','GZMH',    #CD8 Tem
          'XCL1','CAPG','NR4A1',    #Trm
          'CD160','KIR2DL4','KLRC2', #IEL
          'HAVCR2','CXCL13','PDCD1', #Tex
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL',  #naive T
          'FOXP3','CTLA4','TNFRSF18', #CD4 Treg
          'CCL5',   #Tems
          'IL17A','CCR6','CTSH',  #Th17
          'IFNG', #Th1-like
          'CD4','CD8A','CD8B','NKG7','TOP2A','MKI67'
)
marker<-c('GZMK','CXCR4','CST7','GZMH',    #CD8 Tex
          'CD69','EGR1','RUNX3','NR4A1',   #Trm
          'HAVCR2','LAG3','PDCD1','CTLA4', #Tex
          'HSPA1B','HSPB1','HSPA1A',
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL','TCF7',  #naive T
          'FOXP3','IL2RA','TNFRSF18', #CD4 Treg
          'CCL5',   #Tems
          'LTB','IL7R',  #Tcm
          'IL17A','CCR6','CTSH',  #Th17
          'IFNG', #Th1-like
          'GZMA','GZMB','GNLY','CD160',  #NK eff
          'S100A2','S100A8','S100A9','S100A11'
)
DotPlot(Hua.NKT, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))



new.cluster.ids <- c('Cytotoxic CD8 T cells','Naive CD4 T cells','Effector NK cells','Reg T cells','Naive CD4 T cells',
                     'Effector memory CD8 T cells','Central memory T cells','Effector memory CD8 T cells','Exhausted CD8 T cells','Effector NK cells','Central memory T cells',
                     'Cytotoxic CD8 T cells','Central memory T cells','Cycling T cells','Naive CD4 T cells','Circulating NK cells')
names(new.cluster.ids) <- levels(Hua.NKT)
Hua.NKT <- RenameIdents(Hua.NKT, new.cluster.ids)
Hua.NKT@meta.data[["Cell_subtype"]]<-as.character(Hua.NKT@active.ident)
Idents(Hua.NKT)<-Hua.NKT@meta.data[["Cell_subtype"]]
DimPlot(Hua.NKT, reduction = "tsne",raster=FALSE)
saveRDS(Hua.NKT, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/E-MTAB-12305.NKT.anno.rds"))


# Myeloid-cell subclustering
Hua.myeloid<-subset(Hua.sce,Cell_type =='Macrophages')
Hua.myeloid <- RunPCA(Hua.myeloid, features = VariableFeatures(object = Hua.myeloid))
Hua.myeloid <- RunHarmony(Hua.myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.myeloid, ndims = 50, reduction = "harmony")
Hua.myeloid <- FindNeighbors(Hua.myeloid, dims = 1:30,reduction="harmony")
Hua.myeloid <- FindClusters(object = Hua.myeloid,resolution = seq(0.1,1,by=0.1))
clustree(Hua.myeloid@meta.data, prefix = "RNA_snn_res.")
Hua.myeloid <- FindClusters(Hua.myeloid, resolution = 0.8)
Hua.myeloid <- RunUMAP(Hua.myeloid,reduction="harmony", dims = 1:30)
Hua.myeloid <- RunTSNE(Hua.myeloid,reduction="harmony", dims = 1:30)
DimPlot(Hua.myeloid, reduction = "tsne",raster=FALSE,label = T)

marker<-c('CD163','MRC1','CD14',
          'CD1C','LILRA4','CLEC10A')
marker<-c('C1QA','MARCO','APOE','CXCL10',
          'THSB1','FCN1','VCAN','S100A12',
          'GNLY','IFNG','NKG7','GZMA',
          'CSTA','CSTD','S100A8','S100A2','S100A6',
          'HSPH1','KLF2','KLF4',
          'STMN1','TOP2A','MKI67',
          'TMSB4X','TMSB10','CSTB'
          )
marker<-c('XCR1','CLEC9A','BATF3',
          'CD1C','FCER1A',
          'CLEC4C','LILRA4','LILRA5',
          'GZMB','CXCR3','CXCR4','BCL11A','RUNX2')
DotPlot(Hua.myeloid, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('C1-Ma-C1QA','C4-Ma-GZMA','C2-Ma-THBS1','C3-DC2','C2-Ma-THBS1','C1-Ma-C1QA','C8-Ma-SPP1',
                     'C6-Ma-MAFB','C7-Ma-MKI67','C3-DC2','C4-Ma-GZMA','C2-Ma-THBS1','C5-Ma-S100A8','C9-DC1','C6-Ma-MAFB')
names(new.cluster.ids) <- levels(Hua.myeloid)
Hua.myeloid <- RenameIdents(Hua.myeloid, new.cluster.ids)
Hua.myeloid@meta.data[["Cell_subtype"]]<-as.character(Hua.myeloid@active.ident)
Hua.myeloid@meta.data[["Cell_subtype"]]<-factor(Hua.myeloid@meta.data[["Cell_subtype"]],
                                                levels=c('C1-Ma-C1QA','C2-Ma-THBS1','C4-Ma-GZMA','C5-Ma-S100A8','C6-Ma-MAFB',
                                                         'C7-Ma-MKI67','C8-Ma-SPP1','C3-DC2','C9-DC1'))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
DimPlot(Hua.myeloid, reduction = "umap",raster=FALSE)
saveRDS(Hua.myeloid, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/E-MTAB-12305.myeloid.anno.rds"))


# Fibroblast subclustering
Hua.fib<-subset(Hua.sce,Cell_type == 'Fibroblasts')
Hua.fib <- RunPCA(Hua.fib, features = VariableFeatures(object = Hua.fib))
Hua.fib <- RunHarmony(Hua.fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.fib, ndims = 50, reduction = "harmony")
Hua.fib <- FindNeighbors(Hua.fib, dims = 1:30,reduction="harmony")
Hua.fib <- FindClusters(object = Hua.fib,resolution = seq(0.1,1,by=0.1))
clustree(Hua.fib@meta.data, prefix = "RNA_snn_res.")
Hua.fib <- RunUMAP(Hua.fib,reduction="harmony", dims = 1:30)
Hua.fib <- RunTSNE(Hua.fib,reduction="harmony", dims = 1:30)
DimPlot(Hua.fib, reduction = "tsne",raster=FALSE,label = T)
saveRDS(Hua.fib, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/E-MTAB-12305.fib.rds"))

####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Hua.NKT <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/E-MTAB-12305.NKT.anno.rds"))
Idents(Hua.NKT)<-Hua.NKT@meta.data[["Cell_subtype"]]
Hua.NKT<-JoinLayers(Hua.NKT)
Hua.NKT_markers <- FindAllMarkers(Hua.NKT, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.NKT_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Front.Immunol_2022_NKT_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.myeloid <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Front.Immunol_2022/E-MTAB-12305.myeloid.anno.rds"))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
Hua.myeloid<-JoinLayers(Hua.myeloid)
Hua.myeloid_markers <- FindAllMarkers(Hua.myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Front.Immunol_2022_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 06: Li et al., Journal of Medical Virology (2023), S-BSST1035
# Original script: single-cell/Li_J.Med.Virol_2023/S-BSST1035.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/matrix"),full.names = T)
sceList <- lapply(mat,function(x){
  counts = read.csv(x)
  counts = subset(counts,X!='0')
  colnames(counts)<-counts[1,]
  counts<-counts[-1,]
  rownames(counts)<-counts[,1]
  counts<-counts[,-1]
  sce = CreateSeuratObject(counts = counts,
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})


samples<-c('AD1','AD2','AD3','SCC1','SCC2','SCC3')
type<-c('AD','AD','AD','SCC','SCC','SCC')
orig <- NULL
group <- NULL
sce.doubletFinder.list<-list()
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  group = rep(type[i],ncol(sceList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["group"]] = group
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce <- subset(sce, subset = nFeature_RNA >=200 & nCount_RNA >= 800 & percent.mt <= 20)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
  sce.doubletFinder.list[[i]]<-sce
  print(paste0(samples[i]," Finished!"))
}

Hua.sce<-sce.doubletFinder.list
for(i in 1:length(Hua.sce)){
  Hua.sce[[i]][["pANN"]]<-Hua.sce[[i]]@meta.data %>% select(contains('pANN'))
  Hua.sce[[i]][["DF.classifications"]]<-Hua.sce[[i]]@meta.data %>% select(contains('DF.classifications'))
  #CC.sce[[i]]@meta.data %>% select(contains('DF.classifications_0.25'))<-NULL
  #CC.sce[[i]]@meta.data %>% select(contains('pANN_0.25'))<-NULL
}

Hua.sce <- merge(Hua.sce[[1]],
                 y = Hua.sce[-1],
                 add.cell.ids = samples)
c <- grep("pANN_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
c <- grep("DF.classifications_0.25",colnames(Hua.sce@meta.data))
Hua.sce@meta.data <- Hua.sce@meta.data[,-c]
DimPlot(Hua.sce, reduction = 'pca',group.by = 'DF.classifications', label = TRUE,raster=FALSE)
Hua.sce <- subset(Hua.sce, DF.classifications== "Singlet")

# TODO: This original path writes into the Li_Commun.Biol_2022 dataset directory; verify before running.
saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Commun.Biol_2022/E-MTAB-11948.sce.rds"))
Hua.sce <- NormalizeData(Hua.sce, normalization.method = "LogNormalize")
Hua.sce <- FindVariableFeatures(Hua.sce, selection.method = "vst", nfeatures = 2000)
Hua.sce <- ScaleData(Hua.sce, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.sce))
#Hua.sce <- SCTransform(Hua.sce, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.sce <- RunPCA(Hua.sce, features = VariableFeatures(object = Hua.sce))
Hua.sce <- RunHarmony(Hua.sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.sce, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.sce, ndims = 50, reduction = "harmony")
Hua.sce <- FindNeighbors(Hua.sce, dims = 1:20,reduction="harmony")
Hua.sce <- FindClusters(object = Hua.sce,resolution = seq(0.1,1,by=0.1))
clustree(Hua.sce@meta.data, prefix = "RNA_snn_res.")
Hua.sce <- FindClusters(Hua.sce, resolution = 0.8)
Hua.sce <- RunUMAP(Hua.sce,reduction="harmony", dims = 1:20)
Hua.sce <- RunTSNE(Hua.sce,reduction="harmony", dims = 1:20)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
DotPlot(Hua.sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))
DimPlot(Hua.sce, reduction = "tsne", label = TRUE,raster=FALSE)

new.cluster.ids <- c('Epithelial cells','Epithelial cells','Epithelial cells','Epithelial cells','NK/T cells','Epithelial cells',
                     'NK/T cells','Epithelial cells','Myeloid cells','Epithelial cells','Epithelial cells','Fibroblasts',
                     'Fibroblasts','B cells','NK/T cells','Neutrophils','Myeloid cells','NK/T cells','Plasma cells',
                     'Epithelial cells','Mast cells','Epithelial cells','Endothelial cells','Epithelial cells')
names(new.cluster.ids) <- levels(Hua.sce)
Hua.sce <- RenameIdents(Hua.sce, new.cluster.ids)
Hua.sce@meta.data[["Cell_type"]]<-as.character(Hua.sce@active.ident)
Hua.sce@meta.data[["Cell_type"]]<-factor(Hua.sce@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblasts','Endothelial cells','NK/T cells',
                                                  'Neutrophils','Myeloid cells','B cells','Plasma cells','Mast cells'))
Idents(Hua.sce)<-Hua.sce@meta.data[["Cell_type"]]
DimPlot(Hua.sce, reduction = "umap",raster=FALSE)

saveRDS(Hua.sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.sce.anno.rds"))


# Epithelial-cell subclustering
Hua.epi<-subset(Hua.sce,Cell_type =='Epithelial cells')
#Hua.epi <- NormalizeData(Hua.epi, normalization.method = "LogNormalize")
#Hua.epi <- FindVariableFeatures(Hua.epi, selection.method = "vst", nfeatures = 3000)
#Hua.epi <- ScaleData(Hua.epi, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.epi))
#Hua.epi <- SCTransform(Hua.epi, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.epi <- RunPCA(Hua.epi, features = VariableFeatures(object = Hua.epi))
#Hua.epi <- RunHarmony(Hua.epi,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
#DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.epi, ndims = 50, reduction = "pca")
Hua.epi <- FindNeighbors(Hua.epi, dims = 1:30,reduction="pca")
Hua.epi <- FindClusters(object = Hua.epi,resolution = seq(0.1,1,by=0.1))
clustree(Hua.epi@meta.data, prefix = "RNA_snn_res.")
Hua.epi <- FindClusters(Hua.epi, resolution = 0.3)
Hua.epi <- RunUMAP(Hua.epi,reduction="pca", dims = 1:30)
Hua.epi <- RunTSNE(Hua.epi,reduction="pca", dims = 1:30)
DimPlot(Hua.epi, reduction = "umap",raster=FALSE,label = T)

marker<-c('EPCAM', 'MUC5B', 'KRT8', 'KRT18', 'AQP3',
          'KRT6A', 'TP63', 'CDKN2A', 'KRT13', 'DSG3')
marker<-c('KRT14','FABP4','MKI67','HLA-DRB5','MUC5B','VIM','COL18A1','GSTM3','ATP2A3','IL19')
DotPlot(Hua.epi, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('C0','C3','C5','C1','C4','C6','C2',
                     'C9','C8','C3','C5','C0','C1','C7')
names(new.cluster.ids) <- levels(Hua.epi)
Hua.epi <- RenameIdents(Hua.epi, new.cluster.ids)
Hua.epi@meta.data[["Cell_subtype"]]<-as.character(Hua.epi@active.ident)
Hua.epi@meta.data[["Cell_subtype"]]<-factor(Hua.epi@meta.data[["Cell_subtype"]],
                                            levels=c('C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'))
Idents(Hua.epi)<-Hua.epi@meta.data[["Cell_subtype"]]
DimPlot(Hua.epi, reduction = "umap",raster=FALSE)

saveRDS(Hua.epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.epi.anno.rds"))


# T-cell subclustering
Hua.T<-subset(Hua.sce,Cell_type =='NK/T cells')
#Hua.T <- NormalizeData(Hua.T, normalization.method = "LogNormalize")
#Hua.T <- FindVariableFeatures(Hua.T, selection.method = "vst", nfeatures = 3000)
#Hua.T <- ScaleData(Hua.T, vars.to.regress = "percent.mt",features = VariableFeatures(object = Hua.T))
#Hua.T <- SCTransform(Hua.T, vars.to.regress = "percent.mt", verbose = FALSE)
Hua.T <- RunPCA(Hua.T, features = VariableFeatures(object = Hua.T))
Hua.T <- RunHarmony(Hua.T,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.T, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.T, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.T, ndims = 50, reduction = "harmony")
Hua.T <- FindNeighbors(Hua.T, dims = 1:30,reduction="harmony")
Hua.T <- FindClusters(object = Hua.T,resolution = seq(0.1,1,by=0.1))
clustree(Hua.T@meta.data, prefix = "RNA_snn_res.")
Hua.T <- FindClusters(Hua.T, resolution = 1)
Hua.T <- RunUMAP(Hua.T,reduction="harmony", dims = 1:30)
Hua.T <- RunTSNE(Hua.T,reduction="harmony", dims = 1:30)
DimPlot(Hua.T, reduction = "umap",raster=FALSE,label = T)

marker<-c('GZMK','CXCR4','CST7','GZMH',    #CD8 Tem
          'XCL1','CAPG','NR4A1',    #Trm
          'CD160','KIR2DL4','KLRC2', #IEL
          'HAVCR2','CXCL13','PDCD1', #Tex
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL',  #naive T
          'FOXP3','CTLA4','TNFRSF18', #CD4 Treg
          'CCL5',   #Tems
          'IL17A','CCR6','CTSH',  #Th17
          'IFNG', #Th1-like
          'CD4','CD8A','CD8B','NKG7','TOP2A','MKI67'
)
marker<-c('GZMK','CXCR4','CST7','GZMH',    #CD8 Tex
          'CD69','EGR1','RUNX3','NR4A1',   #Trm
          'HAVCR2','LAG3','PDCD1','CTLA4', #Tex
          'HSPA1B','HSPB1','HSPA1A',
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL','TCF7',  #naive T
          'FOXP3','IL2RA','TNFRSF18', #CD4 Treg
          'CCL5',   #Tems
          'LTB','IL7R',  #Tcm
          'IL17A','CCR6','CTSH',  #Th17
          'IFNG', #Th1-like
          'GZMA','GZMB','GNLY','CD160',  #NK eff
          'S100A2','S100A8','S100A9','S100A11'
)
marker<-c('CD8A','NKG7','IL7R','TNFRSF4')
marker<-c('GZMK','CXCR4','CXCR2','CX3CR1',
          'PDCD1','TIGIT','CTLA4','HAVCR2','LAG3','CD274','HLA-DPA1','HLA-DRA','HLA-DRB1',
          'MKI67','TOP2A','CCNB1')
DotPlot(Hua.T, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('Tissue resident CD8 T cells','Central memory CD4 T cells','CD160+ NK cells','Naive CD4 T cells','Tissue resident CD8 T cells','Treg CD4 T cells','Central memory CD4 T cells',
                     'Effector memory CD8 T cells','Proliferative cells','Exhausted CD8 T cells','Central memory CD4 T cells','Central memory CD4 T cells','Tissue resident NK cells')
names(new.cluster.ids) <- levels(Hua.T)
Hua.T <- RenameIdents(Hua.T, new.cluster.ids)
Hua.T@meta.data[["Cell_subtype"]]<-as.character(Hua.T@active.ident)
Idents(Hua.T)<-Hua.T@meta.data[["Cell_subtype"]]
DimPlot(Hua.T, reduction = "tsne",raster=FALSE)

saveRDS(Hua.T, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.T.anno.rds"))


# Fibroblast subclustering
Hua.fib<-subset(Hua.sce,Cell_type == 'Fibroblasts')
Hua.fib <- RunPCA(Hua.fib, features = VariableFeatures(object = Hua.fib))
Hua.fib <- RunHarmony(Hua.fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.fib, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.fib, ndims = 50, reduction = "harmony")
Hua.fib <- FindNeighbors(Hua.fib, dims = 1:30,reduction="harmony")
Hua.fib <- FindClusters(object = Hua.fib,resolution = seq(0.1,1,by=0.1))
clustree(Hua.fib@meta.data, prefix = "RNA_snn_res.")
Hua.fib <- FindClusters(object = Hua.fib,resolution = 1)
Hua.fib <- RunUMAP(Hua.fib,reduction="harmony", dims = 1:30)
Hua.fib <- RunTSNE(Hua.fib,reduction="harmony", dims = 1:30)
DimPlot(Hua.fib, reduction = "tsne",raster=FALSE,label = T)

saveRDS(Hua.fib, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.fib.rds"))


# Myeloid-cell subclustering
Hua.myeloid<-subset(Hua.sce,Cell_type =='Myeloid cells')
Hua.myeloid <- RunPCA(Hua.myeloid, features = VariableFeatures(object = Hua.myeloid))
Hua.myeloid <- RunHarmony(Hua.myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.myeloid, ndims = 50, reduction = "harmony")
Hua.myeloid <- FindNeighbors(Hua.myeloid, dims = 1:30,reduction="harmony")
Hua.myeloid <- FindClusters(object = Hua.myeloid,resolution = seq(0.1,1,by=0.1))
clustree(Hua.myeloid@meta.data, prefix = "RNA_snn_res.")
Hua.myeloid <- FindClusters(Hua.myeloid, resolution = 1)
Hua.myeloid <- RunUMAP(Hua.myeloid,reduction="harmony", dims = 1:30)
Hua.myeloid <- RunTSNE(Hua.myeloid,reduction="harmony", dims = 1:30)
DimPlot(Hua.myeloid, reduction = "umap",raster=FALSE,label = T)

marker<-c('CD163', 'FCGR3A', 'CD68', 
          'HLA‐A', 'HLA‐DQA2', 'CD74','FN1',
          'CCL18', 'CCL13', 'CXCL5','CXCL14',
          'LCN2','SLP1','CLDN10','AGR3','ESPN',
          'KRT16','FABP5','S100A2','S100A9','S100A8',
          'SFRP1','PDLIM3','FAT1','EN1','SOD3',
          #'STMN1','CCNB1','CCNA2','CDKN3','TOP2A','MKI67',
          'AQP5','FGFR2','SCGB3A1','TFF3','BPIFB1',
          'CD1C', 'FCER1A',  'CLEC10A',   #C3‐CD1C cDC2
          'LAMP3', 'CD80', 'CD83','CCR7','CCL19','CCL21', #LAMP3+ DCs
          'TOP2A','MKI67'
          )
DotPlot(Hua.myeloid, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('C1-Ma-FCGR3A','C1-Ma-FCGR3A','C1-Ma-FCGR3A','C7-Ma6','C7-Ma6','C3-CD1C-cDC2','C4-Ma4',
                     'C0-Ma-CD74','C8-Ma7','Cycling cells','C1-Ma-FCGR3A','C2-Ma-CD68','C9-LAMP3 DCs')
names(new.cluster.ids) <- levels(Hua.myeloid)
Hua.myeloid <- RenameIdents(Hua.myeloid, new.cluster.ids)
Hua.myeloid@meta.data[["Cell_subtype"]]<-as.character(Hua.myeloid@active.ident)
#Hua.myeloid@meta.data[["Cell_subtype"]]<-factor(Hua.myeloid@meta.data[["Cell_subtype"]],
#                                                levels=c('MAIT','Trm','Naive','Tex','IEL','Tem','TemRA'))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
DimPlot(Hua.myeloid, reduction = "umap",raster=FALSE)
saveRDS(Hua.myeloid, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.myeloid.anno.rds"))



####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Hua.T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.T.anno.rds"))
Idents(Hua.T)<-Hua.T@meta.data[["Cell_subtype"]]
Hua.T<-JoinLayers(Hua.T)
Hua.T_markers <- FindAllMarkers(Hua.T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_J.Med.Virol_2023_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.epi.anno.rds"))
Idents(Hua.epi)<-Hua.epi@meta.data[["Cell_subtype"]]
Hua.epi<-JoinLayers(Hua.epi)
Hua.epi_markers <- FindAllMarkers(Hua.epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_J.Med.Virol_2023_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.myeloid <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_J.Med.Virol_2023/S-BSST1035.myeloid.anno.rds"))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
Hua.myeloid<-JoinLayers(Hua.myeloid)
Hua.myeloid_markers <- FindAllMarkers(Hua.myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_J.Med.Virol_2023_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 07: Li et al., Molecular Therapy Nucleic Acids (2021), GSE168652
# Original script: single-cell/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652"),full.names = T)
sceList <- lapply(mat,function(x){
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samples<-list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652"))
sceList[[1]][["orig.ident"]] = samples[1]
sceList[[1]][["percent.mt"]] <- PercentageFeatureSet(sceList[[1]], pattern = "^MT-")
sceList[[2]][["orig.ident"]] = samples[2]
sceList[[2]][["percent.mt"]] <- PercentageFeatureSet(sceList[[2]], pattern = "^MT-")
sce <- merge(sceList[[1]],
             y = sceList[-1],
             add.cell.ids = samples)
VlnPlot(sce,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
sce <- subset(sce, subset = nCount_RNA >= 200 & percent.mt <= 10)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- FindNeighbors(sce, dims = 1:30,reduction="pca")
sce <- FindClusters(object = sce,resolution = seq(0.1,1,by=0.1))
clustree(sce@meta.data, prefix = "RNA_snn_res.")
sce <- FindClusters(sce, resolution = 0.5)
sce <- RunUMAP(sce,reduction="pca", dims = 1:30)
sce <- RunTSNE(sce,reduction="pca", dims = 1:30)
DimPlot(sce, reduction = "tsne", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2','SUSD2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
DotPlot(sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))
DimPlot(sce, reduction = "tsne", label = TRUE,raster=FALSE)

new.cluster.ids <- c('Epithelial cells','Smooth muscle cells','Epithelial cells','Epithelial cells','Endostromal cells','Endothelial cells','Fibroblast',
                     'Smooth muscle cells','Endostromal cells','Lymphocyte','Epithelial cells','Epithelial cells','Macrophages','Endostromal cells')
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
sce@meta.data[["Cell_type"]]<-as.character(sce@active.ident)
sce@meta.data[["Cell_type"]]<-factor(sce@meta.data[["Cell_type"]],
                                         levels=c('Epithelial cells','Fibroblast','Smooth muscle cells','Endostromal cells',
                                                  'Endothelial cells','Lymphocyte','Macrophages'))
Idents(sce)<-sce@meta.data[["Cell_type"]]
DimPlot(sce, reduction = "umap",raster=FALSE)
saveRDS(sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652.sce.anno.rds"))

sce<-sce[rownames(sceList[[2]]),]
epi<-subset(sce,Cell_type=='Epithelial cells')
epi <- NormalizeData(epi, normalization.method = "LogNormalize", scale.factor = 10000)
epi <- FindVariableFeatures(epi, selection.method = "vst", nfeatures = 2000)
epi <- ScaleData(epi, vars.to.regress = "percent.mt")
epi <- RunPCA(epi, features = VariableFeatures(object = epi))
epi <- FindNeighbors(epi, dims = 1:30,reduction="pca")
epi <- FindClusters(object = epi,resolution = seq(0.1,1,by=0.1))
clustree(epi@meta.data, prefix = "RNA_snn_res.")
epi <- FindClusters(epi, resolution = 0.1)
epi <- RunUMAP(epi,reduction="pca", dims = 1:30)
epi <- RunTSNE(epi,reduction="pca", dims = 1:30)
DimPlot(epi, reduction = "tsne", label = TRUE,raster=FALSE)

marker<-c('MKI67','CCNB1','TOP2A',
          'CEACAM6','CD55','GJB6',
          'DEFB1','KRT5','CXCL10','WARS','LY6D','CRABP2','CXCL9','PHOV','C10orf99','PYCARD',
          'SOX2','ALDHA1')
DotPlot(epi, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('C1','C3','C0','C8')
names(new.cluster.ids) <- levels(epi)
epi <- RenameIdents(epi, new.cluster.ids)
epi@meta.data[["Cell_subtype"]]<-as.character(epi@active.ident)
epi@meta.data[["Cell_subtype"]]<-factor(epi@meta.data[["Cell_subtype"]],
                                     levels=c('C0','C1','C3','C8'))
Idents(epi)<-epi@meta.data[["Cell_subtype"]]
DimPlot(epi, reduction = "umap",raster=FALSE)
saveRDS(epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652.epi.anno.rds"))



ECs<-subset(sce,Cell_type=='Endothelial cells')
ECs <- NormalizeData(ECs, normalization.method = "LogNormalize", scale.factor = 10000)
ECs <- FindVariableFeatures(ECs, selection.method = "vst", nfeatures = 2000)
ECs <- ScaleData(ECs, vars.to.regress = "percent.mt")
ECs <- RunPCA(ECs, features = VariableFeatures(object = ECs))
ECs <- FindNeighbors(ECs, dims = 1:30,reduction="pca")
ECs <- FindClusters(object = ECs,resolution = seq(0.1,1,by=0.1))
clustree(ECs@meta.data, prefix = "RNA_snn_res.")
ECs <- FindClusters(ECs, resolution = 1)
ECs <- RunUMAP(ECs,reduction="pca", dims = 1:30)
ECs <- RunTSNE(ECs,reduction="pca", dims = 1:30)
DimPlot(ECs, reduction = "tsne", label = TRUE,raster=FALSE)

marker<-c('ELN','EGR1',
          'LRG1','ICAM1',
          'IGFBP3','ARL15',
          'EMCN','FOS')
VlnPlot(ECs,features = marker)
DotPlot(ECs, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

saveRDS(ECs, file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652.ECs.rds"))


####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Li_Mol.Ther.Nucleic.Acids_2021/GSE168652.epi.anno.rds"))
Idents(epi)<-epi@meta.data[["Cell_subtype"]]
epi_markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Li_Mol.Ther.Nucleic.Acids_2021_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 08: Lin et al., EBioMedicine (2023)
# Original script: single-cell/Lin_EBioMedicine_2023/forDB.R
# =============================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(DoubletFinder)
library(clustree)
library(glmGamPoi)
library(ggsci)

####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Lin_epi <- readRDS(file.path(local_project_dir, "宫颈癌空转/ebiomedicine_scRNA_rds/癌细胞/cancer_raw_merge.rds"))
Lin_epi <- RunTSNE(Lin_epi,reduction="pca", dims = 1:20)
DimPlot(Lin_epi, reduction = "tsne",label = TRUE,raster=FALSE)
anno <- c("0"="C1", "1"="C5", "2"="C3", "3"="C6",
          "4"="C4", "5"="C2", "6"="C2", "7"="C8",
          "8"="C7", "9"="C3", "10"="C7", "11"="C4", 
          "12"="C9", "13"="C10","14"="C8","15"="C4")
Lin_epi[['cell_subtype']] = unname(anno[Lin_epi@meta.data$seurat_clusters])
Idents(Lin_epi)<-Lin_epi@meta.data[["cell_subtype"]]
Lin_epi_markers <- FindAllMarkers(Lin_epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Lin_epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Lin_EBioMedicine_2023_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

saveRDS(Lin_epi,file.path(local_project_dir, "宫颈癌数据库/单细胞/Lin_EBioMedicine_2023/EBioMedicine.epi.rds"))

Lin_T <- readRDS(file.path(local_project_dir, "宫颈癌空转/ebiomedicine_scRNA_rds/T细胞/t_scRNA.rds"))
DimPlot(Lin_T, reduction = "tsne",label = TRUE,raster=FALSE)
Lin_T[['cell_subtype']] <- Lin_T@active.ident
Lin_T_markers <- FindAllMarkers(Lin_T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Lin_T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Lin_EBioMedicine_2023_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)


saveRDS(Lin_T,file.path(local_project_dir, "宫颈癌数据库/单细胞/Lin_EBioMedicine_2023/EBioMedicine.T.rds"))


Lin_myeloid <- readRDS(file.path(local_project_dir, "宫颈癌空转/ebiomedicine_scRNA_rds/髓系细胞/Result.RDS"))
DimPlot(Lin_myeloid, reduction = "tsne",label = TRUE,raster=FALSE)
anno <- c("0"="Macro-c1", "1"="Macro-c3", "2"="Macro-c2", "3"="Macro-c4",
          "4"="DC-c1", "5"="Macro-c5", "6"="Other", "7"="DC-c2",
          "8"="Macro-c6", "9"="Other", "10"="Cycling", "11"="Other", 
          "12"="Cycling", "13"="Macro-c7", "14"="Other", "15"="DC-c3",
          "16"="DC-c4","17"="DC-c5")
Lin_myeloid[['cell_subtype']] = unname(anno[Lin_myeloid@meta.data$seurat_clusters])
Idents(Lin_myeloid)<-factor(Lin_myeloid@meta.data[["cell_subtype"]],levels=unique(Lin_myeloid@meta.data[["cell_subtype"]]))
Lin_myeloid_markers <- FindAllMarkers(Lin_myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Lin_myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Lin_EBioMedicine_2023_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

saveRDS(Lin_myeloid,file.path(local_project_dir, "宫颈癌数据库/单细胞/Lin_EBioMedicine_2023/EBioMedicine.myeloid.rds"))

Lin_fib <- readRDS(file.path(local_project_dir, "宫颈癌空转/ebiomedicine_scRNA_rds/成纤维细胞/Result.RDS"))
DimPlot(Lin_fib, reduction = "tsne",label = TRUE,raster=FALSE)
Lin_fib@meta.data$cell_subtype<-gsub('Other...','Other',Lin_fib@meta.data$cell_type)
Lin_fib <- subset(Lin_fib,meta_anno=='Fibroblast')
anno <- c("Fib-C1"="myoCAFs", "Fib-C2"="dCAFs", "Fib-C3"="dCAFs",
          "Fib-C4"="vCAFs", "Fib-C5"="vCAFs", "Fib-C6"="mCAFs", "Fib-C7"="mCAFs",
          "Fib-C8"="mCAFs", "Fib-C9"="mCAFs", "Fib-C10"="myoCAFs", "Fib-C11"="myoCAFs","Fib-C12"="myoCAFs")
Lin_fib[['cell_type']] = unname(anno[Lin_fib@meta.data$cell_type])
Idents(Lin_fib)<-factor(Lin_fib@meta.data[["cell_type"]],levels=unique(Lin_fib@meta.data[["cell_type"]]))
Lin_fib_markers <- FindAllMarkers(Lin_fib, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Lin_fib_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Lin_EBioMedicine_2023_fib_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

saveRDS(Lin_fib,file.path(local_project_dir, "宫颈癌数据库/单细胞/Lin_EBioMedicine_2023/EBioMedicine.fib.rds"))

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))


Lin_endo <- readRDS(file.path(local_project_dir, "宫颈癌空转/ebiomedicine_scRNA_rds/所有细胞/cc_comb2.filtered.RDS"))
Lin_endo@meta.data$group <- "NA"
Lin_endo@meta.data[grep("^Ade",Lin_endo@meta.data$orig.ident), ]$group <- "Ade"
Lin_endo@meta.data[grep("^P", Lin_endo@meta.data$orig.ident), ]$group <- "Scc"

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
saveRDS(Lin_endo,file.path(local_project_dir, "宫颈癌数据库/单细胞/Lin_EBioMedicine_2023/EBioMedicine.rds"))

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
saveRDS(Lin_endo,file.path(local_project_dir, "宫颈癌数据库/单细胞/Lin_EBioMedicine_2023/EBioMedicine.endo.rds"))

Lin_endo_markers <- FindAllMarkers(Lin_endo, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Lin_endo_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Lin_EBioMedicine_2023_endo_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 09: Qu et al., Cancer Communications (2023), GSE197461
# Original script: single-cell/Qu_Cancer.Commun_2023/GSE197461.R
# =============================================================================

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
# Sample integration, clustering, and annotation
mat <- list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461"),full.names = T)
sceList <- lapply(mat,function(x){
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samples<-list.files(file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461"))
type<-c('ADC','ADC','ADC','ADC','ADC','SCC','SCC','SCC')
orig <- NULL
group <- NULL
for(i in 1:length(sceList)){
  sce<-sceList[[i]]
  orig = rep(samples[i],ncol(sceList[[i]]))
  group = rep(type[i],ncol(sceList[[i]]))
  sceList[[i]][["orig.ident"]] = orig
  sceList[[i]][["group"]] = group
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
}

sce <- merge(sceList[[1]],
             y = sceList[-1],
             add.cell.ids = samples)
VlnPlot(sce,features = c('nCount_RNA','nFeature_RNA','percent.mt'))
sce <- subset(sce, subset = nFeature_RNA >= 200 & percent.mt <= 20)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst")
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- RunHarmony(sce,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(sce, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(sce, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(sce, ndims = 50, reduction = "harmony")
sce <- FindNeighbors(sce, dims = 1:10,reduction="harmony")
sce <- FindClusters(object = sce,resolution = seq(0.1,1,by=0.1))
clustree(sce@meta.data, prefix = "RNA_snn_res.")
sce <- FindClusters(sce, resolution = 0.8)
sce <- RunUMAP(sce,reduction="harmony", dims = 1:10)
sce <- RunTSNE(sce,reduction="harmony", dims = 1:10)
DimPlot(sce, reduction = "tsne", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','NCAM1','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9','FCN1',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','WFDC2', 'KRT8','KRT18','KRT19', 'KRT5','TP63','KLF5',   #epithelial cells
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2','SUSD2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2')                     #mast cells
DotPlot(sce, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))
DimPlot(sce, reduction = "tsne", label = TRUE,raster=FALSE)

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
DimPlot(sce, reduction = "umap",raster=FALSE)
saveRDS(sce, file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.sce.anno.rds"))


# Myeloid-cell subclustering
Hua.myeloid<-subset(sce,Cell_type =='Myeloid cells')
Hua.myeloid <- RunPCA(Hua.myeloid, features = VariableFeatures(object = Hua.myeloid))
Hua.myeloid <- RunHarmony(Hua.myeloid,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.myeloid, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.myeloid, ndims = 50, reduction = "harmony")
Hua.myeloid <- FindNeighbors(Hua.myeloid, dims = 1:30,reduction="harmony")
Hua.myeloid <- FindClusters(object = Hua.myeloid,resolution = seq(0.1,1,by=0.1))
clustree(Hua.myeloid@meta.data, prefix = "RNA_snn_res.")
Hua.myeloid <- FindClusters(Hua.myeloid, resolution = 1)
Hua.myeloid <- RunUMAP(Hua.myeloid,reduction="harmony", dims = 1:30)
Hua.myeloid <- RunTSNE(Hua.myeloid,reduction="harmony", dims = 1:30)
DimPlot(Hua.myeloid, reduction = "tsne",raster=FALSE,label = T)

marker<-c('XCR1','CLEC9A',
          'CD207','CD1C','FCER1A',
          'CLEC4C','LILRA4')
marker<-c('CSF3R',      #Neutrophils
          'FCN1','FCGR2A','S100A8','S100A9',     #Monocytes
          'LYZ','CD14','CD68','CD163','MS4A6A','C1QB','C1QA',  #Macrophages
          'IRF8','CD1C','LAMP3', 'GPR183','PLD4',      #DC
          'TOP2A','MKI67')
DotPlot(Hua.myeloid, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))

new.cluster.ids <- c('Macrophages','Monocytes','Dendritic cells','Macrophages','Macrophages','Monocytes','Macrophages',
                     'Proliferation','Macrophages','Neutrophils','Dendritic cells','Dendritic cells','Dendritic cells')
names(new.cluster.ids) <- levels(Hua.myeloid)
Hua.myeloid <- RenameIdents(Hua.myeloid, new.cluster.ids)
Hua.myeloid@meta.data[["Cell_subtype"]]<-as.character(Hua.myeloid@active.ident)
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
DimPlot(Hua.myeloid, reduction = "umap",raster=FALSE)
saveRDS(Hua.myeloid, file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.myeloid.anno.rds"))




# T-cell subclustering
Hua.NKT<-subset(sce,Cell_type =='T cells')
Hua.NKT <- RunPCA(Hua.NKT, features = VariableFeatures(object = Hua.NKT))
Hua.NKT <- RunHarmony(Hua.NKT,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.NKT, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.NKT, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.NKT, ndims = 50, reduction = "harmony")
Hua.NKT <- FindNeighbors(Hua.NKT, dims = 1:20,reduction="harmony")
Hua.NKT <- FindClusters(object = Hua.NKT,resolution = seq(0.1,1,by=0.1))
clustree(Hua.NKT@meta.data, prefix = "RNA_snn_res.")
Hua.NKT <- FindClusters(Hua.NKT, resolution = 0.6)
Hua.NKT <- RunUMAP(Hua.NKT,reduction="harmony", dims = 1:20)
DimPlot(Hua.NKT, reduction = "umap",raster=FALSE,label = T)

marker<-c('CD4','CD8A','CD8B','NKG7','TOP2A','MKI67')
DotPlot(Hua.NKT, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))



new.cluster.ids <- c('CD4','CD8','CD4','CD8','CD8','CD4','NK','CD8','CD4','CD8','NK','CD8','CD4')
names(new.cluster.ids) <- levels(Hua.NKT)
Hua.NKT <- RenameIdents(Hua.NKT, new.cluster.ids)
Hua.NKT@meta.data[["Cell_subtype"]]<-as.character(Hua.NKT@active.ident)
Idents(Hua.NKT)<-Hua.NKT@meta.data[["Cell_subtype"]]
DimPlot(Hua.NKT, reduction = "umap",raster=FALSE)

#CD4
Hua.CD4<-subset(Hua.NKT,Cell_subtype =='CD4')
Hua.CD4 <- RunPCA(Hua.CD4, features = VariableFeatures(object = Hua.CD4))
Hua.CD4 <- RunHarmony(Hua.CD4,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.CD4, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.CD4, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.CD4, ndims = 50, reduction = "harmony")
Hua.CD4 <- FindNeighbors(Hua.CD4, dims = 1:20,reduction="harmony")
Hua.CD4 <- FindClusters(object = Hua.CD4,resolution = seq(0.1,1,by=0.1))
clustree(Hua.CD4@meta.data, prefix = "RNA_snn_res.")
Hua.CD4 <- FindClusters(Hua.CD4, resolution = 0.7)
Hua.CD4 <- RunUMAP(Hua.CD4,reduction="harmony", dims = 1:20)
Hua.CD4 <- RunTSNE(Hua.CD4,reduction="harmony", dims = 1:20)
DimPlot(Hua.CD4, reduction = "umap",raster=FALSE,label = T)

marker<-c('FOXP3','CTLA4','TNFRSF18', #CD4 Treg
          'CCL5','GZMK','CXCR4',   #Tems
          'LEF1','CCR7','SELL',   #naive
          'IL17A','CCR6','CTSH',  #Th17
          'CXCL13','PDCD1','IFNG'  #Th1-like
)
marker<-c('FOXP3','CTLA4','IL2RA', #CD4 Treg
          'CXCL13','PDCD1','IFNG','CD200',  #Th1-like
          'CCL5','GZMK','CXCR4',   #Tems
          'LEF1','CCR7','SELL',   #naive
          'TOP2A','MKI67'
)
DotPlot(Hua.CD4, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))


new.cluster.ids <- c('CD4_naive','CD4_Treg','CD4_Treg','CD4_Tem','CD4_Tem','CD4_Th1-like','CD4_naive','CD4_naive','CD4_Treg','CD4_naive','CD4_naive','CD4_Tem')
names(new.cluster.ids) <- levels(Hua.CD4)
Hua.CD4 <- RenameIdents(Hua.CD4, new.cluster.ids)
Hua.CD4@meta.data[["Cell_subtype"]]<-as.character(Hua.CD4@active.ident)
Hua.CD4@meta.data[["Cell_subtype"]]<-factor(Hua.CD4@meta.data[["Cell_subtype"]],
                                            levels=c('Tem','Th1-like','Th17','Naive','Treg'))
Idents(Hua.CD4)<-Hua.CD4@meta.data[["Cell_subtype"]]
DimPlot(Hua.CD4, reduction = "umap",raster=FALSE)

saveRDS(Hua.CD4, file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.CD4T.anno.rds"))

#CD8
Hua.CD8<-subset(Hua.NKT,Cell_subtype =='CD8')
Hua.CD8 <- RunPCA(Hua.CD8, features = VariableFeatures(object = Hua.CD8))
Hua.CD8 <- RunHarmony(Hua.CD8,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.CD8, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.CD8, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.CD8, ndims = 50, reduction = "harmony")
Hua.CD8 <- FindNeighbors(Hua.CD8, dims = 1:20,reduction="harmony")
Hua.CD8 <- FindClusters(object = Hua.CD8,resolution = seq(0.1,1,by=0.1))
clustree(Hua.CD8@meta.data, prefix = "RNA_snn_res.")
Hua.CD8 <- FindClusters(Hua.CD8, resolution = 0.7)
Hua.CD8 <- RunUMAP(Hua.CD8,reduction="harmony", dims = 1:20)
Hua.CD8 <- RunTSNE(Hua.CD8,reduction="harmony", dims = 1:20)
DimPlot(Hua.CD8, reduction = "umap",raster=FALSE,label = T)

marker<-c('GZMK','CXCR4','CST7',    #CD8 Tem
          'XCL1','CAPG','NR4A1',    #Trm
          'CD160','KIR2DL4','KLRC2', #IEL
          'HAVCR2','CXCL13','PDCD1', #Tex
          'SLC4A10','KLRB1','ZBTB16', #MAIT
          'CX3CR1','FCGR3A','FGFBP2', #TemRA
          'LEF1','CCR7','SELL',  #naive T
          'TOP2A','MKI67'
)
DotPlot(Hua.CD8, features = marker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))+scale_color_gradientn(values = seq(0,1,0.2),colours = viridis_pal(option = 'H')(100))


new.cluster.ids <- c('CD8_Trm','CD8_Tem','CD8_Tex','CD8_Tex','CD8_Cycling','CD8_Trm','CD8_Trm','CD8_Trm','CD8_Trm','CD8_MAIT','CD8_Cycling','CD8_Tex')
names(new.cluster.ids) <- levels(Hua.CD8)
Hua.CD8 <- RenameIdents(Hua.CD8, new.cluster.ids)
Hua.CD8@meta.data[["Cell_subtype"]]<-as.character(Hua.CD8@active.ident)
Idents(Hua.CD8)<-Hua.CD8@meta.data[["Cell_subtype"]]
DimPlot(Hua.CD8, reduction = "umap",raster=FALSE)

saveRDS(Hua.CD8, file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.CD8T.anno.rds"))



# Epithelial-cell subclustering
Hua.epi<-subset(sce,Cell_type =='Epithelial cells')
Hua.epi <- RunPCA(Hua.epi, features = VariableFeatures(object = Hua.epi))
Hua.epi <- RunHarmony(Hua.epi,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Hua.epi, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Hua.epi, ndims = 50, reduction = "harmony")
Hua.epi <- FindNeighbors(Hua.epi, dims = 1:30,reduction="harmony")
Hua.epi <- FindClusters(object = Hua.epi,resolution = seq(0.1,1,by=0.1))
clustree(Hua.epi@meta.data, prefix = "RNA_snn_res.")
Hua.epi <- FindClusters(Hua.epi, resolution = 1)
Hua.epi <- RunUMAP(Hua.epi,reduction="harmony", dims = 1:30)
Hua.epi <- RunTSNE(Hua.epi,reduction="harmony", dims = 1:30)
DimPlot(Hua.epi, reduction = "tsne",raster=FALSE,label = T)

saveRDS(Hua.epi, file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.epi.rds"))


####Marker gene
load(file.path(user_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

Hua.CD4 <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.CD4T.anno.rds"))
Idents(Hua.CD4)<-Hua.CD4@meta.data[["Cell_subtype"]]
Hua.CD4<-JoinLayers(Hua.CD4)
Hua.CD4_markers <- FindAllMarkers(Hua.CD4, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.CD4_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Qu_Cancer.Commun_2023_CD4T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.CD8 <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.CD8T.anno.rds"))
Idents(Hua.CD8)<-Hua.CD8@meta.data[["Cell_subtype"]]
Hua.CD8<-JoinLayers(Hua.CD8)
Hua.CD8_markers <- FindAllMarkers(Hua.CD8, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.CD8_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Qu_Cancer.Commun_2023_CD8T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Hua.myeloid <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Qu_Cancer.Commun_2023/GSE197461.myeloid.anno.rds"))
Idents(Hua.myeloid)<-Hua.myeloid@meta.data[["Cell_subtype"]]
Hua.myeloid<-JoinLayers(Hua.myeloid)
Hua.myeloid_markers <- FindAllMarkers(Hua.myeloid, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Hua.myeloid_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Qu_Cancer.Commun_2023_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 10: Zhang et al., EBioMedicine (2023)
# Original script: single-cell/Zhang_EBioMedicine_2023/forDB.R
# =============================================================================

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
#########Huang Ebiomedicine
Huang_scRNA <- readRDS(file.path(user_project_dir, "宫颈癌空转/Huang_scRNA.rds"))
Huang_scRNA <- RunTSNE(Huang_scRNA,reduction="pca", dims = 1:20)
DimPlot(Huang_scRNA, reduction = "tsne",label = TRUE,raster=FALSE,group.by = 'cell_type')
saveRDS(Huang_scRNA,file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_scRNA.rds"))


Huang_epi <- readRDS(file.path(local_project_dir, "宫颈癌空转/分析/Huang_epi.rds"))
Huang_epi <- RunTSNE(Huang_epi,reduction="pca", dims = 1:20)
DimPlot(Huang_epi, reduction = "tsne",label = TRUE,raster=FALSE,group.by = 'epi_specific1')
saveRDS(Huang_epi,file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_epi.rds"))

Huang_T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_tcell.rds"))
DimPlot(Huang_T, reduction = "tsne",label = TRUE,raster=FALSE,group.by = 'annotation')

Huang_fib <- subset(Huang_scRNA,cell_type=='Fibroblasts')
Huang_fib <- NormalizeData(Huang_fib, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_fib <- FindVariableFeatures(Huang_fib, selection.method = "vst", nfeatures = 3000)
Huang_fib <- ScaleData(Huang_fib, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_fib))
Huang_fib <- RunPCA(Huang_fib, features = VariableFeatures(object = Huang_fib))
Huang_fib <- RunHarmony(Huang_fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Huang_fib, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Huang_fib, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Huang_fib, ndims = 50, reduction = "harmony")
Huang_fib <- FindNeighbors(Huang_fib, dims = 1:20,reduction="harmony")
Huang_fib <- FindClusters(object = Huang_fib,resolution = seq(0.1,1,by=0.1))
clustree(Huang_fib@meta.data, prefix = "RNA_snn_res.")
Huang_fib <- RunUMAP(Huang_fib,reduction="harmony", dims = 1:20)
Huang_fib <- RunTSNE(Huang_fib,reduction="harmony", dims = 1:20)
DimPlot(Huang_fib, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Huang_fib,file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_fib.rds"))


Huang_endo <- subset(Huang_scRNA,cell_type=='Endothelial cell')
Huang_endo <- NormalizeData(Huang_endo, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_endo <- FindVariableFeatures(Huang_endo, selection.method = "vst", nfeatures = 3000)
Huang_endo <- ScaleData(Huang_endo, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_endo))
Huang_endo <- RunPCA(Huang_endo, features = VariableFeatures(object = Huang_endo))
Huang_endo <- RunHarmony(Huang_endo,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Huang_endo, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Huang_endo, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Huang_endo, ndims = 50, reduction = "harmony")
Huang_endo <- FindNeighbors(Huang_endo, dims = 1:20,reduction="harmony")
Huang_endo <- FindClusters(object = Huang_endo,resolution = seq(0.1,1,by=0.1))
clustree(Huang_endo@meta.data, prefix = "RNA_snn_res.")
Huang_endo <- RunUMAP(Huang_endo,reduction="harmony", dims = 1:20)
Huang_endo <- RunTSNE(Huang_endo,reduction="harmony", dims = 1:20)
DimPlot(Huang_endo, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Huang_endo,file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_endo.rds"))


Huang_macro <- subset(Huang_scRNA,cell_type=='Macrophages')
Huang_macro <- NormalizeData(Huang_macro, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_macro <- FindVariableFeatures(Huang_macro, selection.method = "vst", nfeatures = 3000)
Huang_macro <- ScaleData(Huang_macro, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_macro))
Huang_macro <- RunPCA(Huang_macro, features = VariableFeatures(object = Huang_macro))
Huang_macro <- RunHarmony(Huang_macro,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Huang_macro, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Huang_macro, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Huang_macro, ndims = 50, reduction = "harmony")
Huang_macro <- FindNeighbors(Huang_macro, dims = 1:20,reduction="harmony")
Huang_macro <- FindClusters(object = Huang_macro,resolution = seq(0.1,1,by=0.1))
clustree(Huang_macro@meta.data, prefix = "RNA_snn_res.")
Huang_macro <- RunUMAP(Huang_macro,reduction="harmony", dims = 1:20)
Huang_macro <- RunTSNE(Huang_macro,reduction="harmony", dims = 1:20)
DimPlot(Huang_macro, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Huang_macro,file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_macro.rds"))


Huang_B <- subset(Huang_scRNA,cell_type=='B cells')
Huang_B <- NormalizeData(Huang_B, normalization.method = "LogNormalize", scale.factor = 10000)
Huang_B <- FindVariableFeatures(Huang_B, selection.method = "vst", nfeatures = 3000)
Huang_B <- ScaleData(Huang_B, vars.to.regress = "percent.mt",features = VariableFeatures(object = Huang_B))
Huang_B <- RunPCA(Huang_B, features = VariableFeatures(object = Huang_B))
Huang_B <- RunHarmony(Huang_B,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
DimPlot(Huang_B, group.by = 'orig.ident', reduction = "pca",raster=FALSE)
DimPlot(Huang_B, group.by = 'orig.ident', reduction = "harmony",raster=FALSE)
ElbowPlot(Huang_B, ndims = 50, reduction = "harmony")
Huang_B <- FindNeighbors(Huang_B, dims = 1:20,reduction="harmony")
Huang_B <- FindClusters(object = Huang_B,resolution = seq(0.1,1,by=0.1))
clustree(Huang_B@meta.data, prefix = "RNA_snn_res.")
Huang_B <- RunUMAP(Huang_B,reduction="harmony", dims = 1:20)
Huang_B <- RunTSNE(Huang_B,reduction="harmony", dims = 1:20)
DimPlot(Huang_B, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(Huang_B,file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_B.rds"))



####Marker gene
subcelltype_marker_list<-list()
Huang_epi <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_epi.anno.rds"))
Idents(Huang_epi)<-Huang_epi@meta.data[["epi_specific1"]]
Huang_epi_markers <- FindAllMarkers(Huang_epi, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Huang_epi_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Zhang_EBioMedicine_2023_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Huang_T <- readRDS(file.path(local_project_dir, "宫颈癌数据库/单细胞/Zhang_EBioMedicine_2023/Huang_tcell.anno.rds"))
Idents(Huang_T)<-Huang_T@meta.data[["annotation"]]
Huang_T_markers <- FindAllMarkers(Huang_T, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Huang_T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Zhang_EBioMedicine_2023_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(local_project_dir, "宫颈癌数据库/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 11: Sandoval et al., Cancer Research (2026), GSE297041
# Original script: new-data/single-cell/Sandoval_Cancer.Res_2026(GSE297041)/GSE297041.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(clustree)
seu<-readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Sandoval_Cancer.Res_2026(GSE297041)/GSE297041_CESC_18_scRNA_rmdoublet_0.2_cluster_ident_without_anchor.rds"))
seu <- RunTSNE(seu,dims = 1:30)
seu$group<-seu$category
DimPlot(seu, reduction = "umap", label = TRUE,raster=FALSE)
new.cluster.ids <- c('Immune','Immune','Immune','Immune','Tumor','Immune','Tumor','Tumor','Tumor','Tumor','Immune','Tumor','Stroma',
                     'Immune','Tumor','Stroma','Immune','Stroma','Tumor','Tumor','Tumor','Tumor','Immune','Immune','Immune')
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu@meta.data[["cell_type"]]<-seu@active.ident
seu@meta.data[["cell_type"]]<-factor(seu@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Sandoval_Cancer.Res_2026(GSE297041)/GSE297041.sce.anno.rds"))

# Immune-cell subclustering
Immune<-subset(seu,subset=cell_type =='Immune')
Immune <- SCTransform(Immune, vars.to.regress = c('nCount_RNA', 'percent.mito'),variable.features.n = 2000)
Immune <- RunPCA(Immune, features = VariableFeatures(object = Immune))
# Immune-cell subclustering <- RunHarmony(Immune,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Immune <- FindNeighbors(Immune, dims = 1:30,reduction="pca")
Immune <- FindClusters(object = Immune,resolution = seq(0.1,1,by=0.1))
clustree(Immune@meta.data, prefix = "RNA_snn_res.")
Immune <- FindClusters(Immune, resolution = 0.5)
Immune <- RunUMAP(Immune,reduction="pca", dims = 1:30)
Immune <- RunTSNE(Immune,dims = 1:30)
Immune$group<-Immune$category
DimPlot(Immune, reduction = "umap", label = TRUE,raster=FALSE)
marker<-c('CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCN1','VCAN','FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183','CD74',
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
marker<-c('HCAR2','HCAR3','CXCR2','CCL3','CCL4','SPP1','TREM2')
DotPlot(Immune, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))


new.cluster.ids <- c('CD8+ Tcells','HCAR2/3+ Neutrophils','CD8+ Tcells','Monocytes','Treg','CXCR2+ Neutrophils','SPP1+ Macrophages',
                     'CD4+ Tcells','DCs','TREM2+ Macrophages','CCL3/4+ Neutrophils','Plasma cells','CD8+ Tcells','B cells',
                     'Unknown','Proliferating Tcells','Plasma cells','Mast','SPP1+ Macrophages','Plasma cells','pDC','Monocytes')
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
Immune@meta.data[["cell_subtype"]]<-Immune@active.ident
Immune@meta.data[["cell_subtype"]]<-factor(Immune@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Immune, reduction = "umap", label = T,raster=FALSE)
saveRDS(Immune,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Sandoval_Cancer.Res_2026(GSE297041)/GSE297041.immune.anno.rds"))


# Export subtype marker genes
load(file.path(project_dir, "单细胞/subcelltype_marker_list.RData"))
Immune <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Sandoval_Cancer.Res_2026(GSE297041)/GSE297041.immune.anno.rds"))
Idents(Immune)<-Immune@meta.data[["cell_subtype"]]
Immune_markers <- FindAllMarkers(Immune, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Immune_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Sandoval_Cancer.Res_2026_immune_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 12: Yuan et al., Frontiers in Immunology (2026), GSE308792
# Original script: new-data/single-cell/Yuan_Front.Immunol_2026(GSE308792)/GSE308792.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(clustree)
samplename <- list.dirs(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Yuan_Front.Immunol_2026(GSE308792)/Results/"),recursive=F,full.names =F)
files<-paste0(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Yuan_Front.Immunol_2026(GSE308792)/Results/"),samplename,'/',samplename,'_RSEC_MolsPerCell_MEX')
seuList <- lapply(files,function(x){ 
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
names(seuList)<-samplename
group<-rep(c('Para-tumor','SCC'),times=3)
VlnPlot(seuList[[6]],features = c('nCount_RNA','nFeature_RNA'))

for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplename[i],ncol(seuList[[i]]))
  Group = rep(group[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["group"]] = Group
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
  sweep.stats[order(sweep.stats$BCreal),]
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
  print(paste0(samplename[i]," Finished!"))
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = samplename)
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
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "RNA_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('T cell','T cell','Epithelial cell','Epithelial cell','Epithelial cell','Epithelial cell',
                     'Epithelial cell','Epithelial cell','Epithelial cell','Fibroblast','Macrophages','Plasma cell',
                     'Endothelial cell','B cell','Epithelial cell','pDC','Epithelial cell','Epithelial cell',
                     'Epithelial cell','Mast cell')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Yuan_Front.Immunol_2026(GSE308792)/GSE308792.sce.anno.rds"))


#T
Tcell<-subset(seu_merge,subset=cell_type == 'T cell')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:30,reduction="harmony")
Tcell <- FindClusters(object = Tcell,resolution = seq(0.1,1,by=0.1))
clustree(Tcell@meta.data, prefix = "RNA_snn_res.")
Tcell <- FindClusters(Tcell, resolution = 0.3)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:30)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:30)
DimPlot(Tcell, reduction = "umap", label = TRUE,raster=FALSE)
marker<-c('CD4','CD8A','CD8B','TOP2A','MKI67','IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7','NCAM1','FCGR3A','FCGR3B')
DotPlot(Tcell, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('CD8 Tex','CD4 Treg','CD4 Naive T','CD8 Naive T','CD4 Naive T','CD8 Tex','CD16+ NK',
                     'CD4 Tph','Memory T cells','CD56+ NK','Unknown')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids))[c(2:4,6:8,1,5,9)])
DimPlot(Tcell, reduction = "umap", label = T,raster=FALSE)
T_marker<-FindAllMarkers(JoinLayers(Tcell),only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
saveRDS(Tcell,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Yuan_Front.Immunol_2026(GSE308792)/GSE308792.T.anno.rds"))


# Export subtype marker genes
load(file.path(project_dir, "单细胞/subcelltype_marker_list.RData"))
Tcell <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Yuan_Front.Immunol_2026(GSE308792)/GSE308792.T.anno.rds"))
Idents(Tcell)<-Tcell@meta.data[["cell_subtype"]]
Tcell<-JoinLayers(Tcell)
T_markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Yuan_Front.Immunol_2026_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 13: Wu et al., Communications Biology (2025), E-MTAB-15983
# Original script: new-data/single-cell/Wu_Commun.Biol_2025/E-MTAB-15983.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(clustree)
files<-list.files(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Wu_Commun.Biol_2025/E-MTAB-15983"),full.names = T)
seuList <- lapply(files,function(x){
  seu = CreateSeuratObject(counts = Read10X_h5(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samplename<-gsub('_sample_feature_bc_matrix.h5','',basename(files))
group<-c(rep('AC',8),rep('ASC',2),rep('SCC',6))
names(seuList)<-samplename
# TODO: Sample 1 is skipped here, but downstream code expects DoubletFinder fields for every sample.
for(i in 2:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplename[i],ncol(seuList[[i]]))
  Group = rep(group[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["Type"]] = Group
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
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  # Estimate the homotypic doublet proportion     
  homotypic.prop <- modelHomotypic(sce$seurat_clusters) 
  # Estimate the expected doublet count
  DoubletRate = ncol(sce)*8*1e-6 
  nExp_poi <- round(DoubletRate *nrow(sce@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) # Adjust for homotypic doublets
  # TODO: This replaces the homotypic-adjusted estimate calculated immediately above.
nExp_poi.adj<-ncol(sce)*ncol(sce)*7.6*1e-6 
  # Identify doublets with the selected pK
  sce <- doubletFinder(sce, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
  seuList[[i]]<-sce
  print(paste0(samplename[i]," Finished!"))
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}
DimPlot(seuList[[2]], group.by = 'DF.classifications', label = TRUE,raster=FALSE)

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = samplename)
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
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "SCT_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E',                           #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','CSF3R','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'GPR183','PLD4',              #dendritic cells (DCs)
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'MKI67','TOP2A',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2'
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('NK/T','Epithelial','NK/T','Fibroblast','Myeloid','NK/T','Endothelial','Epithelial','B','CAF','Epithelial',
                     'NK/T','Plasma','Mast cell','Myeloid','Epithelial','Myeloid','Epithelial','Fibroblast','Fibroblast','Myeloid',
                     'Epithelial','Epithelial','Endothelial','Myeloid')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = F,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Wu_Commun.Biol_2025/E-MTAB-15983.sce.anno.rds"))

Epithelial<-subset(Epithelial,subset=cell_type=='Epithelial')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:30,reduction="harmony")
Epithelial <- FindClusters(object = Epithelial,resolution = seq(0.1,1,by=0.1))
clustree(Epithelial@meta.data, prefix = "RNA_snn_res.")
Epithelial <- FindClusters(Epithelial, resolution = 0.2)
Epithelial <- RunUMAP(Epithelial,reduction="harmony", dims = 1:30)
Epithelial <- RunTSNE(Epithelial,reduction="harmony", dims = 1:30)
DimPlot(Epithelial, reduction = "tsne", label = TRUE,raster=FALSE)
DotPlot(Epithelial, features = c('MKI67','KRT5','S100A9','KRT6A'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
sample_ratio<-melt(table(Epithelial@meta.data$Type,paste0('Epi',Epithelial@meta.data$seurat_clusters)))
colnames(sample_ratio)<-c('Stage','Cell_type','Count')
ggplot(data = sample_ratio, aes(x = Stage, y = Count, fill = Cell_type)) +
  geom_bar(stat = "identity", width=0.8,aes(group=Cell_type),position="fill")+
  #scale_fill_manual(values = color)+
  theme_classic()+
  labs(x="",y="Cell proportion")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
Epithelial$Cell_subtype<-factor(paste0('Epi',Epithelial@meta.data$seurat_clusters),levels = paste0('Epi',levels(Epithelial@meta.data$seurat_clusters)))
Idents(Epithelial)<-Epithelial$Cell_subtype
saveRDS(Epithelial,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Wu_Commun.Biol_2025/E-MTAB-15983.epi.anno.rds"))

Immune<-subset(seu_merge,subset=cell_type %in% c('NK/T','Myeloid','B','Plasma'))
Immune <- NormalizeData(Immune, normalization.method = "LogNormalize", scale.factor = 10000)
Immune <- FindVariableFeatures(Immune, selection.method = "vst", nfeatures = 2000)
Immune <- ScaleData(Immune, vars.to.regress = "percent.mt")
Immune <- RunPCA(Immune, features = VariableFeatures(object = Immune))
Immune <- RunHarmony(Immune,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Immune <- FindNeighbors(Immune, dims = 1:30,reduction="harmony")
Immune <- FindClusters(object = Immune,resolution = seq(0.1,1,by=0.1))
clustree(Immune@meta.data, prefix = "RNA_snn_res.")
Immune <- FindClusters(Immune, resolution = 0.6)
Immune <- RunUMAP(Immune,reduction="harmony", dims = 1:30)
Immune <- RunTSNE(Immune,reduction="harmony", dims = 1:30)
DimPlot(Immune, reduction = "umap", label = TRUE,raster=FALSE)
marker<-c('CD4','CD8A','CD8B','NKG7','NCAM1','FCGR3A',
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',
          'CD1C','LAMP3', 'GPR183','PLD4', 
          'JCHAIN','CD38', 'IGLC2','IGHG4',
          'MS4A1','CD79B','IGKC','CD79A','MZB1',
          'TOP2A','MKI67')
DotPlot(Immune, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))


new.cluster.ids <- c('CD8 T','CD4 Naive','CD8 T','CD4 Treg','CD163 Mac','CD163 Mac','B','CD4 Tph','NK','Proliferating NK/CD8 T','Plasma','Neu',
                     'CD8 T','Mono/Macro','Mono/Macro','Proliferating Mono/Mac','Mono/Macro','Mono/Macro','Mono/Macro','CD4 Treg','pDC')
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
Immune@meta.data[["cell_subtype"]]<-Immune@active.ident
Immune@meta.data[["cell_subtype"]]<-factor(Immune@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Immune, reduction = "umap", label = F,raster=FALSE)
saveRDS(Immune,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Wu_Commun.Biol_2025/E-MTAB-15983.immune.anno.rds"))


# Export subtype marker genes
load(file.path(project_dir, "单细胞/subcelltype_marker_list.RData"))
subcelltype_marker_list[144:147]<-NULL
Epithelial <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Wu_Commun.Biol_2025/E-MTAB-15983.epi.anno.rds"))
Idents(Epithelial)<-Epithelial@meta.data[["Cell_subtype"]]
Epithelial<-JoinLayers(Epithelial)
Epithelial_markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Epithelial_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Wu_Commun.Biol_2025_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Immune <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Wu_Commun.Biol_2025/E-MTAB-15983.immune.anno.rds"))
Idents(Immune)<-Immune@meta.data[["cell_subtype"]]
Immune<-JoinLayers(Immune)
Immune_markers <- FindAllMarkers(Immune, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Immune_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Wu_Commun.Biol_2025_immune_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 14: Hyeon et al., Molecular Cancer (2025), GSE279998
# Original script: new-data/single-cell/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(stringr)
library(clustree)
sample_dirs <- list.dirs(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998_RAW"),recursive=F)
for(dir in sample_dirs){
  files <- list.files(dir, full.names = TRUE)
  # Standardize the matrix filename
  mat <- files[str_detect(basename(files), "matrix\\.mtx")]
  if(length(mat) == 1){
    file.rename(mat, file.path(dir, "matrix.mtx.gz"))
  }
  # Standardize the feature filename
  fea <- files[str_detect(basename(files), "features\\.tsv")]
  if(length(fea) == 1){
    file.rename(fea, file.path(dir, "features.tsv.gz"))
  }
  # Standardize the barcode filename
  bar <- files[str_detect(basename(files), "barcodes\\.tsv")]
  if(length(bar) == 1){
    file.rename(bar, file.path(dir, "barcodes.tsv.gz"))
  }
}
seuList <- lapply(sample_dirs,function(x){ 
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samplename<-basename(sample_dirs)
group<-c(rep('ADC',3),'ADS',rep('SCC',12))
names(seuList)<-samplename

for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplename[i],ncol(seuList[[i]]))
  Group = rep(group[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["Type"]] = Group
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
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
  seuList[[i]]<-sce
  print(paste0(samplename[i]," Finished!"))
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}
DimPlot(seuList[[2]], group.by = 'DF.classifications', label = TRUE,raster=FALSE)

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = samplename)
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
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "RNA_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:24)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:24)
seu_merge$group<-seu_merge$Type
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Epithelial','CD8 T cell','Treg','CD4 T cell','NK cell','CD8 T cell','Plasma cell','CAF','CD8 T cell',
                     'Macrophage','B cell','CAF/SMC','Endothelial','Plasma cell','Neutrophil','Proliferating T cell',
                     'Mast','Epithelial','pDC','Plasma cell','Plasma cell')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.sce.anno.rds"))

# Epithelial-cell subclustering
Epithelial<-subset(seu_merge,subset=cell_type=='Epithelial')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:50,reduction="harmony")
Epithelial <- FindClusters(object = Epithelial,resolution = seq(0.1,1,by=0.1))
clustree(Epithelial@meta.data, prefix = "RNA_snn_res.")
Epithelial <- FindClusters(Epithelial, resolution = 0.3)
Epithelial <- RunUMAP(Epithelial,reduction="harmony", dims = 1:50)
Epithelial <- RunTSNE(Epithelial,reduction="harmony", dims = 1:50)
Epithelial$group<-Epithelial$Type
DimPlot(Epithelial, reduction = "umap", label = TRUE,raster=FALSE)
DotPlot(Epithelial, features = c('KRT1','KRT3','GADD45B','AKR1C2','PGK1','CXCL10','MUC5AC','CEACAM5','TOP2A','MKI67'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('AKR1C2','GADD45B','CEACAM5','KRT1/3','Proliferating','GADD45B','Unknown','MUC5AC','PGK1','Unknown','Unknown','CXCL10')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Epithelial, reduction = "umap", label = T,raster=FALSE)
saveRDS(Epithelial,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.epi.anno.rds"))


#CAF
Fib<-subset(seu_merge,subset=cell_type=='CAF')
Fib <- NormalizeData(Fib, normalization.method = "LogNormalize", scale.factor = 10000)
Fib <- FindVariableFeatures(Fib, selection.method = "vst", nfeatures = 2000)
Fib <- ScaleData(Fib, vars.to.regress = "percent.mt")
Fib <- RunPCA(Fib, features = VariableFeatures(object = Fib))
Fib <- RunHarmony(Fib,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Fib <- FindNeighbors(Fib, dims = 1:50,reduction="harmony")
Fib <- FindClusters(object = Fib,resolution = seq(0.1,1,by=0.1))
clustree(Fib@meta.data, prefix = "RNA_snn_res.")
Fib <- FindClusters(Fib, resolution = 0.7)
Fib <- RunUMAP(Fib,reduction="harmony", dims = 1:50)
Fib <- RunTSNE(Fib,reduction="harmony", dims = 1:50)
Fib$group<-Fib$Type
DimPlot(Fib, reduction = "umap", label = TRUE,raster=FALSE)
DotPlot(Fib, features = c('PDGFD','CXCL1','CXCL10','CCL21','PLN','RGS5','DES','TOP2A','MKI67'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('PDGFD','CCN5','CXCL1','PDGFD','PLN','PDGFD','CCL21','Unknown','DES','PDGFD','RGS5','CXCL10','Unknown','Unknown')
names(new.cluster.ids) <- levels(Fib)
Fib <- RenameIdents(Fib, new.cluster.ids)
Fib@meta.data[["cell_subtype"]]<-Fib@active.ident
Fib@meta.data[["cell_subtype"]]<-factor(Fib@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Fib, reduction = "umap", label = T,raster=FALSE)
Fib_marker<-FindAllMarkers(JoinLayers(Fib),only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
saveRDS(Fib,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.fib.anno.rds"))


# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type %in% c('Macrophage','pDC','Mast'))
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:27,reduction="harmony")
Mye <- FindClusters(object = Mye,resolution = seq(0.1,1,by=0.1))
clustree(Mye@meta.data, prefix = "RNA_snn_res.")
Mye <- FindClusters(Mye, resolution = 0.2)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:27)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:27)
Mye$group<-Mye$Type
DimPlot(Mye, reduction = "umap", label = TRUE,raster=FALSE)
DotPlot(Mye, features = c('CD1C','LAMP3','CLEC9A',
                          'C1QB','FCGR3B','FCN1'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('C1QB+ Mac','FCN1+ Mac','Mast','CD1C+ DC','C1QB+ Mac','pDC','LAMP3+ DC','CLEC9A+ Mac','Proliferating','C1QB+ Mac')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Mye, reduction = "umap", label = T,raster=FALSE)
saveRDS(Mye,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.myeloid.anno.rds"))

#T
Tcell<-subset(seu_merge,subset=cell_type %in% c('CD8 T cell','Treg','CD4 T cell','NK cell'))
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:21,reduction="harmony")
Tcell <- FindClusters(object = Tcell,resolution = seq(0.1,1,by=0.1))
clustree(Tcell@meta.data, prefix = "RNA_snn_res.")
Tcell <- FindClusters(Tcell, resolution = 0.3)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:21)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:21)
Tcell$group<-Tcell$Type
DimPlot(Tcell, reduction = "umap", label = TRUE,raster=FALSE)
DotPlot(Tcell, features = c('IL7R', 'GZMK','IFIT3','LAG3','TOP2A','MKI67'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('IL7R+ T','GZMK+ T','LAG3+ T','Treg','FCER1G+ NK','IL7R+ T','FCGR3A+ NK','CD8+ T','IFIT3+ T')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Tcell, reduction = "umap", label = T,raster=FALSE)
saveRDS(Tcell,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.T.anno.rds"))


# Export subtype marker genes
load(file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))
Epithelial <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.epi.anno.rds"))
Idents(Epithelial)<-Epithelial@meta.data[["cell_subtype"]]
Epithelial<-JoinLayers(Epithelial)
Epithelial_markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Epithelial_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Hyeon_Mol.Cancer_2025_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Fib <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.fib.anno.rds"))
Idents(Fib)<-Fib@meta.data[["cell_subtype"]]
Fib<-JoinLayers(Fib)
Fib_markers <- FindAllMarkers(Fib, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Fib_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Hyeon_Mol.Cancer_2025_fib_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Tcell <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.T.anno.rds"))
Idents(Tcell)<-Tcell@meta.data[["cell_subtype"]]
Tcell<-JoinLayers(Tcell)
T_markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Hyeon_Mol.Cancer_2025_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Mye <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.myeloid.anno.rds"))
Idents(Mye)<-Mye@meta.data[["cell_subtype"]]
Mye<-JoinLayers(Mye)
Mye_markers <- FindAllMarkers(Mye, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Mye_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Hyeon_Mol.Cancer_2025_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 15: Liu et al., Journal of Experimental and Clinical Cancer Research (2023), SCP1950
# Original script: new-data/single-cell/Liu_J.Exp.Clin.Cancer.Res_2023/SCP1950.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(stringr)
library(clustree)
sample_dirs<-list.dirs(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Liu_J.Exp.Clin.Cancer.Res_2023/SCP1950"),recursive=F)
for(dir in sample_dirs){
  files <- list.files(dir, full.names = TRUE)
  # Standardize the matrix filename
  mat <- files[str_detect(basename(files), "matrix\\.mtx")]
  if(length(mat) == 1){
    file.rename(mat, file.path(dir, "matrix.mtx.gz"))
  }
  # Standardize the feature filename
  fea <- files[str_detect(basename(files), "genes\\.tsv")]
  if(length(fea) == 1){
    file.rename(fea, file.path(dir, "features.tsv.gz"))
  }
  # Standardize the barcode filename
  bar <- files[str_detect(basename(files), "barcodes\\.tsv")]
  if(length(bar) == 1){
    file.rename(bar, file.path(dir, "barcodes.tsv.gz"))
  }
}
seuList <- lapply(sample_dirs,function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samplename<-basename(sample_dirs)
group<-c(rep('I',4),rep('II',3))
names(seuList)<-samplename

for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplename[i],ncol(seuList[[i]]))
  Group = rep(group[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["Type"]] = Group
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
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
  seuList[[i]]<-sce
  print(paste0(samplename[i]," Finished!"))
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = samplename)
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
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "RNA_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.5)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)
seu_merge$group<-seu_merge$Type
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
marker<- c(
  "NELL2", "ITK", "IL7R",
  "SIGLEC1", "FPR3", "MSR1",
  "TOP2A", "ASPM", "CENPF",
  "DCC", "CD38", "POU2AF1",
  "CTLA4", "IL2RA", "F5",
  "KLRC1", "GNLY", "NCR1",
  "BLK", "FCRL1", "MS4A1",
  "SLC24A4", "FLT3", "ZNF366",
  "CLEC4C",
  "CLDN10-AS1", "TOX3", "PROM1",
  "ADAMTS2", "LAMA2", "ABI3BP",
  "RRM2", "EXO1", "SKA3",
  "FLT1", "PCDH17", "VWF",
  "NRXN1", "ADAMTSL1", "PPP2R2B",
  "CPA3", "MS4A2", "TPSAB1",
  "KRT15", "KRT17", "KRT5",
  "IL1RN", "GPRC5A", "SPRR1B",
  "FLT4", "PROX1", "CD34",
  "NTRK2", "PTPRZ1", "CHL1",
  "SPARC", "COL1A1", "COL6A2",
  "RGS5", "ABCC9", "ADGRF5",
  "FOSB", "ATF3", "ITGB4"
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Basal cells','Mesenchymal stem cells','Basal cells','Cancer stem cells','Cancer stem cells',
                     'CD8 T cells','Proliferating Epithelial','Macrophage','Epithelial cells','T cells','naive B cells',
                     'Treg','myofibroblast','Cancer stem cells','mature B cells','Epithelial cells','Endothelial cells',
                     'Pericytes','T cells','CD141+CLEC9A+ DC','γδT','Lymphatic endothelial cells')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Liu_J.Exp.Clin.Cancer.Res_2023/SCP1950.sce.anno.rds"))

# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type %in% c('Macrophage','CD141+CLEC9A+ DC'))
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:30,reduction="harmony")
Mye <- FindClusters(object = Mye,resolution = seq(0.1,1,by=0.1))
clustree(Mye@meta.data, prefix = "RNA_snn_res.")
Mye <- FindClusters(Mye, resolution = 0.6)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:30)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:30)
DimPlot(Mye, reduction = "umap", label = TRUE,raster=FALSE)
DotPlot(Mye, features = c('MS4A6A', 'CD163', 'CD163L1',
                          'PARD3','EGFR','SMAD3',
                          'SLC16A10','SLC11A1','CTSL',
                          'ADAM19', 'HDAC9', 'MCOLN2','SLC8A1', 'RUNX3', 'LCP1',
                          'SEMA3A', 'ESRRG','IL18'),group.by = 'RNA_snn_res.0.1')+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('C0-Res','','','C3-DC','C2-M2','C1-TAM','','C4-M1','','','','C3-DC','','','')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Mye, reduction = "umap", label = T,raster=FALSE)
# TODO: This original output path points to the Hyeon dataset rather than SCP1950.
saveRDS(Mye,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Hyeon_Mol.Cancer_2025(GSE279998)/GSE279998.myeloid.anno.rds"))

# =============================================================================
# Dataset 16: Dai et al., Cell Reports Medicine (2024), GSE236738
# Original script: new-data/single-cell/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(stringr)
library(clustree)
sample_dirs <- list.dirs(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738_RAW"),recursive=F)
for(dir in sample_dirs){
  files <- list.files(dir, full.names = TRUE)
  # Standardize the matrix filename
  mat <- files[str_detect(basename(files), "matrix\\.mtx")]
  if(length(mat) == 1){
    file.rename(mat, file.path(dir, "matrix.mtx.gz"))
  }
  # Standardize the feature filename
  fea <- files[str_detect(basename(files), "features\\.tsv")]
  if(length(fea) == 1){
    file.rename(fea, file.path(dir, "features.tsv.gz"))
  }
  # Standardize the barcode filename
  bar <- files[str_detect(basename(files), "barcodes\\.tsv")]
  if(length(bar) == 1){
    file.rename(bar, file.path(dir, "barcodes.tsv.gz"))
  }
}
seuList <- lapply(sample_dirs,function(x){
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samplename<-basename(sample_dirs)
group<-rep(c('Pre','Post'),3)
names(seuList)<-samplename

for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplename[i],ncol(seuList[[i]]))
  Group = rep(group[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["Type"]] = Group
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
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
  seuList[[i]]<-sce
  print(paste0(samplename[i]," Finished!"))
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = samplename)
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
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "RNA_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.4)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:50)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:50)
seu_merge$group<-seu_merge$Type
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('T cell','T cell','Monocyte','T cell','T cell','Plasma cell','Fibroblast','Tumor cell','Tumor cell',
                     'Tumor cell','T cell','Tumor cell','Neutrophil','Mast cell','B cell','Fibroblast','T cell','Tumor cell',
                     'Blood cell','Endothelial cell','pDC')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.sce.anno.rds"))

# Epithelial-cell subclustering
Epithelial<-subset(seu_merge,subset=cell_type=='Tumor cell')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:50,reduction="harmony")
Epithelial <- FindClusters(object = Epithelial,resolution = seq(0.1,1,by=0.1))
clustree(Epithelial@meta.data, prefix = "RNA_snn_res.")
Epithelial <- FindClusters(Epithelial, resolution = 0.3)
Epithelial <- RunUMAP(Epithelial,reduction="harmony", dims = 1:50)
Epithelial <- RunTSNE(Epithelial,reduction="harmony", dims = 1:50)
Epithelial$group<-Epithelial$Type
DimPlot(Epithelial, reduction = "umap", label = TRUE,raster=FALSE,split.by = 'Type')
Epithelial@meta.data[["cell_subtype"]]<-paste0('Epi',Epithelial$seurat_clusters)
# TODO: These levels come from broad cell-type labels and may convert the Epi labels to NA.
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
Idents(Epithelial)<-Epithelial$cell_subtype
DimPlot(Epithelial, reduction = "umap", label = T,raster=FALSE)
saveRDS(Epithelial,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.epi.anno.rds"))


# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type %in% c('Monocyte','pDC','Neutrophil'))
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:50,reduction="harmony")
Mye <- FindClusters(object = Mye,resolution = seq(0.1,1,by=0.1))
clustree(Mye@meta.data, prefix = "RNA_snn_res.")
Mye <- FindClusters(Mye, resolution = 0.4)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:50)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:50)
Mye$group<-Mye$Type
DimPlot(Mye, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'FCN1','VCAN',
          'TOP2A','MKI67')
DotPlot(Mye, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Macrophage','Monocyte','Neutrophil','Macrophage','Monocyte','Macrophage','Unknown','Unknown','pDC')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Mye, reduction = "umap", label = T,raster=FALSE)
saveRDS(Mye,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.myeloid.anno.rds"))

#T
Tcell<-subset(seu_merge,subset=cell_type == 'T cell')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:50,reduction="harmony")
Tcell <- FindClusters(object = Tcell,resolution = seq(0.1,1,by=0.1))
clustree(Tcell@meta.data, prefix = "RNA_snn_res.")
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:50)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:50)
Tcell$group<-Tcell$Type
DimPlot(Tcell, reduction = "umap", label = TRUE,raster=FALSE)

T_marker<-FindAllMarkers(JoinLayers(Tcell),only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
DotPlot(Tcell, features = c('CD4','CD8A','CD8B', 
                            'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7','NCAM1','FCGR3A','FCGR3B',
                            'IL2RA','FOXP3','BATF',
                            'CCR7','TCF7','BACH2','IL7R','FOXP1',   #naive
                            'TNFAIP3','FOS','FTH1','JUNB','CD69',  #Central Memory T (Tcm)
                            'CXCL13','MAF','PDCD1','CD40LG','TOX',    #Tph(CXCL13)
                            'CCL5','CCL4','CST7','PRF1','GZMH',  #CTL
                            'HAVCR2','LAG3','TIGIT',   #Tex
                            'GZMK','AOAH','PITPNC1', #Effector Memory T（Tem）
                            'ZNF683','ITGAE', #Trm
                            'MKI67','TOP2A','BIRC5'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('CD4 Tna','CD8 T','CD4 Treg','CD56+ NK','CD8 T','CD16+ NK','CD8 T','CD4 Th1-like','Proliferating T','CD8 Tex',
                     'Unknown','Unknown','CD56+ NK','Unknown')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Tcell, reduction = "umap", label = T,raster=FALSE)
saveRDS(Tcell,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.T.anno.rds"))


# Export subtype marker genes
load(file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))
Epithelial <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.epi.anno.rds"))
Idents(Epithelial)<-Epithelial@meta.data[["cell_subtype"]]
Epithelial<-JoinLayers(Epithelial)
Epithelial_markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Epithelial_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Dai_Cell.Rep.Med_2024_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Tcell <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.T.anno.rds"))
Idents(Tcell)<-Tcell@meta.data[["cell_subtype"]]
Tcell<-JoinLayers(Tcell)
T_markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Dai_Cell.Rep.Med_2024_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)


Mye <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Dai_Cell.Rep.Med_2024(GSE236738)/GSE236738.myeloid.anno.rds"))
Idents(Mye)<-Mye@meta.data[["cell_subtype"]]
Mye<-JoinLayers(Mye)
Mye_markers <- FindAllMarkers(Mye, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Mye_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Dai_Cell.Rep.Med_2024_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 17: Peng et al., eLife (2025), SRP567748
# Original script: new-data/single-cell/Peng_eLife_2025/SRP567748.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(stringr)
library(clustree)
sample_dirs <- list.dirs(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748"),recursive=F)
seuList <- lapply(paste0(sample_dirs,'/filtered_feature_bc_matrix'),function(x){ 
  seu = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
samplename<-basename(sample_dirs)
group<-c('ADC','ADC','ADC','ADC','ADC','SCC','SCC','ADC','SCC','ADC','ADC','SCC','ADC','ADC','ADC')
names(seuList)<-samplename
VlnPlot(seuList[[1]],features = c('nCount_RNA','nFeature_RNA'))

for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplename[i],ncol(seuList[[i]]))
  Group = rep(group[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["group"]] = Group
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
  sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = TRUE)
  seuList[[i]]<-sce
  print(paste0(samplename[i]," Finished!"))
}

for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}
DimPlot(seuList[[2]], group.by = 'DF.classifications', label = TRUE,raster=FALSE)

seu_merge <- merge(seuList[[1]],
                   y = seuList[-1],
                   add.cell.ids = samplename)
pann_columns <- grep("pANN_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -pann_columns]
classification_columns <- grep("DF.classifications_0.25", colnames(seu_merge@meta.data))
seu_merge@meta.data <- seu_merge@meta.data[, -classification_columns]

seu_merge <- subset(seu_merge, subset = DF.classifications== "Singlet")

seu_merge<-readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/seu_merge.rds"))
VlnPlot(seu_merge,features = c('nCount_RNA','nFeature_RNA'),group.by = 'DF.classifications')
seu_merge <- subset(seu_merge, subset = nFeature_RNA <= 8000 & nCount_RNA<=40000)

seu_merge <- NormalizeData(seu_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seu_merge <- FindVariableFeatures(seu_merge)
seu_merge <- ScaleData(seu_merge)
seu_merge <- RunPCA(seu_merge, features = VariableFeatures(object = seu_merge))
seu_merge <- RunHarmony(seu_merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
seu_merge <- FindNeighbors(seu_merge, dims = 1:30,reduction="harmony")
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "RNA_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.4)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:30)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:30)
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('T cell','Neutrophil','T cell','Epithelial cell','Epithelial cell','Macrophage','Epithelial cell',
                     'Epithelial cell','Plasma cell','Fibroblast','Epithelial cell','Epithelial cell','T cell',
                     'Epithelial cell','Endothelial cell','B cell','Fibroblast','Mast cell','Epithelial cell','T cell')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.sce.anno.rds"))

# Epithelial-cell subclustering
Epithelial<-subset(seu_merge,subset=cell_type=='Epithelial cell')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial)
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Epithelial <- FindNeighbors(Epithelial, dims = 1:30,reduction="pca")
Epithelial <- FindClusters(object = Epithelial,resolution = seq(0.1,1,by=0.1))
clustree(Epithelial@meta.data, prefix = "RNA_snn_res.")
Epithelial <- FindClusters(Epithelial, resolution = 0.3)
Epithelial <- RunUMAP(Epithelial,reduction="pca", dims = 1:30)
Epithelial <- RunTSNE(Epithelial,reduction="pca", dims = 1:30)
DimPlot(Epithelial, reduction = "umap", label = TRUE,raster=FALSE)
marker <- c("NTS", "KRT14", "LY6D", "MT1X", "CCL20",
            "IGLC2", "IGKC", "IGLC3", "IGHG1", "IGHG4",
            "TMPRSS11E", "KRT10", "MMP7", "FTH1", "GPNMB",
            "TFF2", "AGR2", "TFF1", 
            "CCL5", "CCL4", "NKG7", "SRGN", "GZMB",
            "TMC5", "RNF213", "ANKRD36C", "ARGLU1", "MALAT1",
            "CAPS", "TFF3", "CLDN3", "FCGBP", "AGR3",
            "SCGB3A1", "SCGB1D2", "SCGB2A2", "BPIF1", "SCGB2A1",
            "SST", "FXYD4", "C3", "PCSK1", "F3",
            "CYSTM1","MUC5AC","PDK4", "SLC26A3", "CA9", "SPINK1", "LCN2",
            "REG1A", "PLA2G2A", "TCN1", "MEG3", "LYPD2",
            "RRAD", "MT1E", "SLC40A1", "BTG2", "C2orf88",
            'TOP2A','MKI67'
)
DotPlot(Epithelial, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Epi_02_IGLC2','Epi_01_NTS','Epi_03_TMPRSS11E','Epi_05_CCL5','Epi_03_TMPRSS11E','Epi_05_CCL5',
                     'Epi_02_IGLC2','Epi_10_CYSTM1','Epi_01_NTS','Epi_04_TFF2','Epi_06_TMC5','Epi_06_TMC5','Epi_09_SST',
                     'Epi_02_IGLC2','Epi_07_CAPS','Epi_12_RRAD','Epi_06_TMC5','Epi_02_IGLC2','Epi_08_SCGB3A1','Epi_07_CAPS')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids)))
DimPlot(Epithelial, reduction = "umap", label = T,raster=FALSE)
saveRDS(Epithelial,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.epi.anno.rds"))


# Neutrophil subclustering
Neu<-subset(seu_merge,subset=cell_type =='Neutrophil')
Neu@assays[["RNA"]]@meta.data<-data.frame()
Neu <- NormalizeData(Neu, normalization.method = "LogNormalize", scale.factor = 10000)
Neu <- FindVariableFeatures(Neu, selection.method = "vst", nfeatures = 2000)
Neu <- ScaleData(Neu)
Neu <- RunPCA(Neu, features = VariableFeatures(object = Neu))
Neu <- RunHarmony(Neu,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Neu <- FindNeighbors(Neu, dims = 1:30,reduction="harmony")
Neu <- FindClusters(object = Neu,resolution = seq(0.1,1,by=0.1))
clustree(Neu@meta.data, prefix = "RNA_snn_res.")
Neu <- FindClusters(Neu, resolution = 0.5)
Neu <- RunUMAP(Neu,reduction="harmony", dims = 1:30)
Neu <- RunTSNE(Neu,reduction="harmony", dims = 1:30)
DimPlot(Neu, reduction = "umap", label = TRUE,raster=FALSE)
marker<-c('IFITM2','S100A8','LST1','S100A12','CFTM2',
          'CCL5','KRT5','ELF3','JUN','KRT13',
          'CD66','CD11','MT1X','MT2A','MT1G','MT1H','MT1F')
DotPlot(Neu, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Neu_01','Neu_01','Neu_01','Neu_01','Neu_01','Neu_01','Neu_01','Neu_02','Neu_03')
names(new.cluster.ids) <- levels(Neu)
Neu <- RenameIdents(Neu, new.cluster.ids)
Neu@meta.data[["cell_subtype"]]<-Neu@active.ident
Neu@meta.data[["cell_subtype"]]<-factor(Neu@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Neu, reduction = "umap", label = T,raster=FALSE)
saveRDS(Neu,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.neu.anno.rds"))

#T
Tcell<-subset(seu_merge,subset=cell_type == 'T cell')
Tcell@assays[["RNA"]]@meta.data<-data.frame()
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell)
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:30,reduction="harmony")
Tcell <- FindClusters(object = Tcell,resolution = seq(0.1,1,by=0.1))
clustree(Tcell@meta.data, prefix = "RNA_snn_res.")
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:30)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:30)
DimPlot(Tcell, reduction = "umap", label = TRUE,raster=FALSE)
marker<-c('CD8A','TIGIT','PDCD1','LAG3','HAVCR2',
          'CD3D','CD3E',
          'CD4','FOXP3','IL2RA',
          'CD69','CD3G')
DotPlot(Tcell, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Exhausted T','Exhausted T','Treg','Naive T','Cytotoxic T','Tph','Proliferating T','Cytotoxic T',
                     'Exhausted T','Cytotoxic T','Exhausted T','CD16+ NK','Activated T','CD56+ NK')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Tcell, reduction = "umap", label = T,raster=FALSE)
saveRDS(Tcell,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.T.anno.rds"))



# B-cell and plasma-cell subclustering
Plasma<-subset(seu_merge,subset=cell_type %in% c('B cell','Plasma cell'))
Plasma@assays[["RNA"]]@meta.data<-data.frame()
Plasma <- NormalizeData(Plasma, normalization.method = "LogNormalize", scale.factor = 10000)
Plasma <- FindVariableFeatures(Plasma, selection.method = "vst", nfeatures = 2000)
Plasma <- ScaleData(Plasma)
Plasma <- RunPCA(Plasma, features = VariableFeatures(object = Plasma))
Plasma <- RunHarmony(Plasma,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Plasma <- FindNeighbors(Plasma, dims = 1:30,reduction="harmony")
Plasma <- FindClusters(object = Plasma,resolution = seq(0.1,1,by=0.1))
clustree(Plasma@meta.data, prefix = "RNA_snn_res.")
Plasma <- FindClusters(Plasma, resolution = 0.6)
Plasma <- RunUMAP(Plasma,reduction="harmony", dims = 1:30)
Plasma <- RunTSNE(Plasma,reduction="harmony", dims = 1:30)
DimPlot(Plasma, reduction = "umap", label = TRUE,raster=FALSE)
marker <- c("IGHA2", "IGHG2", "IGLC2", "IGHA1", "IGHG3",
            "HLA-DRA", "MS4A1", "HLA-DPB1", "HLA-DRB1","CD37", 
            "KRT17", "FABP5","S100A2", "S100A8", "S100A9",
            "CCL5", "CXCL13", "IL32", "CD3D","CCL4", 
            "CXCL8", "FTH1", "G0S2","C15orf48", "SOD2", 
            "HMGB2","HIST1H1B", "HIST1H4C","MKI67", "TUBA1B")
DotPlot(Plasma, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Plasma/B_02','Plasma/B_01','Plasma/B_01','Plasma/B_01','Plasma/B_01','Plasma/B_04','Plasma/B_02',
                     'Plasma/B_03','Plasma/B_05','Plasma/B_06','Plasma/B_01','Plasma/B_05')
names(new.cluster.ids) <- levels(Plasma)
Plasma <- RenameIdents(Plasma, new.cluster.ids)
Plasma@meta.data[["cell_subtype"]]<-Plasma@active.ident
Plasma@meta.data[["cell_subtype"]]<-factor(Plasma@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids)))
DimPlot(Plasma, reduction = "umap", label = T,raster=FALSE)
saveRDS(Plasma,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.B.anno.rds"))


# Export subtype marker genes
load(file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))
Epithelial <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.epi.anno.rds"))
Idents(Epithelial)<-Epithelial@meta.data[["cell_subtype"]]
Epithelial<-JoinLayers(Epithelial)
Epithelial_markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Epithelial_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Peng_eLife_2025_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Tcell <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.T.anno.rds"))
Idents(Tcell)<-Tcell@meta.data[["cell_subtype"]]
Tcell<-JoinLayers(Tcell)
T_markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Peng_eLife_2025_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Bcell <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.B.anno.rds"))
Idents(Bcell)<-Bcell@meta.data[["cell_subtype"]]
Bcell<-JoinLayers(Bcell)
B_markers <- FindAllMarkers(Bcell, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
B_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Peng_eLife_2025_B_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Neu <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Peng_eLife_2025/SRP567748.neu.anno.rds"))
Idents(Neu)<-Neu@meta.data[["cell_subtype"]]
Neu<-JoinLayers(Neu)
neu_markers <- FindAllMarkers(Neu, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
neu_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Peng_eLife_2025_neu_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

# =============================================================================
# Dataset 18: Cao et al., Journal for ImmunoTherapy of Cancer (2025), treatment cohort
# Original script: new-data/single-cell/Cao_J.Immunother.Cancer_2025/Treatment.R
# =============================================================================

library(Seurat)
library(DoubletFinder)
library(harmony)
library(dplyr)
library(clustree)
samplenames<-list.files(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/matrix文件"))
samplenames<-gsub('_副本','',samplenames)
filename<-paste0(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/matrix文件/"),samplenames,'_副本/',samplenames,'/count/sample_filtered_feature_bc_matrix')
seuList <- lapply(filename,function(x){
  sce = CreateSeuratObject(counts = Read10X(x),
                           min.cells = 3,
                           min.features = 200,
                           assay = "RNA")
})
VlnPlot(seuList[[2]],features = c('nCount_RNA','nFeature_RNA'))
names(seuList)<-samplenames
Group<-rep(c('Pre','Post'),times=5)
Response<-c(rep('NMPR',8),rep('MPR',2))
for(i in 1:length(seuList)){
  sce<-seuList[[i]]
  orig = rep(samplenames[i],ncol(seuList[[i]]))
  group = rep(Group[i],ncol(seuList[[i]]))
  response = rep(Response[i],ncol(seuList[[i]]))
  sce[["orig.ident"]] = orig
  sce[["group"]] = group
  sce[["Response"]] = response
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  sce <- subset(sce, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & nCount_RNA <= 45000 & percent.mt <= 25)
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
  sce <- ScaleData(sce, vars.to.regress = "percent.mt")
  #sce <- SCTransform(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  sce <- RunPCA(sce, features = VariableFeatures(object = sce))
  sce <- FindNeighbors(sce, dims = 1:20,reduction="pca")
  sce <- FindClusters(sce, resolution = 1)
  sce <- RunUMAP(sce,reduction="pca", dims = 1:20)
  # TODO: The parameter sweep uses PCs 1:15, while doubletFinder below uses PCs 1:20.
sweep.res.list <- paramSweep(sce, PCs = 1:15, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal),]
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
                       nExp = nExp_poi.adj, sct = TRUE)
  seuList[[i]]<-sce
  print(paste0(samplenames[i]," Finished!"))
}
for(i in 1:length(seuList)){
  seuList[[i]][["pANN"]]<-seuList[[i]]@meta.data %>% select(contains('pANN'))
  seuList[[i]][["DF.classifications"]]<-seuList[[i]]@meta.data %>% select(contains('DF.classifications'))
}

seu_merge <- merge(seuList[[1]],
                    y = seuList[-1],
                    add.cell.ids = samplenames)
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
seu_merge <- FindClusters(object = seu_merge,resolution = seq(0.1,1,by=0.1))
clustree(seu_merge@meta.data, prefix = "RNA_snn_res.")
seu_merge <- FindClusters(seu_merge, resolution = 0.1)
seu_merge <- RunUMAP(seu_merge,reduction="harmony", dims = 1:15)
seu_merge <- RunTSNE(seu_merge,reduction="harmony", dims = 1:15)
DimPlot(seu_merge, reduction = "umap", label = TRUE,raster=FALSE)

marker<-c('COL1A1','MMP11','DCN','COL6A3','SFRP4','COL1A2','COL12A1','LUM',          #fibroblasts
          'CD2', 'CD3D','CD3E','CD4','CD8A','CD8B', 'IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7',        #NK
          'ITGAX', 'CSF1R', 'FCGR3A',          #myeloid cells
          'FCGR2A','S100A8','S100A9',                     #FCGR2A+monocytes
          'CD14', 'CD68','CD163', 'MS4A6A','C1QB','C1QA',              #macrophages
          'CD1C','LAMP3', 'PLD4',              #dendritic cells (DCs)
          'FCGR3B','PTGS2','CSF3R',
          "CLEC4C","LILRA4",'GPR183',
          'CDKN2A', 'CDH1', 'EPCAM','KRT14','KRT5','KRT6A','WFDC2','TSPAN8','KRT8','KRT18','MUC5B','KRT19',    #epithelial cells
          'CHGA', 'CHGB', 'SYP', 'NCAM1','NRXN1', 'INSM1','ASCL1','ASCL2',
          'ACTA2', 'RGS5','MYH11',              #smooth muscle cells
          'CLDN5','CDH5', 'EMCN','PECAM1', 'PCDH17','KDR','A2M','VWF','ENG','RAMP2',             #ECs
          'JCHAIN','CD38', 'IGLC2','IGHG4',                     #plasma cells
          'MS4A1','CD79B','IGKC','CD79A','MZB1',                     #B cells
          'CPA3','TPSAB1','KIT', 'IL1RL1','MS4A2',                     #mast cells
          'HBB','HBA1','HBA2','TOP2A','MKI67'
)
DotPlot(seu_merge, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Epithelial','T/NK','Epithelial','Fibroblast','Myeloid','B/Plasma','B/Plasma','Endothelial','Fibroblast','Mast',
                     'T/NK','Epithelial','Epithelial','Epithelial')
names(new.cluster.ids) <- levels(seu_merge)
seu_merge <- RenameIdents(seu_merge, new.cluster.ids)
seu_merge@meta.data[["cell_type"]]<-seu_merge@active.ident
seu_merge@meta.data[["cell_type"]]<-factor(seu_merge@meta.data[["cell_type"]],levels=unique(new.cluster.ids))
DimPlot(seu_merge, reduction = "umap", label = T,raster=FALSE)
saveRDS(seu_merge,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.sce.anno.rds"))


# Epithelial-cell subclustering
Epithelial <- subset(seu_merge,subset=cell_type=='Epithelial')
Epithelial <- NormalizeData(Epithelial, normalization.method = "LogNormalize", scale.factor = 10000)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, vars.to.regress = "percent.mt")
Epithelial <- RunPCA(Epithelial, features = VariableFeatures(object = Epithelial))
Epithelial <- RunHarmony(Epithelial,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
# NOTE: Harmony is calculated above, but the original downstream epithelial workflow uses PCA.
Epithelial <- FindNeighbors(Epithelial, dims = 1:15,reduction="pca")
Epithelial <- FindClusters(object = Epithelial,resolution = seq(0.1,1,by=0.1))
clustree(Epithelial@meta.data, prefix = "RNA_snn_res.")
Epithelial <- FindClusters(Epithelial, resolution = 0.2)
Epithelial <- RunUMAP(Epithelial,reduction="pca", dims = 1:15)
Epithelial <- RunTSNE(Epithelial,reduction="pca", dims = 1:15)
DimPlot(Epithelial, reduction = "umap", label = TRUE,raster=FALSE)
DotPlot(Epithelial, features = c('SPP1','EPCAM','MYC','CEACAM5','MUC5B','MKI67','NEAT1','TP63','SERPINA1'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('EP0_SPP1','EP8_SERPINA1','EP0_SPP1','EP3_CEACAM5','Unknown','EP1_EPCAM','EP0_SPP1','EP6_NEAT1',
                     'EP3_CEACAM5','EP2_MYC','EP3_CEACAM5','EP4_MUC5B','EP0_SPP1','EP0_SPP1','EP0_SPP1','EP6_NEAT1')
names(new.cluster.ids) <- levels(Epithelial)
Epithelial <- RenameIdents(Epithelial, new.cluster.ids)
Epithelial@meta.data[["cell_subtype"]]<-Epithelial@active.ident
Epithelial@meta.data[["cell_subtype"]]<-factor(Epithelial@meta.data[["cell_subtype"]],levels=sort(unique(new.cluster.ids)))
DimPlot(Epithelial, reduction = "umap", label = T,raster=FALSE)
saveRDS(Epithelial,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.epi.anno.rds"))


# Myeloid-cell subclustering
Mye<-subset(seu_merge,subset=cell_type =='Myeloid')
Mye <- NormalizeData(Mye, normalization.method = "LogNormalize", scale.factor = 10000)
Mye <- FindVariableFeatures(Mye, selection.method = "vst", nfeatures = 2000)
Mye <- ScaleData(Mye, vars.to.regress = "percent.mt")
Mye <- RunPCA(Mye, features = VariableFeatures(object = Mye))
Mye <- RunHarmony(Mye,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Mye <- FindNeighbors(Mye, dims = 1:15,reduction="harmony")
Mye <- FindClusters(object = Mye,resolution = seq(0.1,1,by=0.1))
clustree(Mye@meta.data, prefix = "RNA_snn_res.")
Mye <- FindClusters(Mye, resolution = 0.4)
Mye <- RunUMAP(Mye,reduction="harmony", dims = 1:15)
Mye <- RunTSNE(Mye,reduction="harmony", dims = 1:15)
DimPlot(Mye, reduction = "umap", label = TRUE,raster=FALSE)
marker <- c("IL1B", "FCN1", "VCAN",
            "SPP1", "TREM2", "MMP9",
            "CCL18", "FOLR2", "SLC40A1",
            "C1QC", "C1QB", "C1QA",
            "CD1E", "CLEC10A", "CD1C",
            "LAMP3", "CCR7", "FSCN1",
            "GZMB", "JCHAIN", "LILRA4","CLEC4C",'GPR183',
            'FCGR3B','PTGS2','CSF3R','CXCL8','TOP2A','MKI67')
DotPlot(Mye, features = marker)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('Macro_CCL18','Macro_IL1B','Macro_IL1B','Macro_SPP1','Unknown','Macro_C1QC','Unknown','Macro_SPP1','cDC_CD1C',
                     'pDC_GZMB','Macro_MKI67','Macro_IL1B','Macro_SPP1','mDC_LAMP3','Macro_SPP1')
names(new.cluster.ids) <- levels(Mye)
Mye <- RenameIdents(Mye, new.cluster.ids)
Mye@meta.data[["cell_subtype"]]<-Mye@active.ident
Mye@meta.data[["cell_subtype"]]<-factor(Mye@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Mye, reduction = "umap", label = T,raster=FALSE)
saveRDS(Mye,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.myeloid.anno.rds"))

#T
Tcell<-subset(seu_merge,subset=cell_type =='T/NK')
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell, vars.to.regress = "percent.mt")
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
Tcell <- RunHarmony(Tcell,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
Tcell <- FindNeighbors(Tcell, dims = 1:15,reduction="harmony")
Tcell <- FindClusters(object = Tcell,resolution = seq(0.1,1,by=0.1))
clustree(Tcell@meta.data, prefix = "RNA_snn_res.")
Tcell <- FindClusters(Tcell, resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="harmony", dims = 1:15)
Tcell <- RunTSNE(Tcell,reduction="harmony", dims = 1:15)
DimPlot(Tcell, reduction = "umap", label = TRUE,raster=FALSE)
marker<-c('CD4','CD8A','CD8B','TOP2A','MKI67','IL2RA','FOXP3','BATF',                          #T cells
          'GNLY','XCL1','KLRB1', 'NCR1', 'FGFBP2','KLRC1','KLRF1','NKG7','NCAM1','FCGR3A','FCGR3B')
DotPlot(Tcell, features = c('TCF7','ZNF683','PDCD1','CTLA4','KLRC2'))+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 45,hjust = 1,vjust=1))
new.cluster.ids <- c('CD4_Tcm','CD4_Treg','CD8_CTL','CD4_Tn','CD8_Tex_PDCD1','CD8_Trm_ZNF683','CD16+_NK','CD8_Tn',
                     'CD56+_NK','Proliferating_T','CD8_Tn','Proliferating_T',
                     'CD16+_NK','Proliferating_T','Unknown','CD4_Treg')
names(new.cluster.ids) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, new.cluster.ids)
Tcell@meta.data[["cell_subtype"]]<-Tcell@active.ident
Tcell@meta.data[["cell_subtype"]]<-factor(Tcell@meta.data[["cell_subtype"]],levels=unique(new.cluster.ids))
DimPlot(Tcell, reduction = "umap", label = T,raster=FALSE)
saveRDS(Tcell,file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.T.anno.rds"))

# Export subtype marker genes
load(file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))
Epithelial <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.epi.anno.rds"))
Idents(Epithelial)<-Epithelial@meta.data[["cell_subtype"]]
Epithelial<-JoinLayers(Epithelial)
Epithelial_markers <- FindAllMarkers(Epithelial, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Epithelial_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_J.Immunother.Cancer_2025_epi_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

Tcell <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.T.anno.rds"))
Idents(Tcell)<-Tcell@meta.data[["cell_subtype"]]
Tcell<-JoinLayers(Tcell)
T_markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
T_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_J.Immunother.Cancer_2025_T_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)


Mye <- readRDS(file.path(project_dir, "返修/新纳入数据/单细胞/数据集/Cao_J.Immunother.Cancer_2025/Treatment.myeloid.anno.rds"))
Idents(Mye)<-Mye@meta.data[["cell_subtype"]]
Mye<-JoinLayers(Mye)
Mye_markers <- FindAllMarkers(Mye, only.pos = TRUE, min.pct = 0.25, recorrect_umi = FALSE, logfc.threshold = 0.25)
Mye_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) -> deg_top
deg_top <- split(deg_top$gene, deg_top$cluster)
names(deg_top)<-paste0('Cao_J.Immunother.Cancer_2025_myeloid_',names(deg_top))
subcelltype_marker_list<-c(subcelltype_marker_list,deg_top)

save(subcelltype_marker_list,file=file.path(project_dir, "返修/新纳入数据/单细胞/subcelltype_marker_list.RData"))

