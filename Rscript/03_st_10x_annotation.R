# Purpose: Integrate, cluster, and manually annotate the selected 10x spatial transcriptomics samples.
# Input: preprocessed spatial Seurat object with an SCT assay.
# Output: annotated spatial Seurat object.
# NOTE: The selected samples, clustering parameters, and positional annotation vectors are preserved.

local_project_dir <- Sys.getenv("ICCOR_LOCAL_PROJECT_DIR", unset = "ST")

library(Seurat)
library(harmony)
library(clustree)
library(ggplot2)
library(scales)
dir<-list.dirs('spaceranger')
names(dir)<-basename(dir)
CC <- list()
for(i in 1:length(dir)){
  CC[[i]] <-Seurat::Load10X_Spatial(data.dir = dir[i])
  CC[[i]]@meta.data$orig.ident <-names(dir)[i]
  CC[[i]] <- SCTransform(CC[[i]], assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
}
CC.merge <- merge(CC[[1]], y=CC[-1],add.cell.ids = names(dir))
CC.merge <- RunPCA(CC.merge, assay = "SCT", verbose = FALSE)
CC.merge <- RunHarmony(CC.merge,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
CC.merge <- FindNeighbors(CC.merge, reduction = "harmony", dims = 1:20)
CC.merge <- FindClusters(CC.merge, verbose = FALSE,resolution = 0.5)
CC.merge <- RunUMAP(CC.merge, reduction = "harmony", dims = 1:20)
CC.merge <- RunTSNE(CC.merge, reduction = "harmony", dims = 1:20)
new.cluster.ids <- c('Stromal','B/Plasma','Epithelial','Epithelial','Stromal','Epithelial','B/Plasma','Stromal',
                     'Epithelial','Epithelial','Epithelial','Endothelial','Unknown','Epithelial','Epithelial','Macrophages')
names(new.cluster.ids) <- levels(CC.merge)
CC.merge <- RenameIdents(CC.merge, new.cluster.ids)
CC.merge@meta.data[["cell_type"]]<-CC.merge@active.ident
DimPlot(CC.merge, reduction = "umap", label = TRUE,raster=FALSE)
saveRDS(CC.merge, file.path(local_project_dir, "ST_10X.anno.rds"))





