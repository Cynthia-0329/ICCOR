# Purpose: Integrate GEO, WPP, and TCGA-CESC bulk RNA-seq counts, correct batch effects, and normalize expression.
# Input: featureCounts tables, WPP counts, TCGA-CESC STAR counts, and sample batch metadata.
# Output: combined normalized expression matrices, dataset-specific matrices

project_dir <- Sys.getenv("ICCOR_PROJECT_DIR")
dataset_output_dir <- getwd()

library(dplyr)
library(openxlsx)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sva)
library(DESeq2)
library(RColorBrewer)
library(tibble)
# GEO counts
gene_count_merge<-data.frame()
for (i in list.files(file.path(project_dir, "bulkRNAseq"))){
  data_set<-strsplit(i,'\\.')[[1]][1]
  gene_count<-read.table(paste0(file.path(project_dir, "bulkRNAseq"),i),fill=TRUE,header=T)
  rownames(gene_count)<-gene_count$Geneid
  gene_count_merge<-cbind(gene_count_merge,gene_count[,7:length(gene_count)])
}
colnames(gene_count_merge)[2:length(gene_count_merge)]<-gsub('..Aligned_Bam.','',colnames(gene_count_merge)[2:length(gene_count_merge)]) %>%
  gsub('SE.','',x=.) %>%
  gsub('PE.','',x=.) %>%
  gsub('.sorted.bam','',x=.)
gene_count_merge$Geneid<-NULL

# TCGA-CESC counts
query <- GDCquery(
  project = "TCGA-CESC", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
cesc_se <- GDCprepare(query)
TCGA_count <- assay(cesc_se, "unstranded")
TCGA_count <- as.data.frame(TCGA_count[!grepl('_PAR_Y', rownames(TCGA_count)),])
rownames(TCGA_count)<-unlist(lapply(strsplit(rownames(TCGA_count),'\\.'),function(x){x[[1]]}))
colnames(TCGA_count)<-lapply(colnames(TCGA_count),function(x){substr(x,1,16)})

gene_id<-intersect(intersect(rownames(gene_count_merge),rownames(WPP_count)),rownames(TCGA_count))
gene_count_merge<-gene_count_merge[gene_id,]
TCGA_count<-TCGA_count[gene_id,]

gene_count_merge<-cbind(TCGA_count,gene_count_merge)
batch_info<-read.xlsx(file.path("Sample_Batch.xlsx"))
gene_count_merge<-gene_count_merge[,batch_info$Sample]

# Correct batch effects across cohorts
gene_count_combat <- ComBat_seq(counts = as.matrix(gene_count_merge),batch = batch_info$Batch)

sample_info <- data.frame(
  batch = batch_info$Batch,
  condition = batch_info$Type,
  row.names = colnames(gene_count_combat)
)
dds <- DESeqDataSetFromMatrix(countData = gene_count_combat,
                              colData = sample_info,
                              design = ~ 1)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts<-as.data.frame(normalized_counts) %>% rownames_to_column('gene_name')
write.table(normalized_counts,file=file.path(project_dir, "Bulk_RNAseq_combat_norm_expr.tsv"),
            col.names = T,sep='\t',quote = F,row.names = F)