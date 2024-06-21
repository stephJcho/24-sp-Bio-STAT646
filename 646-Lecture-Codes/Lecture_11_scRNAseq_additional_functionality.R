getwd()

library(Seurat)
library(ggplot2)
library(tidyverse)

#############################################################
# Read in a data matrix from a public repository: example 1
#############################################################
# Mouse Cell Atlas 2.0 adult mouse liver:
# https://bis.zju.edu.cn/MCA/dpline.html?tissue=Liver

count=read.csv('data/Adult-Liver_dge.csv', row.names = 1)
#can just use read.csv cuz it's already been processed
dim(count)
#dim() for sanity check, and also could check if data is not corrupt
count[1:5,1:5]
#count: cell x genes matrix
annotation=read.csv('data/Adult-Liver_barcodes_anno.csv', row.names = 1)
head(annotation)
#check if the count's col names and annot's row names match
colnames(annotation)=c('celltype','cluster')
all(colnames(count)==rownames(annotation))
genes=read.csv('data/Adult-Liver_gene.csv', row.names = 1)
rownames(count)=genes[,1]
rownames(count)
rm(genes)

mca <- CreateSeuratObject(counts = count, project = "mca", min.cells = 3, min.features = 200,
                          meta.data = annotation)
#last two arguments are very crude initial value, will do qc later
#rest of the steps are similar to pbmc3k example, you do stuff with pbmc object
rm(count); rm(annotation)

mca[["percent.mt"]] <- PercentageFeatureSet(mca, pattern = "^mt-")
# calculating the mitochondria percentage
VlnPlot(mca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
mca <- subset(mca, subset = nFeature_RNA > 300 & nFeature_RNA < 1000 &
                nCount_RNA<1500 & percent.mt < 5)
# subsetting is to remove outlier cells (a qc step)
VlnPlot(mca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# can automatically set threshold by setting range from median
mca <- NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca <- FindVariableFeatures(mca, selection.method = "vst", nfeatures = 2000)
mca <- ScaleData(mca, features = rownames(mca))
mca <- RunPCA(mca, features = VariableFeatures(object = mca))
ElbowPlot(mca)
mca <- RunTSNE(mca, dims = 1:10)

#p1: used tsne
p1=DimPlot(mca, reduction = "tsne", group.by = 'celltype', label = TRUE)
#group.by to use group from original data as a reference
p1

# Now, generate a visualization of gene expression
# But, unlike pbmc3k (based on clusters), it should now be based on cell types

Idents(mca)=mca$celltype
#to decide btw cell type specific expression or cluster specific exp. 
mca_markers <- FindAllMarkers(mca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mca_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(mca, features = top10$gene) + NoLegend()

genes.to.plot=top10$gene[seq(1,nrow(top10),20)]
FeaturePlot(mca, features = genes.to.plot, ncol = 3)

#############################################################
# Read in a data matrix from a public repository: example 2
#############################################################
# # Ignore the following script, which is used to generate the example.
# reference <- SeuratDisk::LoadH5Seurat("data/pbmc_multimodal_citeseq.h5seurat")
# table(reference$donor)
# reference = reference[,reference$donor %in% paste0('P',1)]
# table(reference$celltype.l1)
# table(reference$celltype.l2)
# celltype=reference@meta.data[,c('celltype.l1','celltype.l2')]
# count=reference@assays$SCT@data
# pbmc.ref <- CreateSeuratObject(counts = count, project = "pbmc.ref",
#                           min.cells = 3, min.features = 200,
#                           meta.data = celltype)
# saveRDS(pbmc.ref, file='data/pbmc.ref.rds')

# Load data
rm(list = ls()) # Remove all data
gc() # Clear RAM
pbmc.3k=readRDS('data/pbmc3k_final.rds') # From Lecture_9_PBMC_3K_Tutorial.Rmd
pbmc.ref=readRDS('data/pbmc.ref.rds') # Reference PBMC data with annotated cell types.

# Separate data analysis
pbmc.3k <- NormalizeData(pbmc.3k, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.3k <- FindVariableFeatures(pbmc.3k, selection.method = "vst", nfeatures = 2000)
pbmc.3k <- ScaleData(pbmc.3k, features = rownames(pbmc.3k))
pbmc.3k <- RunPCA(pbmc.3k, features = VariableFeatures(object = pbmc.3k))
ElbowPlot(pbmc.3k)
pbmc.3k <- RunTSNE(pbmc.3k, dims = 1:10)

# now we work with two separate datasets

#details of dataset
table(pbmc.ref$celltype.l1)
table(pbmc.ref$celltype.l2)
# assigning the markers (or labeling) in pbmc3k was manual
# usually, we can find the reference and just import it, just like this example

# no qc step at all in this file; advised to do it during midterm

pbmc.ref <- NormalizeData(pbmc.ref, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc.ref <- FindVariableFeatures(pbmc.ref, selection.method = "vst", nfeatures = 2000)
pbmc.ref <- ScaleData(pbmc.ref, features = rownames(pbmc.ref))
pbmc.ref <- RunPCA(pbmc.ref, features = VariableFeatures(object = pbmc.ref))
ElbowPlot(pbmc.ref)
pbmc.ref <- RunTSNE(pbmc.ref, dims = 1:10)

p1=DimPlot(pbmc.ref, group.by='celltype.l1', label = TRUE)
p2=DimPlot(pbmc.ref, group.by='celltype.l2', label = TRUE)
p1 + p2

# Cell type label transfer
# find anchors between reference and query
anchors <- FindTransferAnchors(
  reference = pbmc.ref,
  query = pbmc.3k,
  dims = 1:10,
  reference.reduction = "pca",
  recompute.residuals = FALSE
)

# Transfer cell type labels and protein data from the reference to the query
# project the query data onto the umap structure of the reference
pbmc.3k <- MapQuery(
  anchorset = anchors,
  query = pbmc.3k,
  reference = pbmc.ref,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2"
  ),
  reference.reduction = "pca"
)

p3=DimPlot(pbmc.3k, group.by='predicted.celltype.l1', label = TRUE)
p4=DimPlot(pbmc.3k, group.by='predicted.celltype.l2', label = TRUE)
p3 + p4
# same result as pbmc3k file, but the labeling is done automatically

# not mandatory for exam, but does help
par(mfrow=c(1,2))
hist(pbmc.3k$predicted.celltype.l1.score) # Sometimes you may want to apply cell QC to this
hist(pbmc.3k$predicted.celltype.l2.score)
par(mfrow=c(1,1))

# Data integration
# https://satijalab.org/seurat/archive/v3.0/pbmc_alignment.html
pbmc.3k$datasource='3k'
pbmc.ref$datasource='ref'
pbmc.3k$celltype.l1=pbmc.3k$predicted.celltype.l1
pbmc.3k$celltype.l2=pbmc.3k$predicted.celltype.l2

pbmc.anchors <- FindIntegrationAnchors(object.list = list(pbmc.3k, pbmc.ref), dims = 1:10)
pbmc.combined <- IntegrateData(anchorset = pbmc.anchors, dims = 1:10)
pbmc.combined

DefaultAssay(pbmc.combined) <- "RNA" # Without integration (i.e., concatenating cells)
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)

p5 <- DimPlot(pbmc.combined, reduction = "umap", group.by = "datasource", label=TRUE)
p6 <- DimPlot(pbmc.combined, reduction = "umap", group.by='celltype.l1', label=TRUE)
p5+p6

DefaultAssay(pbmc.combined) <- "integrated" # With CCA integration
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindNeighbors(pbmc.combined, reduction = "pca", dims = 1:20)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.5)

p7 <- DimPlot(pbmc.combined, reduction = "umap", group.by = "datasource", label=TRUE)
p8 <- DimPlot(pbmc.combined, reduction = "umap", group.by='celltype.l1', label=TRUE)
p7+p8

