## Disclaimer: This file has been modified from the original R file (provided by Dr. Yuchao Jiang);
## such as notes taken from the lecture, annotations, extra links and omitting of redundant lines of code.

data(iris)
X=as.matrix(iris[,1:4])
X2=apply(X,2,scale)
dim(X)

#############################
# PCA
#############################
par(mfrow=c(2,2))

# prcomp
pca1=prcomp(X)
plot(pca1$x[,1], pca1$x[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
#legend('topright', col=c(1,2,3), legend=unique(iris$Species), pch=c(16,16,16))
title('prcomp')

# prcomp without scaling
pca2=prcomp(X, center = F, scale=F)
plot(pca2$x[,1], pca2$x[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
legend('topleft', col=c(1,2,3), legend=unique(iris$Species), pch=c(16,16,16))
title('prcomp no scaling')

# SVD
pca3=X2%*%(svd(X2)$v)
plot(pca3[,1], pca3[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
#legend('bottomright', col=c(1,2,3), legend=unique(iris$Species), pch=c(16,16,16))
title('SVD')

# SVD with scaling
pca4=X%*%(svd(X)$v)
plot(pca4[,1], pca4[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
legend('topleft', col=c(1,2,3), legend=unique(iris$Species), pch=c(16,16,16))
title('SVD no scaling')

#############################
# Efficient PCA
#############################
library(Seurat)
pbmc=readRDS('pbmc3k_final.rds')
pbmc=UpdateSeuratObject(pbmc)
dim(pbmc@assays$RNA@scale.data)
X.sc=t(pbmc@assays$RNA@scale.data)
dim(X.sc) # cell by gene scaled matrix
# pbmc.pc1=X.sc%*%(svd(X.sc, nu=10, nv=10)$v) # This will take VERY long
library(irlba) # faster implementation of PCA
pbmc.pc2=prcomp_irlba(X.sc, n=10)
par(mfrow=c(1,1))
plot(pbmc.pc2$x[,1],pbmc.pc2$x[,2], col=pbmc$seurat_clusters, pch=16,
     xlab='PC1', ylab='PC2')

#############################
# PCA, sparse PCA, and kernal PCA
#############################
par(mfrow=c(1,3))
# prcomp
pca1=prcomp(X)
pca1$rotation
plot(pca1$x[,1], pca1$x[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
#legend('topright', col=c(1,2,3), legend=unique(iris$Species), pch=c(16,16,16))
title('PCA')

library(sparsepca)
spca=robspca(X2, alpha=0.01)
spca$loadings
plot(spca$scores[,1],spca$scores[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
title('sparse PCA')

library(kernlab)
kpca=kpca(X2)
plot(kpca@pcv[,1],kpca@pcv[,2], col=iris$Species, pch=16, xlab='PC1', ylab='PC2')
title('kernel PCA')
# nonlinear projection on the original space

#############################
# tSNE
#############################
library(tsne)
par(mfrow=c(1,3))
tsne_data=tsne(X,k=3,initial_dims=30)
plot(tsne_data[,1],tsne_data[,2],col=iris$Species, pch=16, xlab='tSNE1', ylab='tSNE2')
title('tSNE')

#############################
# UMAP
#############################
library(umap)
iris.umap = umap(X)
plot(iris.umap$layout[,1], iris.umap$layout[,2], pch=16, xlab='UMAP1', ylab='UMAP2', col=iris$Species)
title('UMAP')

#############################
# Python UMAP
#############################
library(reticulate)
reticulate::use_python('/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11', required = TRUE)
iris.umap_learn <- umap(X, method="umap-learn")
plot(iris.umap_learn$layout[,1], iris.umap_learn$layout[,2], pch=16, xlab='UMAP1', ylab='UMAP2', col=iris$Species)
title('UMAP learn')

a=Sys.time()
pbmc.umap1=umap(as.matrix(pbmc.pc2$x))
Sys.time()-a

a=Sys.time()
pbmc.umap2=umap(as.matrix(pbmc.pc2$x), method="umap-learn")
Sys.time()-a

par(mfrow=c(1,2))
plot(pbmc.pc2$x[,1],pbmc.pc2$x[,2], col=pbmc$seurat_clusters, 
     pch=16, cex=0.5, xlab='PC1', ylab='PC2', main='PCA')
plot(pbmc.umap2$layout[,1], pbmc.umap2$layout[,2], col=pbmc$seurat_clusters, 
     pch=16, cex=0.5, xlab='UMAP1', ylab='UMAP2', main='UMAP')

