
###############################################
# p-values under the null
###############################################

# null distribution of independent p-values 
pval=rep(NA,10000)
for(i in 1:length(pval)){
  x1=rnorm(10,0,1)
  x2=rnorm(10,0,1)
  pval[i]=t.test(x1,x2)$p.value
}
hist(pval,100)

###############################################
# p-values from parametric vs non-parametric testing
###############################################

# Observed data
set.seed(1)
x1=rnorm(10,-1,1)
x2=rnorm(10,1,1)
boxplot(x1,x2)

qqnorm(x1); qqline(x1)
qqnorm(x2); qqline(x2)

t.test(x1, x2)$p.value #parametric
wilcox.test(x1, x2)$p.value #nonparametric

test.stat=function(x1,x2){
  # t.stat=(mean(x1)-mean(x2))/sqrt(var(x1)/10+var(x2)/10) # t-statistic
  t.stat=mean(x1)-mean(x2)
  return(t.stat)
}
t.stat.obs=test.stat(x1,x2)

# permutation-based
n=10000
t.stat.samp=rep(NA,10000)
t.stat.samp[1]=t.stat.obs # usually we put the observed statistics as the first sampled one
set.seed(2)
for(i in 2:n){
  if(i %% 1000 ==0) cat(i,'\t')
  temp=c(x1,x2)[sample(1:20)]
  t.stat.samp[i]=test.stat(temp[1:10],temp[11:20])
}
#1st 10: 1st group; 2nd 10: 2nd group
hist(t.stat.samp,100)
abline(v=t.stat.obs,col=2,lty=2) #the stitched line; from unpermutated sample

sum(abs(t.stat.samp)>=abs(t.stat.obs))/n

# bootstrap-based
n=10000
t.stat.samp=rep(NA,10000)
t.stat.samp[1]=t.stat.obs
set.seed(2)
for(i in 2:n){
  if(i %% 1000 ==0) cat(i,'\t')
  temp=c(x1,x2)[sample(1:20, replace=T)]
  t.stat.samp[i]=test.stat(temp[1:10],temp[11:20])
}
hist(t.stat.samp,100)
abline(v=t.stat.obs,col=2,lty=2)

sum(abs(t.stat.samp)>=abs(t.stat.obs))/n


###############################################
# Differential expression
###############################################
# http://bioconductor.org/packages/release/data/experiment/html/pasilla.html
# Load in dataset
library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]

head(cts,2)
coldata

rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# Construct DESeqDataSet
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
mcols(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]

# Run differential expression
dds <- DESeq(dds)
res <- results(dds)
res


################################################
# Explore p-values and different multiple testint correction methods
################################################

hist(res$pvalue,100) # nominal p-value
hist(res$padj,100) # adjusted p-value, returned by DESeq2

alpha=0.05 # significance level
pval=res$pvalue
sum(is.na(pval))
pval=pval[!is.na(pval)]
m=length(pval) # number of testing (i.e., genes after QC)

# nominal (PCER) (No correction)
sum(pval < alpha)

# Bonferroni (FWER)
sum(pval < alpha/m)
sum(pval*m < alpha)

pval.bonferroni1=pmin(1,pval*m)
pval.bonferroni2=p.adjust(p = pval, method='bonferroni')

head(cbind(pval.bonferroni1, pval.bonferroni2))
tail(cbind(pval.bonferroni1, pval.bonferroni2))

plot(pval.bonferroni1, pval.bonferroni2)
max(abs(pval.bonferroni1-pval.bonferroni2))

# Benjamini-Hochberg (FDR): my implementation
plot(sort(pval),alpha*(1:m)/m, type='l')
abline(a=0, b=1, col=2)

plot(sort(pval),alpha*(1:m)/m, type='l', xlim=c(0,0.01), ylim=c(0,0.01))
abline(a=0, b=1, col=2)

max(which(alpha*(1:m)/m-sort(pval)>=0))
padj=(sort(pval)*m/(1:m))[rank(pval)]

padj2=p.adjust(pval, method='BH')
plot(padj, padj2); abline(a=0,b=1, col='red', lty=2)
max(abs(padj-padj2))


# multtest package: BH FDR control
library(multtest)
padj2=mt.rawp2adjp(pval,"BH")$adjp
head(padj2)

padj2=padj2[rank(pval),2]
sum(padj2<alpha)

max(abs(padj-padj2)) 
plot(padj, padj2); abline(a=0, b=1, col=2)

# q-value developed by John Storey
library(qvalue)
qval=qvalue(p=pval)$qvalues
sum(qval < alpha)
plot(pval, qval)
plot(padj, qval); abline(a=0,b=1, col=2)
# q-value returns more significant calls while also controlling for FDR



