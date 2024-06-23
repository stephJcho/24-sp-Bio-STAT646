###########################################################################
# Paper: PMID 27667667
# Data link (available in paper):
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061
###########################################################################

# From E-MTAB-5061.idf.txt
# Format of the datafile ‘pancreas_refseq_rpkms_counts_3514sc.txt’: 
# The file contains both the normalized rpkm values and the raw read counts 
# for each sample. Columns correspond to samples and rows to genes. 
# The first line of the file (starting with: #samples) contains the sample IDs 
# to be used as column labels for both the rpkm and counts.
# The columns of the rpkm and the counts have the same order with the sample IDs.
# Columns 1:3514 correspond to rpkm values, Columns 3515:7028 correspond to read counts.
# Rows 1:26179 correspond to data for RefSeq genes, Rows 26180:26271 correspond 
# to data for the 92 external RNA spike-in controls (ERCCs), 
# Row 26272 (last) contains data for ‘eGFP’.

cell.barcodes <- scan(text = readLines("data/E-MTAB-5061/pancreas_refseq_rpkms_counts_3514sc.txt", 1), 
                      what = "", quiet = TRUE)[-1] # 3514 cell barcodes
count=read.table('data/E-MTAB-5061/pancreas_refseq_rpkms_counts_3514sc.txt',sep='\t')
gene.meta=count[,1:2] # First two columns are gene symbols and IDs
count=count[,-(1:2)]
count=as.matrix(count[,3515:7028]) # Columns 3515:7028 correspond to read counts.
rownames(count)=gene.meta[,1]
colnames(count)=cell.barcodes
rm(gene.meta, cell.barcodes)

cell.meta=read.csv('data/E-MTAB-5061/E-MTAB-5061.sdrf.txt',sep='\t')
table(cell.meta$Characteristics..individual.) # Six healthy individuals, four T2D patients
cell.meta=cell.meta[,1:5] # Remove unnecessary meta info for the cells

# Hint:
# 1. Need to reorder the cells to make them match between count and cell.meta
# 2. Need to select the healthy individuals (and remove the T2D patients)

###########################################################################
# Paper: PMID 27667667
# Data link (available in paper):
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061
###########################################################################

baron1=read.csv('data/GSE84133_RAW/GSM2230757_human1_umifm_counts.csv.gz')
baron2=read.csv('data/GSE84133_RAW/GSM2230758_human2_umifm_counts.csv.gz')
baron3=read.csv('data/GSE84133_RAW/GSM2230759_human3_umifm_counts.csv.gz')
baron4=read.csv('data/GSE84133_RAW/GSM2230760_human4_umifm_counts.csv.gz')

# Check that the columns are the same before combining rows
all(colnames(baron1)==colnames(baron2))
all(colnames(baron1)==colnames(baron3))
all(colnames(baron1)==colnames(baron4))

baron=rbind(baron1, baron2, baron3, baron4)
rm(baron1,baron2,baron3,baron4)
cell.meta=baron[,1:3]
colnames(cell.meta)[1]='cell'
cell.meta$individual=substr(cell.meta$cell, start=1, stop=6)
count=as.matrix(baron[,-(1:3)])
rm(baron)
rownames(count)=cell.meta$cell
rownames(cell.meta)=cell.meta$cell
colnames(cell.meta)[3]='celltype'

