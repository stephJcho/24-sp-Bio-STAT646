## Disclaimer: This file has been modified from the original R file (provided by Dr. Yuchao Jiang);
## such as notes taken from the lecture, annotations, extra links and omitting of redundant lines of code.

# Help
?plot.lm
methods(plot) # plot.lm is a hidden/un-exported function -- but it does exist.
# It's common to hide custom methods for generic functions.

# Packages
# install.packages('ggplot2')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges", "org.Hs.eg.db", "biomaRt",
                       "clusterProfiler", "DOSE", "org.Hs.eg.db", "GO.db",
                       "GOSemSim", "enrichplot"), force = TRUE)
library(ggplot2)

################################################################################
################################################################################
####  Bioconductor: Analysis and comprehension of high-throughput genomic data
################################################################################
################################################################################

# Biostrings: DNA, RNA, Amino Acids (AA) representation and manipulation
library(Biostrings)
dna <- DNAStringSet(c("AACTCC", "CTGCA")) #two DNA sequences; first has 6 different nucleotides, other: 5
dna
complement(dna) #gives complement of each nucleotide (in terms of base pair; A-T, G-C)
reverse(dna) #just rewrites in reverse order
reverseComplement(dna) #put in reverse and get complement
reverse(complement(dna)) #same result
# A short example DNA sequence and a search for start / stop codons.
seq <- DNAString("aacataatgcagtagaacccatgagccc")
matchPattern("ATG", seq) #finds matches of the argument in the latter argument
matchPattern("TAA", seq)
#why those two base pairs?
#ATG (AUG in RNA) is the 'start codon' (one of the, in fact, first, amino acids),
#which signals the beginning of protein translation and initiates the formation of the protein chain.
#TAA is the 'stop codon', the last amino acid

# GRanges: The GRanges class represents a collection of genomic ranges that each
# have a single start and end location on the genome. It can be used to store 
# the location of genomic features such as contiguous binding sites, transcripts,
# and exons. These objects can be created by using the GRanges constructor function.
library(GenomicRanges) #can add GC content as metadata
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), #chromosome and their numbers
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)), #start & end position
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)), #which strand (+ or -, complement to each other)
  score = 1:10,
  GC = seq(1, 0, length=10))
gr
seqnames(gr) #which chromosome the seq? comes from
granges(gr) #object, with rows of sequences
ranges(gr)
sort(gr)
mcols(gr) #metadata
mcols(gr)$score
names(gr)
length(gr)

# Problem: gene symbol or annotations for gene types can vary
# But we will learn 3 most used annotations for genes in this course

# 1. Symbols (ex: PTEN [p-ten])
# 2. Ensemble ID (ex: BRCA1)
# 3. Entrez ID (just numbers)

# BioMart: repository of genomic data. You can use BioMart to retrieve genomic
# sequences and annotation of genomic loci. You can also use BioMart to convert 
# between different annotations of the same loci (so that you can align data 
# from different sources). Finally, you can use BioMart to perform enrichment 
# analysis (described below).
library(biomaRt)
help(package = "biomaRt")

# Construct a BioMart dataset consisting of the Ensembl gene documentation for
# humans. Ensembl is a popular genomic data repository.
listDatasets(useMart("ensembl"))
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# use the code below if the previous line doesn't work
mart <- useEnsembl(biomart = "ensembl", mirror = "useast", dataset = "hsapiens_gene_ensembl")

# Extract gene symbols and Ensembl gene id's for each human gene. Actually, these are 
# symbols and id's for all the known coding regions in the human genome. There are only 
# about 20k protein-coding genes in the human genome but many more segments that code 
# for things other than proteins. HGNC is a collection of unique symbols and names for 
# genomic loci in humans. For example, HGNC has a unique symbol for each known protein-
# coding gene. In addition to protein-coding genes, HGNC has symbols / names for "RNA 
# genes" (DNA loci that code for RNA that do not code for proteins but perform other 
# functions) and pseudogenes (genomic DNA sequences that are similar to normal genes but 
# non-functional...thought to be defunct relatives of functional genes). 
listAttributes(mart)[1:10,]
my_results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), mart = mart)
# key for this class: matching symbols of the genes to ensemble id of the genes

head(my_results) # Enter gene symbols into the HGNC website to get their details
# MT: mitochondria; first set of genes, has own set of DNA although it is outside of nucleus
tail(my_results) # Genes have HGNC symbol but does have an Ensembl ID.
# According to Ensembl, these DNA loci are known pseudogenes.
# 70000 smth from this human query, while human have around 20000 genes
# Reason: there are non-traditional, additional sets of genes, ex) rows without symbol or smth

# List the peptide sequences for the exons of MYC gene. An exon is a piece of a gene 
# that is required for coding the gene's product. The exons are separated by introns. 
# The introns are removed prior to gene expression, and the exons are gathered together 
# into a single sequence fragment that is then translated into protein.
#
# "MYC" is an HGNC symbol. The MYC gene plays an important role in the development of 
# cancer: http://en.wikipedia.org/wiki/Myc
seq <- getSequence(id="MYC", type="hgnc_symbol", seqType="peptide", mart = mart)
# peptide/amino acid/protein sequence; can change argument to 'cdna' to get cDNA, etc.
show(seq) #9 different rows, or peptide sequences for this gene
# reason: each rows have different length (alternative slicing)
# related to BRCA1 example from Lec 1
?getSequence

# Entrez ID conversions. Converting one id type to another. The org.Hs.eg.db
# package contains the genome-wide annotation for humans. Below, we will 
# convert Entrez gene id's to gene symbols, and vice-versa.
library(org.Hs.eg.db)

# Here is a list of Entrez gene id's. The last one is a non-existent Entrez ID (37690). 
# This code will convert Entrez ID's to gene symbols.
EIDs <- c("1", "10", "100", "1000", "5982", "37690") # designed so that the result won't give any hits

ls("package:org.Hs.eg.db")
EID_symbols <- mget(EIDs, org.Hs.egSYMBOL, ifnotfound=NA)
EID_symbols <- unlist(EID_symbols)
EID_symbols
EID_symbols <- EID_symbols[!is.na(EID_symbols)]
EID_symbols

# This code will convert gene symbols back to Entrez ID's.
EIDs_rev <- unlist(mget(EID_symbols, org.Hs.egSYMBOL2EG, ifnotfound = NA))
EIDs_rev

# Exercise: finding Ensemble ID and EID(Entrez ID) of BRCA1
unlist(mget('BRCA1', org.Hs.egSYMBOL2EG, ifnotfound = NA)) #EID: 672
my_results[which(my_results$hgnc_symbol=='BRCA1'),] #Ensemble ID: ~~12048 (from last col)
head(my_results) #my_results: data frame with symbol and Ensemble ID
# Some packages take input as EID by default, others as Ensemble ID

# 3rd method of annotation I guess?
# GO annotation. Gene Ontology ("GO") data characterizes genes in terms of 
# their cellular component (where it is located in the cell), molecular function
# (what function in the cell the gene is involved in), and biological process 
# (a collection of molecular functions that define a higher-order biological function).
# 
# Get GO annotations for human genes. Part of this annotation information is GO data,  
# accessible with the org.Hs.egGO* functions described below.

EIDs <- c("1", "10", "100") # Some example genes (Entrez id's).

# We get a list of lists. One list per EID that we supplied. Inside the list for any one 
# EID, we have a variety of entries, one or more for each of cellular component 
# (indicated by a "CC" in the "Ontology" slot), biological process (indicated by a "BP" 
# in the Ontology slot), and molecular function (indicated by a "MF" in the Ontology 
# slot). There is also an Evidence slot, with codes defined in the help file for the 
# org.Hs.egGO function. Finally, there is a GOID slot that provides a GO code that you 
# can look up at the GO website: http://geneontology.org/.

## GO: a group of genes in a specific biological function (or pathway?); a category of genes doing same shit

# The gene with Entrez id 10 is known to be involved in the biological process (BP):
# xenobiotic metabolic process, which comes from "traceable author statement" 
# evidence (GO:0006805). This gene is also known to be involved in the molecular 
# function (MF) arylamine N-acetyltransferase activity (evidence = "inferred from 
# biological aspect of ancestor" (IBA); GO:0004060).
GO <- mget(EIDs, org.Hs.egGO) #a list, each element is the Gene Ontology that nth gene is involved in
length(GO)
GO[[2]] # GO that the second gene (with EID 10) is involved in
# 2nd gene is part of '0006805' GO, which consists of multiple different genes, including EID 10
?org.Hs.egGO # to find what evidence elements are and their codes

# We can also retrieve all the genes involved in a particular GO code. Here are the 
# genes involved in this molecular function: arylamine N-acetyltransferase activity. 
# This corresponds to the GO code GO:0004060. Notice the gene with Entrez id 10 is 
# included. This code was returned above for this gene.
GOterm <- "GO:0004060"
GOgenes <- mget(GOterm, org.Hs.egGO2ALLEGS)

# We can extract the GO material for a given GO id using the GO.db package. This way, we 
# don't have to go to the gene ontology website.
BiocManager::install("GO.db")
library(GO.db)
GO_data = list(GOID(GOterm), Term(GOterm), Synonym(GOterm),
               Secondary(GOterm), Definition(GOterm), Ontology(GOterm))
GO_data

# GO enrichment. This involves using hypothesis tests to assess whether GO id's are 
# over-represented (relative to what might be expected by chance) in a list of 
# selected genes. The most widely used test is the Fisher's Exact Test, which involves 
# a hypergeometric calculation. For example, suppose we begin with 25,000 genes and 
# select 100 as most interesting. Consider GO id GO:0000002, with 23 genes associated 
# with it. It is as if there are 25,000 balls, of which 23 are in GO:0000002. We draw 
# 100 and ask for the probability that we observe K or more genes from GO:0000002.

# We do this for each GO id, resulting in many p-values (18,826 for humans). Those GO 
# id's with p-values below a selected threshold are said to be "enriched" relative to 
# our selected genes.
library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(enrichplot)

# geneList data set which was employed in DOSE package is a good choice.
# It contains IDs of the genes of interest (e.g., differentially expressed genes).
## I can just provide the package the list of the genes (I am interested), and it does calculation (hypergeom), gives p-val and plot
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene # genes I am interested in; n number of genes, in a hypergeometric argument
length(geneList)
length(gene)

# Giving 207 genes from prior analysis

# The enrichGO() function performs the GO Enrichment Analysis on a given vector 
# of genes. Enrichment analysis is an approach for distinguishing a group of genes 
# that are designated to a class of predefined bins according to their functional 
# specifications. The enriched outcome may contain very general terms. To use 
# this function, you have to fill out an argument called “OrgDb”. For that, the
# org.Hs.eg.db package including human genome-wide annotation is required. 
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC", # Cellular Component (CC), Biological Process (BP), Molecular Function (MF)
                pAdjustMethod = "BH", # Benjamini-Hochberg multiple testing correction
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

# The clusterProfiler package is a set of methods specified to analyze and visualize 
# functional profiles like GO of genes and gene clusters. By using this, you will 
# be able to cluster different genes according to their similarities. 
d <- godata('org.Hs.eg.db', ont="CC")
ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
emapplot(ego2)
emapplot_cluster(ego2)

# each circle of the plot is one specific GO term that the 207 genes are (significantly) enriched in
# the latter gives also the clusters 
