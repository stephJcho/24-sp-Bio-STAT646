####
#### Genome-Wide Association Studies
####
devtools::install_github("isglobal-brge/SNPassoc")

library(SNPassoc)

## We will use an example data set from the SNPassoc package. This is a very small subset 
## of a case-control GWAS. There are 35 SNPs, each genotyped on 157 individuals. For each 
## individual, we have gender, arterial blood pressure and "protein levels" (not sure 
## what they mean by this exactly). Two objects are loaded. 'SNPs' contains the genotypes 
## along with the covariate information for each individual. 'SNPs.info.pos' contains 
## names for each SNP, as well as its genomic location.
data(SNPs)
?SNPs
data(SNPs.info.pos)

head(SNPs)
head(SNPs.info.pos)

## Create a 'SNP' object.
Y_SNP <- setupSNP(SNPs, colSNPs = 6:40, info = SNPs.info.pos, sep = "")

## Carry out an association test for the first SNP. This just does logistic regression 
## (response = case / control), assessing the odds ratio associated with 'snp10001' and 
## adjusting for 'sex' and 'blood.pre'. Note that the way we define the 'snp10001' 
## covariate depends on the type of genotypes we assume. The 'snp10001' SNP is either T 
## or C. If the C is dominant, then CT and CC will both give the C trait, while only TT 
## will give the T trait. If C is recessive, then only CC will give the C trait, while 
## CT and TT will both give the T trait. There are a couple of other ways to define this. 
## Finally, we can define genotype as "log additive." With 'snp10001', this would mean 
## coding TT as 0, CT as 1, and CC as 2 (i.e., treating genotype as an ordinal variable). 
## The 'association' function below reports results under all possible genotype 
## definitions by default, but you can specify one or more via the 'model' argument. In 
## the following, the estimated odds ratio associated with a one-unit increase in 
## genotype on the log-additive scale is 0.87 (95% CI: [0.51, 1.49], p = 0.61). This is 
## not statistically significant, indicating that this SNP is not associated with case / 
## control status. Note that the effect estimates and statistical significance varies 
## depending on how we define genotype.
assoc_out <- association(casco ~ sex + blood.pre + snp10001, data = Y_SNP)
assoc_out

## The function 'WGassociation' is like 'assocation', but it is applied to all SNPs at 
## once (so we get p-values for each SNP). Here, we use 'WGassociation', using the log-
## additive genotype model. The 'WGstats' function reports tabulations, OR estimates, and 
## p-values for each of the SNPs. Note that some of the SNPs are "monomorphic"; that is, 
## their genotype does not differ for any of our individuals. We cannot estimate an OR 
## for such SNPs. The 'plot' function shows -log10(p-values) for each SNP.
wg_assoc_out <- WGassociation(casco ~ sex + blood.pre, data = Y_SNP, model = "log-additive")
WGstats(wg_assoc_out)
plot(wg_assoc_out)
pvalues(wg_assoc_out)

## In addition to testing association with a categorical response, we can use a 
## quantitative response. We have a variable called 'protein' in our example data set.
## This does linear regression.
prot_assoc_out <- association(log(protein) ~ sex + blood.pre + snp10001, data = Y_SNP)
prot_assoc_out
pvalues(wg_assoc_out)

## Here is a more interesting example. These are data from the HapMap project, consisting 
## of 9305 SNPs (across 22 chromosomes) for 120 individuals. Sixty of the individuals 
## come from the "central European" population, and sixty of the individuals come from 
## the Yoruba (in Africa) population. Note that this takes a while to run. The plot again 
## shows -log10(p-values), now visualized spatially along the chromosomes.
data(HapMap)
data(HapMap.SNPs.pos)

hapmap_SNP <- setupSNP(HapMap, colSNPs = 3:9307, sort = TRUE, 
  info = HapMap.SNPs.pos, sep = "")
plot(hapmap_SNP$rs6659552)
plot(hapmap_SNP$rs6659552, type=pie)
plotMissing(hapmap_SNP, print.labels.SNPs = FALSE) # Plot a grid showing which genotypes are missing

hapmap_wg_assoc_out <- WGassociation(group, data = hapmap_SNP, model = "log-additive")
plot(hapmap_wg_assoc_out)

hist(pvalues(hapmap_wg_assoc_out)[,2],20, main='Histogram of p-values')

padj.bonferroni=p.adjust(pvalues(hapmap_wg_assoc_out)[,2], method='bonferroni')
padj.BH=p.adjust(pvalues(hapmap_wg_assoc_out)[,2], method='BH')

sum(padj.bonferroni<=0.05, na.rm=TRUE)
sum(padj.BH<=0.05, na.rm=TRUE)



## Now try it yourself
## https://cran.r-project.org/web/packages/SNPassoc/vignettes/SNPassoc.html
data(asthma)
asthma.s <- setupSNP(data=asthma, colSNPs=7:ncol(asthma), sep="")

asthma.s_wg_assoc_out <- WGassociation(casecontrol~gender+age+bmi+smoke,
                                       data = asthma.s, model = "log-additive")


hist(pvalues(asthma.s_wg_assoc_out)[,2],20, main='Histogram of p-values')
padj.bonferroni=p.adjust(pvalues(asthma.s_wg_assoc_out)[,2], method='bonferroni')
padj.BH=p.adjust(pvalues(asthma.s_wg_assoc_out)[,2], method='BH')

sum(padj.bonferroni<=0.1, na.rm=TRUE)
sum(padj.BH<=0.1, na.rm=TRUE)
sum(pvalues(asthma.s_wg_assoc_out)[,2]<=0.05)
