#DE_analysis Corallium rubrum RNAseq

#set working directory
setwd("C:/Users/Sandra/OneDrive - The University of Hong Kong/Doctorado/ICM/Tesis/Chap2/RNA-Seq/DESeq2") #mypc
setwd("C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/DESeq2") #Quim's PC

#what you need:

# Packages
library(digest)
library(XML)
library(RSQLite)
library(DESeq2)
library("pcaExplorer")
library("ComplexHeatmap")
library("gplots")
library("RColorBrewer")
library("gplots") #heatmap.2 function
library("pheatmap")
library("dplyr")
library("ggplot2")

# Load data (count_data)
cts <- as.matrix(read.csv("crubrum.gene.counts.matrix.csv", row.names="SeqName"))
View(cts)

# Load sample info,conditions,treatments,etc (coldata)
coldata <- read.csv("crubrum.coldata.csv", row.names = 1)

# Change decimal to integers
cts=round(cts)

########################## control vs treatment ##################################

# check if all colnames of count matrix are in the rownames of sample info
(all(rownames(coldata) %in% colnames(cts)) || all(colnames(cts) %in% rownames(coldata)))

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts) == rownames(coldata)))

# the variables in the design formula (treatment, day, population, or 
#Individual) are numeric variables with integer values, and the DESeq2 
#package is interpreting them as continuous variables. 
#using a numeric variable as a predictor in the model implies that the fold 
#change between groups would increase as the variable increases, 
#which may not be biologically meaningful.

#To fix this, you can convert these variables to factors using the factor() 
#function:

coldata$day <- factor(coldata$day)
coldata$population <- factor(coldata$population)
coldata$Individual <- factor(coldata$Individual)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~treatment)#differences btw Cont vs Treat

# Remove genes which have zero reads in all samples
keep <- rowSums(counts(dds)) >= 10 #suggested in manual
dds <- dds[keep,]

# in case you want to obtain normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
norm_cts_C_vs_treat <- counts(dds, normalized=TRUE)
norm_cts_log <- log2(norm_cts_C_vs_treat) #log transform normalised counts
write.csv(norm_cts_log, file="C:/Users/Joaquim Garrabou/OneDrive - The University of Hong Kong/Doctorado/ICM/Tesis/Chap2/RNA-Seq/crubrum_rnaseq/1.dds~treat/log2_normalised_counts_C-T_crubrum_10cts.csv", row.names=T)

#Set factor levels for comparison
#dds$treatment <- factor(dds$treatment, levels = c("Control", "Treatment"))
dds$treatment <- relevel(dds$treatment, ref = "Control")

#Differential expression
dds <- DESeq(dds)#,parallel
dds

save(dds, file = "C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2/dds.RData")
#load(file="C:/Users/Joaquim Garrabou/OneDrive - The University of Hong Kong/Doctorado/ICM/Tesis/Chap2/RNA-Seq/crubrum_rnaseq/1.dds~treat/dds.RData")

#obtain results from dds object
res <- results(dds, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"))
res #notice the upregulated and downregulated (+ or -, respectively) -
# the log2 fold change and Wald test p value will be for the last variable (treatment)
#in the design formula, and if this is a factor (it is), the comparison will be the 
#last level of this variable over the reference level (so, it's values for treatment)
#filter baseMean genes > 10
res <- res[res$baseMean > 10, ]
summary(res) 

save(res, file = "C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2/res.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

plotDispEsts(dds)

#order res based on ajusted smallest p-values
resOrdered <- res[order(res$padj),]
summary(res)

######
plotMA(res, ylim=c(-8,8)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotCounts(dds, gene = which.min(res$padj), intgroup = "treatment")

mcols(res)$description #Information about which variables and tests were used 
#can be found by calling the function mcols on the results object.

#[1] "mean of normalized counts for all samples"             
#[2] "log2 fold change (MLE): treatment Treatment vs Control"
#[3] "standard error: treatment Treatment vs Control"        
#[4] "Wald statistic: treatment Treatment vs Control"        
#[5] "Wald test p-value: treatment Treatment vs Control"     
#[6] "BH adjusted p-values" 

#extract data from res object
re <- data.frame(resOrdered)
re <- na.omit(re)
write.csv(re,file="C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2/deseq2_C-T_crubrum.csv") #This is the result of the analysis. Name it well.

#extract data which pass an adjusted p value threshold (0.05)
sig=subset(re,re$padj<0.05)
sig2=subset(sig,abs(sig$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig2, file = "C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2/sig2_deseq2_C-T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig2,file="C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2/diff_genes_C-T_crubrum.csv") #save as cvs file


############# Data transformation and visualization ############################
# for visualization or clustering – it might be useful to work with 
#transformed versions of the count data.
#The most useful is the log as you have zeros and non-zeros in others - big diff

#VST and rlog - remove the dependence of the variance on the mean, 
#particularly the high variance of the logarithm of #count data when the mean 
#is low. Both VST and rlog use the experiment-wide trend of variance over mean, 
#in order to transform the data to remove the experiment-wide trend

#VST - variance stabilizing transformations (VST)
#rlog - incorporates a prior on the sample differences and 
#requires fitting a shrinkage term for each sample and each gene which takes time

#Guidelines:

#    VST is suitable for differential expression analysis: 
#VST can be a good choice if you are interested in identifying 
#differentially expressed genes, as it stabilizes the variance of the counts, 
#which can improve the accuracy of differential expression analysis.

#rlog is suitable for exploratory data analysis: rlog can be a good choice if 
#you are interested in exploratory data analysis, such as clustering 
#or visualization, as it can reduce the effects of outliers and make the 
#data more normally distributed.

#Consider the sample size: VST tends to work better for larger sample sizes, 
#while rlog can be more effective for smaller sample sizes.

#Consider the biological variability: VST assumes that the biological 
#variability is independent of the mean expression level, while rlog does
#not make this assumption. If you have reason to believe that the biological variability is dependent on the mean expression level, rlog may be a better choice.

#Consider the downstream analysis: Some downstream analyses, such as gene set 
#enrichment analysis, may have different requirements for input data. 
#It is important to check the requirements of your downstream analysis 
#to ensure that your choice of transformation is appropriate.

#So, we choose rlog

#rld <- vst(dds,blind = FALSE) #faster transformation
rld <- rlog(dds, blind = FALSE)
save(rld, file = "C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2/rld.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/rld.RData")

######################## Heatmap C-T ###########################

#heatmap of the count matrix based on rld values:

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("treatment","population")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances without clustering
sampledist <- dist(t(assay(rld))) # calculate the euclidean distance between the samples
sampledist #gives us an overview over similarities and dissimilarities between samples
sampledistmatrix <- as.matrix(sampledist)
colnames(sampledistmatrix) <- paste(rld$treatment,rld$population,sep="-")
colnames(sampledistmatrix) <- NULL

pheatmap(sampledistmatrix, clustering_distance_rows=sampledist,
           clustering_distance_cols=sampledist)

#test for similarities in gene expression bwn samples. How similar the repliates are for each other
#and whether the samples from a given groups cluster each other or not
#0-1 correlation values. We expect >0.80, otherwise they may be outliers

## heatmap with gene clustering
sig2
topVarGenes <- head(order( sig2$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rld)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

##### PCA C-T #####

head(assay(rld))
plotPCA(rld, intgroup="treatment")
plotPCA(rld, intgroup=c("population"))
plotPCA(rld, intgroup=c("treatment","population"))

pcaExplorer(dds = dds, dst = rld)

plotPCA(rld, intgroup = "treatment",
        ntop = 500, returnData = FALSE)

