###########################################################################
#                                                                         #
########### All variables to play ####### "Top-down" approach #############
#                                                                         #
###########################################################################

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

#remove outlier LOP_C0_1

cts1 <- cts[,-22]
coldata1 <- coldata [-22, ]

# doublecheck if all colnames of count matrix are in the rownames of sample info
(all(rownames(coldata1) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata1)))

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata1)))

dds1 <- DESeqDataSetFromMatrix(countData = cts1, 
                               colData = coldata1, 
                               design = ~day + Individual + population +treatment)#differences btw Cont vs Treat

# Remove genes which have zero reads in all samples
keep <- rowSums(counts(dds1)) >= 10 #suggested in manual
dds1 <- dds1[keep,]

# in case you want to obtain normalized counts
dds1 <- estimateSizeFactors(dds1)
sizeFactors(dds1)
norm_cts_C_vs_treat <- counts(dds, normalized=TRUE)
norm_cts_log <- log2(norm_cts_C_vs_treat) #log transform normalised counts
write.csv(norm_cts_log, file="C:/Users/Joaquim Garrabou/OneDrive - The University of Hong Kong/Doctorado/ICM/Tesis/Chap2/RNA-Seq/crubrum_rnaseq/1.dds~treat/log2_normalised_counts_C-T_crubrum_10cts.csv", row.names=T)

#Set factor levels for comparison
#dds$treatment <- factor(dds$treatment, levels = c("Control", "Treatment"))
dds1$treatment <- relevel(dds1$treatment, ref = "Control")

#Differential expression
dds1 <- DESeq(dds1)#,parallel
dds1

save(dds1, file = "dds1_sampleout.RData")
#load(file="dds1_sampleout.RData")

#obtain results from dds object
res1 <- results(dds1, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"))#obtain the results of the differential expression analysis of two groups of samples ("Treatment" and "Control") 
res1 #notice the upregulated and downregulated (+ or -, respectively) -
# the log2 fold change and Wald test p value will be for the last variable (treatment)
#in the design formula, and if this is a factor (it is), the comparison will be the 
#last level of this variable over the reference level (so, it's values for treatment)
res1 <- res1[res1$baseMean > 10, ]
summary(res1)

save(res1, file="res1_sampleout.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

plotDispEsts(dds1)

#order res based on ajusted smallest p-values
res1Ordered <- res1[order(res1$padj),]
summary(res1)

######
plotMA(res1, ylim=c(-12,12)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotCounts(dds1, gene = which.min(res1$padj), intgroup = "treatment")

mcols(res1)$description #Information about which variables and tests were used 
#can be found by calling the function mcols on the results object.

#[1] "mean of normalized counts for all samples"             
#[2] "log2 fold change (MLE): treatment Treatment vs Control"
#[3] "standard error: treatment Treatment vs Control"        
#[4] "Wald statistic: treatment Treatment vs Control"        
#[5] "Wald test p-value: treatment Treatment vs Control"     
#[6] "BH adjusted p-values" 

#extract data from res object
re1 <- data.frame(res1Ordered)
re1 <- na.omit(re1)
write.csv(re1,file="deseq2_allfactors_C-T_crubrum_sampleout.csv") #This is the result of the analysis. Name it well.

#extract data which pass an adjusted p value threshold (0.05)
sig1=subset(re1,re1$padj<0.05)
sig21=subset(sig1,abs(sig1$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig21, file = "sig2_deseq2_allfactors_C-T_crubrum_sampleout.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig21,file="diff_genes_allfactors_C-T_crubrum_sampleout.csv") #save as cvs file


### Data transformation and visualization ########
#So, we choose rlog

#rld <- vst(dds1,blind = FALSE) #faster transformation
rld1 <- rlog(dds1, blind = FALSE)
save(rld1, file = "rld1_sampleout.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/rld.RData")

### Heatmap C-T #######

#heatmap of the count matrix based on rld values:
select <- order(rowMeans(counts(dds1,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds1)[,c("treatment","day", "Individual", "population")])
pheatmap(assay(rld1)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances without clustering
sampledist <- dist(t(assay(rld1))) # calculate the euclidean distance between the samples
sampledist #gives us an overview over similarities and dissimilarities between samples
sampledistmatrix <- as.matrix(sampledist)
colnames(sampledistmatrix) <- paste(rld1$treatment,rld1$population,sep="-")
colnames(sampledistmatrix) <- NULL

pheatmap(sampledistmatrix, clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist)

#test for similarities in gene expression bwn samples. How similar the repliates are for each other
#and whether the samples from a given groups cluster each other or not
#0-1 correlation values. We expect >0.80, otherwise they may be outliers

## heatmap with gene clustering
sig21
topVarGenes <- head(order( sig21$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rld1)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

##### PCA C-T #####

head(assay(rld1))
plotPCA(rld1, intgroup="treatment")
plotPCA(rld1, intgroup=c("population"))
plotPCA(rld1, intgroup=c("treatment","population", "day"))
plotPCA(rld1, intgroup= c("treatment","Individual"))

pcaExplorer(dds = dds1, dst = rld1)

plotPCA(rld1, intgroup = "treatment",
        ntop = 500, returnData = FALSE)

############################################### 
######## Removing the factor "day" ############
###############################################

#remove new outlier CAS_T0_1 from cts1 and coldata1

cts2 <- cts1[,-10]
coldata2 <- coldata1[-10, ]

# doublecheck if all colnames of count matrix are in the rownames of sample info
(all(rownames(coldata1) %in% colnames(cts1)) || all(colnames(cts1) %in% rownames(coldata1)))

# check if columns of count matrix are in same order as rows of sample info
(all(colnames(cts1) == rownames(coldata1)))

#work with cts1 that does not have an outlier 
dds2 <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata2, design = ~ Individual + population +treatment)#differences btw Cont vs Treat
#we want to account for the effect of individual, population, 
#and treatment in the analysis of the RNA-seq data. This design is 
#commonly used when we want to test the differential expression of genes 
#between two experimental conditions, while accounting for differences between 
#individuals and/or populations.

# Remove genes which have zero reads in all samples
keep <- rowSums(counts(dds2)) >= 10 #suggested in manual
dds2 <- dds2[keep,]

# in case you want to obtain normalized counts
dds2 <- estimateSizeFactors(dds2)
sizeFactors(dds2)
norm_cts_C_vs_treat <- counts(dds, normalized=TRUE)
norm_cts_log <- log2(norm_cts_C_vs_treat) #log transform normalised counts
write.csv(norm_cts_log, file="C:/Users/Joaquim Garrabou/OneDrive - The University of Hong Kong/Doctorado/ICM/Tesis/Chap2/RNA-Seq/crubrum_rnaseq/1.dds~treat/log2_normalised_counts_C-T_crubrum_10cts.csv", row.names=T)

#Set factor levels for comparison
#dds$treatment <- factor(dds$treatment, levels = c("Control", "Treatment"))
dds2$treatment <- relevel(dds2$treatment, ref = "Control")

#Differential expression
dds2 <- DESeq(dds2)#,parallel
dds2

save(dds2, file = "dds2_sampleout.RData")
#load(file="dds2.RData")

#obtain results from dds object
res2 <- results(dds2, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"))#obtain the results of the differential expression analysis of two groups of samples ("Treatment" and "Control") 
res2 #notice the upregulated and downregulated (+ or -, respectively) -
# the log2 fold change and Wald test p value will be for the last variable (treatment)
#in the design formula, and if this is a factor (it is), the comparison will be the 
#last level of this variable over the reference level (so, it's values for treatment)
res2 <- res2[res2$baseMean > 10, ]#remove genes with a baseMean<10
summary(res2)

save(res2, file="res2_sampleout.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

plotDispEsts(dds2)

#order res based on ajusted smallest p-values
res2Ordered <- res2[order(res2$padj),]
summary(res2)

plotMA(res2, ylim=c(-12,12)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotCounts(dds2, gene = which.min(res2$padj), intgroup = "treatment")

mcols(res2)$description #Information about which variables and tests were used 
#can be found by calling the function mcols on the results object.

#[1] "mean of normalized counts for all samples"             
#[2] "log2 fold change (MLE): treatment Treatment vs Control"
#[3] "standard error: treatment Treatment vs Control"        
#[4] "Wald statistic: treatment Treatment vs Control"        
#[5] "Wald test p-value: treatment Treatment vs Control"     
#[6] "BH adjusted p-values" 

#extract data from res object
re2 <- data.frame(res2Ordered)
re2 <- na.omit(re2)
write.csv(re2,file="deseq2_allfactors_noday_C-T_crubrum_sampleout.csv") #This is the result of the analysis. Name it well.

#extract data which pass an adjusted p value threshold (0.05)
sig2=subset(re2,re2$padj<0.05)
sig22=subset(sig2,abs(sig2$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig22, file = "sig2_deseq2_allfactors_noday_C-T_crubrum_sampleout.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig22,file="diff_genes_allfactors_noday_C-T_crubrum_sampleout.csv") #save as cvs file


### Data transformation and visualization #####
#So, we choose rlog

#rld <- vst(dds2,blind = FALSE) #faster transformation
rld2 <- rlog(dds2, blind = FALSE)

save(rld2, file = "rld2_samplesout.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/rld.RData")

##### Heatmap C-T ######

#heatmap of the count matrix based on rld values:
select <- order(rowMeans(counts(dds2,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds2)[,c("treatment", "Individual", "population")])
pheatmap(assay(rld2)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances without clustering
sampledist <- dist(t(assay(rld2))) # calculate the euclidean distance between the samples
sampledist #gives us an overview over similarities and dissimilarities between samples
sampledistmatrix <- as.matrix(sampledist)
colnames(sampledistmatrix) <- paste(rld2$treatment,rld2$population,sep="-")
colnames(sampledistmatrix) <- NULL

pheatmap(sampledistmatrix, clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist)

#test for similarities in gene expression bwn samples. How similar the repliates are for each other
#and whether the samples from a given groups cluster each other or not
#0-1 correlation values. We expect >0.80, otherwise they may be outliers

### heatmap with gene clustering ###
sig22
topVarGenes <- head(order( sig22$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rld2)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

##### PCA C-T #####

head(assay(rld2))
plotPCA(rld2, intgroup="treatment")
plotPCA(rld2, intgroup=c("population"))
plotPCA(rld2, intgroup=c("treatment","population", "day"))
plotPCA(rld2, intgroup="day")

pcaExplorer(dds = dds2, dst = rld2)

plotPCA(rld2, intgroup = "treatment",
        ntop = 500, returnData = FALSE)

###### All factors minus Individual and day: ~population + treatment ######

#Testing the effect of the combination of population and treatment means that you 
#are testing whether the effect of treatment on gene expression is different across 
#different populations. This means that you are interested in knowing whether there 
#is an interaction effect between population and treatment.--> Grouping (pop+treat)

#On the other hand, testing the independent effect of population and treatment means
#that you are testing whether treatment or population alone have an effect on gene 
#expression, regardless of the other variable. This means that you are not 
#interested in the interaction effect between population and treatment, but rather 
#in their individual effects. --> ~population + treatment

#In summary, testing the effect of the combination of population and treatment is
#more specific and allows you to identify genes that are differentially expressed 
#due to the interaction between these two factors, while testing the independent 
#effect of population and treatment allows you to identify genes that are 
#differentially expressed due to the main effect of each factor alone.

#Is the treatment effect different across populations?

ddsMF <- dds2

#Using this design is similar to adding an interaction term, 
#in that it models multiple condition effects which can be easily 
#extracted with results. We have two factors population 
#(with values CAS and LOP) and treatment (with values Control and Treatment), 
#and we want to extract the condition effect specifically for 
#each population. We could use the following approach to obtain, e.g.
#the condition effect for population CAS and LOP:

ddsMF$group <- factor(paste0(ddsMF$population, ddsMF$treatment))
design(ddsMF) <- ~group #we want to model the effect of "treatment", "population", and the
#interaction between "treatment" and "population". This allows us to
#test if the treatment effect differs between the "CAS" and "LOP" 
#populations.

ddsMF_pop <- DESeq(ddsMF)
resultsNames(ddsMF_pop)
#[1] "Intercept"                        "group_CASTreatment_vs_CASControl"
#[3] "group_LOPControl_vs_CASControl"   "group_LOPTreatment_vs_CASControl"

# obtain results using contrast for "CAS" and "LOP" accounting for "treatment" effect
res_CAS_pop <- results(ddsMF_pop, alpha = 0.05, contrast = c("group", "CASTreatment", "CASControl"), test = "Wald")
res_CAS_pop <- res_CAS_pop[res_CAS_pop$baseMean > 10, ] #remove samples with a baseMean < 10
res_CAS_pop

res_LOP_pop <- results(ddsMF_pop, alpha = 0.05, contrast = c("group", "LOPTreatment", "LOPControl"), test = "Wald")
res_LOP_pop <- res_LOP_pop[res_LOP_pop$baseMean > 10, ] #remove samples with a baseMean < 10
res_LOP_pop

save(ddsMF_pop, file = "ddsMF.RData")
#load(file="C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/dds.RData")

#obtain results using contrast for control vs treatment btw pops
res_Cont_CAS_LOP <- results(ddsMF_pop,alpha = 0.05, contrast = c("group", "CASControl", "LOPControl"), test = "Wald") #Diff gene exp btw Controls
res_Cont_CAS_LOP <- res_Cont_CAS_LOP[res_Cont_CAS_LOP$baseMean > 10, ] #remove samples with a baseMean < 10
res_Cont_CAS_LOP

res_treat_CAS_LOP <- results(ddsMF_pop,alpha = 0.05, contrast = c("group", "CASTreatment", "LOPTreatment"), test = "Wald") #Diff gene exp btw treats
res_treat_CAS_LOP <- res_treat_CAS_LOP[res_treat_CAS_LOP$baseMean > 10, ] #remove samples with a baseMean < 10
res_treat_CAS_LOP

#notice the upregulated and downregulated (+ or -, respectively)
summary(res_CAS_pop)
summary(res_LOP_pop)
summary(res_Cont_CAS_LOP)
summary(res_treat_CAS_LOP)

#CAS
save(res_CAS_pop, file = "res_CAS.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

#LOP
save(res_LOP_pop, file = "res_LOP.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

#control
save(res_Cont_CAS_LOP, file = "res_cont_CAS-LOP.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

#treatment
save(res_treat_CAS_LOP, file = "res_treat_CAS-LOP.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

plotDispEsts(ddsMF_pop)

#order res based on ajusted smallest p-values
resOrdered_CAS <- res_CAS_pop[order(res_CAS_pop$padj),]
summary(res_CAS_pop)

resOrdered_LOP <- res_LOP_pop[order(res_LOP_pop$padj),]
summary(res_LOP_pop)

resOrdered_control <- res_Cont_CAS_LOP[order(res_Cont_CAS_LOP$padj),]
summary(res_Cont_CAS_LOP)

resOrdered_treatment <- res_treat_CAS_LOP[order(res_treat_CAS_LOP$padj),]
summary(res_treat_CAS_LOP)

#plotMA
plotMA(res_Cont_CAS_LOP, ylim=c(-8,8)) 
plotMA(res_treat_CAS_LOP, ylim=c(-8,8)) 

plotMA(res_CAS_pop, ylim=c(-8,8)) 
plotMA(res_LOP_pop, ylim=c(-8,8)) 

plotCounts(ddsMF_pop, gene = which.min(res_CAS_pop$padj), intgroup = "group")

#extract data from res object

#CAS
reCAS <- data.frame(resOrdered_CAS)
reCAS <- na.omit(reCAS)
write.csv(reCAS,file="deseq2_CAScontrol-CAStreatment_crubrum.csv") #This is the result of the analysis. Name it well.
#extract data which pass an adjusted p value threshold (0.05)
sigCAS=subset(reCAS,reCAS$padj<0.05)
sig2CAS=subset(sigCAS,abs(sigCAS$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig2CAS, file = "sig2_deseq2_CASC-CAST_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig2CAS,file="diff_genes_CASC-CAST_crubrum.csv") #save as cvs file

#LOP
reLOP <- data.frame(resOrdered_LOP)
reLOP <- na.omit(reLOP)
write.csv(reLOP,file="deseq2_LOPcontrol-LOPtreatment_crubrum.csv") #This is the result of the analysis. Name it well.
#extract data which pass an adjusted p value threshold (0.05)
sigLOP=subset(reLOP,reLOP$padj<0.05)
sig2LOP=subset(sigLOP,abs(sigLOP$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig2LOP, file = "sig2_deseq2_LOPC-LOPT_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig2LOP,file="diff_genes_LOPC-LOPT_crubrum.csv") #save as cvs file

#Controls
reControl <- data.frame(resOrdered_control)
reControl <- na.omit(reControl)
write.csv(reControl,file="deseq2_CAScontrol-LOPcontrol_crubrum.csv") #This is the result of the analysis. Name it well.
#extract data which pass an adjusted p value threshold (0.05)
sigControl=subset(reControl,reControl$padj<0.05)
sig2Control=subset(sigControl,abs(sigControl$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig2Control, file = "sig2_deseq2_CASC-LOPC_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig2Control,file="diff_genes_CASC-LOPC_crubrum.csv") #save as cvs file

#treatments
retre <- data.frame(resOrdered_treatment)
retre <- na.omit(retre)
write.csv(retre,file="deseq2_CAST-LOPT_crubrum.csv") #This is the result of the analysis. Name it well.
#extract data which pass an adjusted p value threshold (0.05)
sigtre=subset(retre,retre$padj<0.05)
sig2tre=subset(sigtre,abs(sigtre$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig2tre, file = "sig2_deseq2_CAST-LOPT_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig2tre,file="diff_genes_CAST-LOPT_crubrum.csv") #save as cvs file

#So, we choose rlog
rldpop <- rlog(ddsMF_pop, blind = FALSE)

save(rldpop, file = "rld.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/rld.RData")

##### Heatmap CAS-LOP #######

#heatmap of the count matrix based on rld values:

select <- order(rowMeans(counts(ddsMF_pop,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsMF_pop)[,c("group", "population", "treatment")])
pheatmap(assay(rldpop)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances without clustering
sampledist <- dist(t(assay(rldpop))) # calculate the euclidean distance between the samples
sampledist #gives us an overview over similarities and dissimilarities between samples
sampledistmatrix <- as.matrix(sampledist)
colnames(sampledistmatrix) <- paste(rldpop$group,rldpop$population,sep="-")
colnames(sampledistmatrix) <- NULL

pheatmap(sampledistmatrix, clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist)

## heatmap with gene clustering
sig2CAS
topVarGenes <- head(order(sig2CAS$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rldpop)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

sig2LOP
topVarGenes <- head(order(sig2LOP$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rldpop)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

##### PCA #####

head(assay(rldpop))
plotPCA(rldpop, intgroup="group")

pcaExplorer(dds = ddsMF_pop, dst = rldpop)

##########################################################
####### Same analysis but using the classic design: ######
##########################################################

#use cts2 without outliers

dds3 <- DESeqDataSetFromMatrix(countData = cts2, colData = coldata2, design = ~population +treatment)
#differences btw Cont vs Treat we want to account for the effect of individual, population, 
#and treatment in the analysis of the RNA-seq data. This design is 
#commonly used when we want to test the differential expression of genes 
#between two experimental conditions, while accounting for differences between 
#individuals and/or populations.

# Remove genes which have zero reads in all samples
keep <- rowSums(counts(dds3)) >= 10 #suggested in manual
dds3 <- dds3[keep,]

# in case you want to obtain normalized counts
dds3 <- estimateSizeFactors(dds3)
sizeFactors(dds2)
norm_cts_C_vs_treat <- counts(dds, normalized=TRUE)
norm_cts_log <- log2(norm_cts_C_vs_treat) #log transform normalised counts
write.csv(norm_cts_log, file="C:/Users/Joaquim Garrabou/OneDrive - The University of Hong Kong/Doctorado/ICM/Tesis/Chap2/RNA-Seq/crubrum_rnaseq/1.dds~treat/log2_normalised_counts_C-T_crubrum_10cts.csv", row.names=T)

#Set factor levels for comparison
#dds$treatment <- factor(dds$treatment, levels = c("Control", "Treatment"))
dds3$treatment <- relevel(dds3$treatment, ref = "Control")

#Differential expression
dds3 <- DESeq(dds3)#,parallel
dds3

save(dds3, file = "dds3.RData")
#load(file="dds2.RData")

#obtain results from dds object
res3 <- results(dds3, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"), test = "Wald")#obtain the results of the differential expression analysis of two groups of samples ("Treatment" and "Control") 
res3 #notice the upregulated and downregulated (+ or -, respectively) -
# the log2 fold change and Wald test p value will be for the last variable (treatment)
#in the design formula, and if this is a factor (it is), the comparison will be the 
#last level of this variable over the reference level (so, it's values for treatment)
res3 <- res3[res3$baseMean > 10, ] #remove samples with a baseMean < 10
summary(res3)

save(res3, file="res3.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/res.RData")

plotDispEsts(dds3)

#order res based on ajusted smallest p-values
res3Ordered <- res3[order(res3$padj),]
summary(res3)

######
plotMA(res3, ylim=c(-8,8)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotCounts(dds3, gene = which.min(res3$padj), intgroup = "treatment")

mcols(res2)$description #Information about which variables and tests were used 
#can be found by calling the function mcols on the results object.

#[1] "mean of normalized counts for all samples"              "log2 fold change (MLE): treatment Treatment vs Control"
#[3] "standard error: treatment Treatment vs Control"         "Wald statistic: treatment Treatment vs Control"        
#[5] "Wald test p-value: treatment Treatment vs Control"      "BH adjusted p-values" 

#extract data from res object
re3 <- data.frame(res3Ordered)
re3<- na.omit(re3)
write.csv(re3,file="deseq2_pop+treat_C-T_crubrum.csv") #This is the result of the analysis. Name it well.

#extract data which pass an adjusted p value threshold (0.05)
sig2=subset(re3,re3$padj<0.05)
sig23=subset(sig2,abs(sig2$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig23, file = "sig2_deseq2_pop+treat_C-T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig23,file="diff_genes_pop+treat_C-T_crubrum.csv") #save as cvs file


#### Data transformation and visualization ####

#So, we choose rlog
#rld <- vst(dds2,blind = FALSE) #faster transformation
rld3 <- rlog(dds3, blind = FALSE)

save(rld3, file = "rld3.RData")
#load(file = "C:/Users/Sandra/Documents/Biologia/Posgrado/PhD/Project/Chapter 2/RNA-seq/crubrum_rnaseq/1.dds~treat/rld.RData")

### Heatmap C-T ###

#heatmap of the count matrix based on rld values:
select <- order(rowMeans(counts(dds3,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds3)[,c("treatment", "Individual", "population", "day")])
pheatmap(assay(rld3)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances without clustering
sampledist <- dist(t(assay(rld3))) # calculate the euclidean distance between the samples
sampledist #gives us an overview over similarities and dissimilarities between samples
sampledistmatrix <- as.matrix(sampledist)
colnames(sampledistmatrix) <- paste(rld3$treatment,rld3$population,sep="-")
colnames(sampledistmatrix) <- NULL

pheatmap(sampledistmatrix, clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist)

#test for similarities in gene expression bwn samples. How similar the repliates are for each other
#and whether the samples from a given groups cluster each other or not
#0-1 correlation values. We expect >0.80, otherwise they may be outliers

## heatmap with gene clustering
sig23
topVarGenes <- head(order(sig23$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rld3)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

##### PCA C-T #####
head(assay(rld3))
plotPCA(rld3, intgroup="treatment")
plotPCA(rld3, intgroup=c("population"))
plotPCA(rld3, intgroup=c("treatment","population"))
plotPCA(rld3, intgroup=c("treatment","population", "Individual"))

pcaExplorer(dds = dds3, dst = rld3)

plotPCA(rld3, intgroup = "treatment",
        ntop = 500, returnData = FALSE)

dev.off()
