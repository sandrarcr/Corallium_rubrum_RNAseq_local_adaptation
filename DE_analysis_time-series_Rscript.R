################################### Time-series ####################################

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

#day: This variable represents the time points or days at which the experiment was conducted. 
#It suggests that you have measured gene expression at different time points and want to account 
#for any time-dependent effects on gene expression.
#treatment: This variable represents the treatment conditions or experimental groups to which your samples belong. 
#It suggests that you want to compare the gene expression between different treatment groups.
#day:treatment: This component represents the interaction between the day and treatment variables. 
#It allows you to assess whether the effect of the treatment on gene expression varies depending on the time point. 
#In other words, it allows you to identify genes that show different responses to treatment at different time points.

dds_time <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata1, design = ~day + treatment + day:treatment) #full model

#factors
dds_time$day <- factor(dds_time$day, levels = c("0","1","2"))

# Remove genes which have at least 10 reads in all samples
keep <- rowSums(counts(dds_time)) >= 10
dds_time <- dds_time[keep,]

######## Run Differential expression Analysis
dds_LRT <- DESeq(dds_time, test="LRT", reduced = ~ day + treatment) #reduced --The LRT examines two 
#models for the counts, a full model with a certain number of terms and a reduced model,
#in which some of the terms of the full model are removed. The test determines if the 
#increased likelihood of the data using the extra terms in the full model is more than 
#expected if those extra terms are truly zero. The LRT is therefore useful for testing 
#multiple terms at once, for example testing 3 or more levels of a factor at once, or
#all interactions between two variables. o

save(dds_LRT, file = "dds_lrt_time.RData")

resultsNames(dds_LRT)
#[1] "Intercept"                      "day_1_vs_0"                     "day_2_vs_0"                    
#[4] "treatment_Treatment_vs_Control" "day1.treatmentTreatment"        "day2.treatmentTreatment"     

#obtain results from dds object
resT1_treat <- results(dds_LRT, alpha = 0.05, name="day1.treatmentTreatment", test="Wald") #the effect of treatment on gene expression is different at day 1 compared to the reference level (day 0).
resT2_treat <- results(dds_LRT, alpha = 0.05, name = "day2.treatmentTreatment", test="Wald") #the effect of treatment on gene expression is different at day 2 compared to the reference level (day 0)
res_time <-  results(dds_LRT, alpha = 0.05, name = "treatment_Treatment_vs_Control", test="Wald")#compares the treatment condition (Treatment) to the control condition (Control), while keeping the day variable constant.
res_intercept <-  results(dds_LRT, alpha = 0.05, name = "Intercept", test="Wald")#compares the treatment condition (Treatment) to the control condition (Control), while keeping the day variable constant.
res_t1 <- results(dds_LRT, alpha = 0.05, name = "day_1_vs_0", test="Wald")
rest_t2 <- results(dds_LRT, alpha = 0.05, name = "day_2_vs_0", test="Wald")

summary(resT1_treat) 
summary(resT2_treat) 

save(resT1_treat, file = "resT1_treat.RData")
save(resT2_treat, file = "resT2_treat.RData")

plotDispEsts(dds_LRT)

#order res based on ajusted smallest p-values
resOrderedT1 <- resT1_treat[order(resT1_treat$padj),]
summary(resT1_treat)

resOrderedT2 <- resT2_treat[order(resT2_treat$padj),]
summary(resT2_treat)

######
plotMA(resT1_treat, ylim=c(-8,8)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotMA(resT2_treat, ylim=c(-8,8)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotCounts(dds_LRT, gene = which.min(resT1_treat$padj), intgroup = "treatment")
plotCounts(dds_LRT, gene = which.min(resT2_treat$padj), intgroup = "treatment")

#extract data from res object
re1 <- data.frame(resOrderedT1)
re1 <- na.omit(re1)
write.csv(re1,file="deseq2_time1_T_crubrum.csv")

#extract data from res object
re2 <- data.frame(resOrderedT2)
re2 <- na.omit(re2)
write.csv(re2,file="deseq2_time2_T_crubrum.csv")

#extract data which pass an adjusted p value threshold (0.05)
sig1=subset(re1,re1$padj<0.05)
sig21=subset(sig1,abs(sig1$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig21, file = "sig2_deseq2_time1-T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig21,file="diff_genes_time1-T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig2=subset(re2,re2$padj<0.05)
sig22=subset(sig2,abs(sig2$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig22, file = "sig2_deseq2_time2-T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig22,file="diff_genes_time2-T_crubrum.csv") #save as cvs file

#visualization

rld <- rlog(dds_LRT, blind = FALSE)
save(rld, file = "rld_timeseries.RData")

######################## Heatmap C-T ###########################

#heatmap of the count matrix based on rld values:

select <- order(rowMeans(counts(dds_LRT,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_LRT)[,c("treatment","day")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap of the sample-to-sample distances without clustering
sampledist <- dist(t(assay(rld))) # calculate the euclidean distance between the samples
sampledist #gives us an overview over similarities and dissimilarities between samples
sampledistmatrix <- as.matrix(sampledist)
colnames(sampledistmatrix) <- paste(rld$treatment,rld$day,sep="-")
colnames(sampledistmatrix) <- NULL

pheatmap(sampledistmatrix, clustering_distance_rows=sampledist,
         clustering_distance_cols=sampledist)

## heatmap with gene clustering
sig21
topVarGenes <- head(order( sig21$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rld)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

## heatmap with gene clustering
sig22
topVarGenes <- head(order( sig22$padj ), decreasing=TRUE , 100)
heatmap.2(assay(rld)[ topVarGenes, ], 
          Colv = FALSE,  scale="row",labRow=TRUE,
          trace="none", col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(800))

dev.off()

##### PCA C-T #####

head(assay(rld))
plotPCA(rld, intgroup="treatment")
plotPCA(rld, intgroup=c("day"))
plotPCA(rld, intgroup=c("treatment","day"))

pcaExplorer(dds = dds_LRT, dst = rld)

plotPCA(rld, intgroup = "treatment",
        ntop = 500, returnData = FALSE)

#model counts for the groups over time 
cru <- plotCounts(dds_LRT, which.min(res_t1$padj), 
                  intgroup = c("day","treatment"), returnData = TRUE)
cru$day <- as.numeric(as.character(cru$day))
ggplot(cru,
       aes(x = day, y = count, color = treatment, group = treatment)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()


#model counts for the groups over time 
cru <- plotCounts(dds_LRT, which.min(res_t2$padj), 
                  intgroup = c("day","treatment"), returnData = TRUE)
cru$day <- as.numeric(as.character(cru$day))
ggplot(cru,
       aes(x = day, y = count, color = treatment, group = treatment)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()

