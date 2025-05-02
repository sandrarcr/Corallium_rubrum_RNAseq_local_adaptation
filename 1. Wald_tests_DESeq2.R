
###################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2025. Does local adaptation influence thermal responses in red coral populations across depth gradients? Transcriptomic insights for effective conservation.
#
# Script written by: Sandra Ramirez
#
#This script was run to obtain significant differentially expressed genes due to a thermal treatment stress (25C) vs a control (18C) for two contrasting populations of octocoral Corallium rubrum across three different time points (T0, T5 and T10). Data sets used under this script can be found in the associated repository
#####################################################################

setwd("~")

#load libraries

library(digest)
library(XML)
library(RSQLite)
library(DESeq2)
library("pcaExplorer")
library("ComplexHeatmap")
library("gplots")
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("dplyr")
library("ggplot2")

### Q1. Wald test to see the overall effect of temperature, regardless of the population, day, interactions
# This will show what genes are DE due to overall temp effect and also what is upregulated and what is downregulated

#load data
data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~population + treatment + day)
# to do a wald test, you cannot have the term temperature and the interaction together in the design - it gives an error

#check the levels of each factor
dds$population
dds$treatment
dds$day

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]

#obtain normalized counts
dds_norm <- estimateSizeFactors(dds)
sizeFactors(dds_norm)
norm_cts_C_vs_T <- counts(dds_norm, normalized=TRUE)
norm_cts_log <- log2(norm_cts_C_vs_T) #log transform normalised counts
write.csv(norm_cts_log, file="log2_normalised_counts_Cont-Treat.csv", row.names=T)

# Explicitly set the factor levels
dds$treatment <- relevel(dds$treatment, ref = "Control")

dds <- DESeq(dds) #run DEseq2
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("treatment","Treatment","Control"), test ="Wald")
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
write.csv(results,file="") #total results

#filter results
sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results, file="")

#**********************************************************************
#Q2. What is the overall effect of temperature for each population?? To make it easier, first subset the data and do the analysis separately for each population

#load data

#shallow
cas_data <- read.csv("CAS_counts.csv",row.names=1)
cas_meta <- read.csv("CAS_metadata.csv",row.names=1)

cas_data <- round(cas_data)

(all(rownames(cas_meta) %in% colnames(cas_data)) || all(colnames(cas_data) %in% rownames(cas_meta)))
(all(colnames(cas_data) == rownames(cas_meta)))

cas_meta$day <- factor(cas_meta$day)

dds <- DESeqDataSetFromMatrix(countData = cas_data, colData = cas_meta, 
                              design = ~ treatment + day)

#check the levels for each factor
dds$treatment
dds$day

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$treatment <- relevel(dds$treatment, ref = "Control")

dds <- DESeq(dds)
resultsNames(dds)

res_cas <- results(dds, alpha = 0.05, contrast = c("treatment","Treatment","Control"), test ="Wald")
res_cas <- res_cas[res_cas$baseMean > 10, ]

summary(res_cas)
res_cas

results_cas <- data.frame(res_cas)
results_cas <- na.omit(results_cas)
write.csv(results_cas,file="") #total results


sig_results_cas <- results_cas[results_cas$padj<0.05, ]
sig_results_cas=subset(sig_results_cas,abs(sig_results_cas$log2FoldChange)>0.3)
write.csv(sig_results_cas,file="",row.names=T)

#mesophotic
lop_data <- read.csv("LOP_counts.csv",row.names=1)
lop_meta <- read.csv("LOP_metadata.csv",row.names=1)

#remove outlier LOP_C0_1

lop_data <- lop_data[,-4]
lop_meta <- lop_meta[-4, ]

lop_data <- round(lop_data)

(all(rownames(lop_meta) %in% colnames(lop_data)) || all(colnames(lop_data) %in% rownames(lop_meta)))
(all(colnames(lop_data) == rownames(lop_meta)))

lop_meta$day <- factor(lop_meta$day)

dds <- DESeqDataSetFromMatrix(countData = lop_data, colData = lop_meta, 
                              design = ~ treatment + day)

#check the levels for each factor
dds$treatment
dds$day

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$treatment <- relevel(dds$treatment, ref = "Control")

dds <- DESeq(dds)
resultsNames(dds)

res_lop <- results(dds, alpha = 0.05, contrast = c("treatment","Treatment","Control"), test ="Wald")
res_lop <- res_lop[res_lop$baseMean > 10, ]
summary(res_lop)
res_lop

results_lop <- data.frame(res_lop)
results_lop <- na.omit(results_lop)
write.csv(results_lop,file="") #total results

sig_results_lop <- results_lop[results_lop$padj<0.05, ]
sig_results_lop=subset(sig_results_lop,abs(sig_results_lop$log2FoldChange)>0.3)
write.csv(sig_results_lop,file="",row.names=T)

#*****************************************************************************************************************
#Q3. How does temperature affects gene expression for each day - for each population

## load data for shallow
cas_data <- read.csv("CAS_counts.csv",row.names=1)
cas_meta <- read.csv("CAS_metadata.csv",row.names=1)

cas_data <- round(cas_data)

(all(rownames(cas_meta) %in% colnames(cas_data)) || all(colnames(cas_data) %in% rownames(cas_meta)))
(all(colnames(cas_data) == rownames(cas_meta)))

cas_meta$day <- factor(cas_meta$day)

dds <- DESeqDataSetFromMatrix(countData = cas_data, colData = cas_meta, 
                              design = ~ factor)

#check the levels for each factor
dds$factor

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#### T0

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "ctrl_day0")

dds <- DESeq(dds)
resultsNames(dds)

#pairwise for treat day0 vs ctrl day0
res_day0 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day0","ctrl_day0"), test ="Wald")
#So the genes are up/ down regulated in treat day0 vs ctrl day0
res_day0 <- res_day0[res_day0$baseMean > 10, ]
summary(res_day0)
res_day0

results_day0 <- data.frame(res_day0)
results_day0 <- na.omit(results_day0)

sig_results_day0 <- results_day0[results_day0$padj<0.05, ]
sig_results_day0=subset(sig_results_day0,abs(sig_results_day0$log2FoldChange)>0.3)

### T5

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "ctrl_day1")

#run dds again
dds <- DESeq(dds)
resultsNames(dds)

#pairwise for treat day1 vs ctrl day1
res_day1 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day1","ctrl_day1"), test ="Wald")
#So the genes are up/ down regulated in treat day1 vs ctrl day1
res_day1 <- res_day1[res_day1$baseMean > 10, ]
summary(res_day1)
res_day1

results_day1 <- data.frame(res_day1)
results_day1 <- na.omit(results_day1)
write.csv(results_day1,file="") #total results

sig_results_day1 <- results_day1[results_day1$padj<0.05, ]
sig_results_day1=subset(sig_results_day1,abs(sig_results_day1$log2FoldChange)>0.3)
write.csv(sig_results_day1,file="",row.names=T)

### T10

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "ctrl_day2")

#run dds again
dds <- DESeq(dds)
resultsNames(dds)

#pairwise for treat day2 vs ctrl day2
res_day2 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","ctrl_day2"), test ="Wald")
#So the genes are up/ down regulated in treat day2 vs ctrl day2
res_day2 <- res_day2[res_day2$baseMean > 10, ]
summary(res_day2)
res_day2

#filter results
results_day2 <- data.frame(res_day2)
results_day2 <- na.omit(results_day2)
write.csv(results_day2,file="") #total results

sig_results_day2 <- results_day2[results_day2$padj<0.05, ]
sig_results_day2=subset(sig_results_day2,abs(sig_results_day2$log2FoldChange)>0.3)
write.csv(sig_results_day2,file="",row.names=T)

#only for treatments

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "treat_day0")

dds <- DESeq(dds)
resultsNames(dds)

#pairwise for day 1 vs 0 (in the Treatment)
res_treat_day1 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day1","treat_day0"), 
                          test ="Wald")
res_treat_day1 <- res_treat_day1[res_treat_day1$baseMean > 10, ]
#So the genes are up/ down regulated in treatment day1 vs treatment day0

summary(res_treat_day1)
res_treat_day1

results_treat_day1 <- data.frame(res_treat_day1)
results_treat_day1 <- na.omit(results_treat_day1)
write.csv(results_treat_day1,file="") #total results

sig_results_treat_day1 <- results_treat_day1[results_treat_day1$padj<0.05, ]
sig_results_treat_day1=subset(sig_results_treat_day1,abs(sig_results_treat_day1$log2FoldChange)>0.3)
write.csv(sig_results_treat_day1,file="")

#pairwise for day 2 vs 0 (in the Treatment)
res_treat_day2 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","treat_day0"), 
                          test ="Wald")
#So the genes are up/ down regulated in treatment day2 vs treatment day0
res_treat_day2 <- res_treat_day2[res_treat_day2$baseMean > 10, ]

summary(res_treat_day2)
res_treat_day2

results_treat_day2 <- data.frame(res_treat_day2)
results_treat_day2 <- na.omit(results_treat_day2)
write.csv(results_treat_day2,file="") #total results

sig_results_treat_day2 <- results_treat_day2[results_treat_day2$padj<0.05, ]
sig_results_treat_day2=subset(sig_results_treat_day2,abs(sig_results_treat_day2$log2FoldChange)>0.3)
write.csv(sig_results_treat_day2,file="")

#pairwise for day 2 vs day1 (in the Treatment)

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "treat_day1")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","treat_day1"), test ="Wald")
#So the genes are up/ down regulated in treatment day2 vs treatment day1
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
write.csv(results,file="") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results,file="")

#load data for mesophotic

lop_data <- read.csv("LOP_counts.csv",row.names=1)
lop_meta <- read.csv("LOP_metadata.csv",row.names=1)

#remove outlier LOP_C0_1

lop_data <- lop_data[,-4]
lop_meta <- lop_meta[-4, ]

lop_data <- round(lop_data)

(all(rownames(lop_meta) %in% colnames(lop_data)) || all(colnames(lop_data) %in% rownames(lop_meta)))
(all(colnames(lop_data) == rownames(lop_meta)))

lop_meta$day <- factor(lop_meta$day)

dds <- DESeqDataSetFromMatrix(countData = lop_data, colData = lop_meta, 
                              design = ~ factor)

#check the levels for each factor
dds$factor

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

#### T0

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "ctrl_day0")

dds <- DESeq(dds)
resultsNames(dds)

#pairwise for treat day0 vs ctrl day0
res_day0 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day0","ctrl_day0"), test ="Wald")
#So the genes are up/ down regulated in treat day0 vs ctrl day0
res_day0 <- res_day0[res_day0$baseMean > 10, ]

summary(res_day0)
res_day0

results_day0 <- data.frame(res_day0)
results_day0 <- na.omit(results_day0)
write.csv(results_day0,file="") #total results

sig_results_day0 <- results_day0[results_day0$padj<0.05, ]
sig_results_day0=subset(sig_results_day0,abs(sig_results_day0$log2FoldChange)>0.3)
write.csv(sig_results_day0,file="",row.names=T)

### T5

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "ctrl_day1")

#run dds again
dds <- DESeq(dds)
resultsNames(dds)

#pairwise for treat day1 vs ctrl day1
res_day1 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day1","ctrl_day1"), test ="Wald")
#So the genes are up/ down regulated in treat day1 vs ctrl day1
res_day1 <- res_day1[res_day1$baseMean > 10, ]

summary(res_day1)
res_day1

results_day1 <- data.frame(res_day1)
results_day1 <- na.omit(results_day1)
write.csv(results_day1,file="") #total results

sig_results_day1 <- results_day1[results_day1$padj<0.05, ]
sig_results_day1=subset(sig_results_day1,abs(sig_results_day1$log2FoldChange)>0.3)
write.csv(sig_results_day1,file="",row.names=T)

### T10

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "ctrl_day2")

#run dds again
dds <- DESeq(dds)
resultsNames(dds)

#pairwise for treat day2 vs ctrl day2
res_day2 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","ctrl_day2"), test ="Wald")
#So the genes are up/ down regulated in treat day2 vs ctrl day2
res_day2 <- res_day2[res_day2$baseMean > 10, ]

summary(res_day2)
res_day2

results_day2 <- data.frame(res_day2)
results_day2 <- na.omit(results_day2)
write.csv(results_day2,file="") #total results

sig_results_day2 <- results_day2[results_day2$padj<0.05, ]
sig_results_day2=subset(sig_results_day2,abs(sig_results_day2$log2FoldChange)>0.3)
write.csv(sig_results_day2,file="",row.names=T)

#only for treatments

# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "treat_day0")

dds <- DESeq(dds)
resultsNames(dds)

#pairwise for day 1 vs 0 (in the Treatment)
res_treat_day1 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day1","treat_day0"), test ="Wald")
#So the genes are up/ down regulated in treatment day1 vs treatment day0
res_treat_day1 <- res_treat_day1[res_treat_day1$baseMean > 10, ]

summary(res_treat_day1)
res_treat_day1

results_treat_day1 <- data.frame(res_treat_day1)
results_treat_day1 <- na.omit(results_treat_day1)
write.csv(results_treat_day1,file="") #total results

sig_results_treat_day1 <- results_treat_day1[results_treat_day1$padj<0.05, ]
sig_results_treat_day1=subset(sig_results_treat_day1,abs(sig_results_treat_day1$log2FoldChange)>0.3)
write.csv(sig_results_treat_day1,file="")

#pairwise for day 2 vs 0 (in the Treatment)
res_treat_day2 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","treat_day0"), test ="Wald")
#So the genes are up/ down regulated in treatment day2 vs treatment day0
res_treat_day2 <- res_treat_day2[res_treat_day2$baseMean > 10, ]

summary(res_treat_day2)
res_treat_day2

results_treat_day2 <- data.frame(res_treat_day2)
results_treat_day2 <- na.omit(results_treat_day2)
write.csv(results_treat_day2,file="") #total results

sig_results_treat_day2 <- results_treat_day2[results_treat_day2$padj<0.05, ]
sig_results_treat_day2=subset(sig_results_treat_day2,abs(sig_results_treat_day2$log2FoldChange)>0.3)
write.csv(sig_results_treat_day2,file="")

#pairwise for day 2 vs day1 (in the Treatment)

#first need to change the reference level
# Explicitly set the factor levels - the reference level
dds$treatment <- relevel(dds$factor, ref = "treat_day1")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","treat_day1"), test ="Wald")
#So the genes are up/ down regulated in treatment day2 vs treatment day1
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
write.csv(results,file="") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results,file="")

#*************************************************************************************
#* Q4. what is the differnece between the two populations at time 0
#* only for control to see the baseline difference between the two populations

#load data
data <- read.csv("counts_cas_lop_ctrl_day0.csv",row.names=1)
meta <- read.csv("coldata_cas_lop_ctrl_day0.csv",row.names=1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~population)

#check the levels of each factor
dds$population


# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$population <- relevel(dds$population, ref = "CAS")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("population","CAS","LOP"), test ="Wald") #L2FC> 0 upregulated in CAS
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)

write.csv(results,file="") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results,file="")

########## Q4 Ctrl vs treat per each pop at day T0 #####

#shallow
data <- read.csv("counts_cas_ctrl_T_day0.csv",row.names=1)
meta <- read.csv("coldata_cas_ctrl_T_day0.csv",row.names=1)


data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$treatment <- factor(meta$treatment)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~treatment)

#check the levels of each factor
dds$treatment

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$treatment <- relevel(dds$treatment, ref = "Control")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("treatment", "Treatment","Control"), test ="Wald")
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
write.csv(results,file="") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results,file="")

#mesophotic
data <- read.csv("counts_lop_ctrl_T_day0.csv",row.names=1)
meta <- read.csv("coldata_lop_ctrl_T_day0.csv",row.names=1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$treatment <- factor(meta$treatment)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~treatment)

#check the levels of each factor
dds$treatment

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Explicitly set the factor levels
dds$treatment <- relevel(dds$treatment, ref = "Control")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("treatment", "Treatment","Control"), test ="Wald")
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
write.csv(results,file="") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results,file="")