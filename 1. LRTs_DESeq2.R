####################################################################
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

#load data
data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

#check potential relevant interaction effects between factors using LRT

#treatment and time

#full model
dds_full <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                                   design = ~ treatment + day + population + treatment:day)

#reduced model - remove the interaction term
dds_reduced <- DESeq(dds_full, test="LRT", reduced = ~ treatment + population + day ) 

res_interaction_effect <- results(dds_reduced)
res_interaction_effect
summary(res_interaction_effect)
results_interaction_effect <- data.frame(res_interaction_effect)
results_interaction_effect <- na.omit(results_interaction_effect)

res_sig_interaction_effect <- results_interaction_effect[results_interaction_effect$padj < 0.05, ]
### 69 genes DE due to interaction effect


#treatment and population

#full model
dds_full <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                                   design = ~ treatment + day + population + treatment:population)

#reduced model - remove the interaction term
dds_reduced <- DESeq(dds_full, test="LRT", reduced = ~ treatment + day + population) 

res_pop_int <- results(dds_reduced)
res_pop_int
summary(res_pop_int)
results_pop_int <- data.frame(res_pop_int)
results_pop_int <- na.omit(results_pop_int)

res_sig_pop_int <- results_pop_int[results_pop_int$padj < 0.05, ]
### 0 genes DE due to pop_int effect

#Interactions don't have a major effect. We proceed by doing the analysis factor by factor without interactions

#### treatment effect #####

#load data - no need to repeat if it is already done.
data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

#full model for all downstream tests
dds_full <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                                   design = ~ treatment + day + population)

# effect of treatment - in reduced model remove treatment
dds_reduced_treat <- DESeq(dds_full, test="LRT", reduced = ~ day + population) 

res_treatment <- results(dds_reduced_treat)
res_treatment
summary(res_treatment)
results_treatment <- data.frame(res_treatment)
results_treatment <- na.omit(results_treatment)

res_sig_treatment <- results_treatment[results_treatment$padj < 0.05, ]
### 1806 genes DE due to treatment effect

##### population effect ######

#same full model 

# effect of population - in reduced model remove population
dds_reduced <- DESeq(dds_full, test="LRT", reduced = ~ treatment + day) 

res_pop <- results(dds_reduced)
res_pop
summary(res_pop)
results_pop <- data.frame(res_pop)
results_pop <- na.omit(results_pop)

res_sig_pop <- results_pop[results_pop$padj < 0.05, ]
### 3192 genes DE due to population effect

##### time effect ######

#same full model

# effect of day - in reduced model remove day
dds_reduced <- DESeq(dds_full, test="LRT", reduced = ~ treatment + population) 

res_day <- results(dds_reduced)
res_day
summary(res_day)
results_day <- data.frame(res_day)
results_day <- na.omit(results_day)

res_sig_day <- results_day[results_day$padj < 0.05, ]
### 40 genes DE due to day effect


# 3192 DEGs (pop effect), 1806 DEGs (treatment effect) and 40 DEGs (day effect) prove that population and temperature had a bigger effect in our data. Proceed for Wald tests statistics.