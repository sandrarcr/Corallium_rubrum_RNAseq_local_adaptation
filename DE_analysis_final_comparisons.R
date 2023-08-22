#First 3 main questions:

# Q1. overall effect of temp regardless of the population or day
# Q2. overall effect of temp regardless of day for each population - what is common and what what is 
#unique to each population?
# Q3. How is treatment affecting gene expression at each time point/ day i.e is there a difference in 
#gene expression from one day to the other? - separately for each population

# All genes affected by temp - so either the results from the LRT (with 1837 DE genes) or from the Wald test (with 2052 DE genes)
# Do all the remaining analysis for each population separately as there is a population effect - so basically all the pairwise comparisons (overall treat vs ctrl regardless of the day; treat day0 vs ctrl day0; treat day1 vs ctrl day1; treat day2 vs ctrl day2)
# And then you can compare the comparisons between the two populations. 

#load packages
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

setwd("path/directory/")

data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

#remove detected outliers? //LOP_C0_1

data <- data[,-22]
meta <- meta[-22, ]

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

#First check if there is actually an interaction effect

dds_full <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                                   design = ~ treatment + day + population + treatment:day)

# Remove genes which have zero reads in all samples
keep <- rowSums(counts(dds_full)) >= 10 #suggested in manual
dds_full <- dds_full[keep,]

#reduced model - remove the interaction term
dds_reduced <- DESeq(dds_full, test="LRT", reduced = ~ treatment + population + day ) 

res_interaction_effect <- results(dds_reduced)
res_interaction_effect
res_interaction_effect <- res_interaction_effect[res_interaction_effect$baseMean > 10, ]
summary(res_interaction_effect)

save(res_interaction_effect, file="/path/directory/filename.RData")#save object

results_interaction_effect <- data.frame(res_interaction_effect)
results_interaction_effect <- na.omit(results_interaction_effect)

write.csv(results_interaction_effect,file="filename.csv") #total results

#filter res by significance
res_sig_interaction_effect <- results_interaction_effect[results_interaction_effect$padj < 0.05, ]
write.csv(res_sig_interaction_effect,file="diff_siggenes_interaction_LRT_effect_crubrum.csv")#save as cvs
#129 with applied filters and outliers removed

### how many DE genes due to temp effect? - so remove all temp terms in the design - this would the 
#both the treatment and the interaction terms

dds_reduced_tmp <- DESeq(dds_full, test="LRT", reduced = ~ population + day ) 

res_tmp <- results(dds_reduced_tmp)
res_tmp
res_tmp <- res_tmp[res_tmp$baseMean > 10, ]
summary(res_tmp)

save(res_tmp, file="/path/directory/filename.RData")#save object

results_tmp <- data.frame(res_tmp)
results_tmp <- na.omit(results_tmp)

write.csv(results_tmp,file="/path/directory/filename.csv") #total results

#filter res by significance 
res_sig_tmp <- results_tmp[results_tmp$padj < 0.05, ]
write.csv(res_sig_tmp,file="filename.csv")#save as cvs
#1781 with applied filters and outliers removed

# So from testing the interaction effect and temperature effect, it seems temperature definitely has a 
#very strong effect.
#* There is a small effect of the interaction between temperature and day (but not very strong - 129 genes 
#due to interaction vs 1781 genes due to temperature)

### Just to check - what happens if the full model is ~ population + day + treatment and the reduced 
#model is ~ population + day or ~population + treatment??? 
### so assuming there is no interaction just checking if day or treatment has a larger effect

data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

#remove outlier LOP_C0_1

data <- data[,-22]
meta <- meta[-22, ]

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

dds_full <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                                   design = ~ treatment + day + population)

# effect of day - in reduced model remove day
dds_reduced <- DESeq(dds_full, test="LRT", reduced = ~ treatment + population) 

res_day <- results(dds_reduced)
res_day
res_day <- res_day[res_day$baseMean > 10, ]
summary(res_day)

save(res_day, file="/path/directory/filename.RData")#save object

results_day <- data.frame(res_day)
results_day <- na.omit(results_day)

write.csv(results_day,file="filename.csv") #total results

#filter res by significance 
res_sig_day <- results_day[results_day$padj < 0.05, ]
write.csv(res_sig_day, file="filename.csv")#save as cvs
#51 genes with all applied filters and outlier removed

# effect of treatment - in reduced model remove treatment
dds_reduced_treat <- DESeq(dds_full, test="LRT", reduced = ~ day + population) 

res_treatment <- results(dds_reduced_treat)
res_treatment
res_treatment <- res_treatment[res_treatment$baseMean > 10, ]
summary(res_treatment)

save(res_treatment, file="filename.RData")#save object

results_treatment <- data.frame(res_treatment)
results_treatment <- na.omit(results_treatment)
write.csv(results_treatment, file="/path/directory/filename.csv")#save as cvs

#filter res by significance 
res_sig_treatment <- results_treatment[results_treatment$padj < 0.05, ]
write.csv(res_sig_treatment, file="filename.csv")#save as cv
#1794 degs

#*so this confirms that actually the effect of "day" is not too much - only 51 genes vs 1794 due to 
#*treatment. So this could also explain the small effect of interaction between day and treatment - day 
#*alone has a small effect so the interaction has a small effect

#*******************************************************************************************************************
### Q1 
# Next - do a wald test to see the overall effect of temperature, regardless of the population, day, 
#interaction. This will show what genes are DE due to overall temp effect and also what is upregulated and what is 
#downregulated

data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

#remove outlier LOP_C0_1

data <- data[,-22]
meta <- meta[-22, ]

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~population + treatment + day)

# to do a wald test, you cannot have the term temperature and the interaction together in the design - it 
#gives an error
# so here we ignore the interaction effect, and just look at the overall effect of the temperature on the 
#data
# This is why getting the overall temp effect using LRT is more accurate/ correct cuz there you can 
#include the interaction effect, but then you cannot differentiate up or down regulated genes

#check the levels of each factor
dds$population
dds$treatment
dds$day

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

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
write.csv(results,file="/path/directory/filename.csv") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
write.csv(sig_results, file="filename.csv")#save as cvs
#1957

#So this answers the first question Q1 - what is the overall effect of temperature regardless of the day 
#and population

#**********************************************************************
#Next Q2. What is the overall effect of temperature for each population?? 
#*To do this you first need to subset the data and do the analysis separately for each population: CAS and LOP

#CAS

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
write.csv(results_cas,file="filename.csv") #total results

sig_results_cas <- results_cas[results_cas$padj<0.05, ]
sig_results_cas=subset(sig_results_cas,abs(sig_results_cas$log2FoldChange)>0.3)
#441 DE after filtering for basemean and lo2fc -- no outliers

write.csv(sig_results_cas,file="filename.csv",row.names=T)

### now do the same for the LOP population

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
write.csv(results_lop,file="filename.csv") #total results

sig_results_lop <- results_lop[results_lop$padj<0.05, ]
sig_results_lop=subset(sig_results_lop,abs(sig_results_lop$log2FoldChange)>0.3)
#1263

write.csv(sig_results_lop,file="filename.csv",row.names=T)

#### Now you can compare the DE results between cas and lop and see what genes are common to each 
#population and what genes are unique to each population (e.g. Venn diagrams).

#*Now moving to Q3. How does temp affects gene expression for each day - for each population

#CAS

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

#### First for day0

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

#No file to save as there are 0 genes

#0 DE genes - I think this makes sense, if for treat day0, they were sampled as soon as the temp is set 
#to 25, then they have had no time to respond to the increase in temp. i.e they are basically the same as 
#control
# This also means that any difference you see in DE genes numbers while doing this comparison but keeping 
#all the samples from both populations is actually due to the inherent difference between the two
#populaitons - something that might be interesting to look at too

### Next treat vs ctrl for day 1 - so we need to change the reference level to ctrl day1 first

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
write.csv(results_day1,file="filename.csv") #total results

sig_results_day1 <- results_day1[results_day1$padj<0.05, ]
sig_results_day1=subset(sig_results_day1,abs(sig_results_day1$log2FoldChange)>0.3)
#348 DE genes - need to filter by basemean and log2fc
#348 once it's filter

write.csv(sig_results_day1,file="filename.csv",row.names=T)

### Next treat vs ctrl for day 2 - so we need to change the reference level to ctrl day2 first

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
write.csv(results_day2,file="filename.csv") #total results

sig_results_day2 <- results_day2[results_day2$padj<0.05, ]
sig_results_day2=subset(sig_results_day2,abs(sig_results_day2$log2FoldChange)>0.3)

write.csv(sig_results_day2,file="filename.csv",row.names=T)
#354 DE

# Next we can look at how the gene expression changes by day within each treatment
# So, we would do this only for the samples at 25 (so the treatment), not for the control. 
#Cuz biologically speaking it doesn't make sense why control day 1 or day 2 should be different from 
#control day 0. Unless there is a specific reason you expect them to be different

# So, doing for the treatment in CAS
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
write.csv(results_treat_day1,file="filename.csv") #total results

sig_results_treat_day1 <- results_treat_day1[results_treat_day1$padj<0.05, ]
sig_results_treat_day1=subset(sig_results_treat_day1,abs(sig_results_treat_day1$log2FoldChange)>0.3)
# 629 DE genes

write.csv(sig_results_treat_day1,file="filename.csv.csv")

#pairwise for day 2 vs 0 (in the Treatment)
res_treat_day2 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","treat_day0"), 
                          test ="Wald")
#So the genes are up/ down regulated in treatment day2 vs treatment day0
res_treat_day2 <- res_treat_day2[res_treat_day2$baseMean > 10, ]

summary(res_treat_day2)
res_treat_day2

results_treat_day2 <- data.frame(res_treat_day2)
results_treat_day2 <- na.omit(results_treat_day2)
write.csv(results_treat_day2,file="filename.csv") #total results

sig_results_treat_day2 <- results_treat_day2[results_treat_day2$padj<0.05, ]
sig_results_treat_day2=subset(sig_results_treat_day2,abs(sig_results_treat_day2$log2FoldChange)>0.3)

write.csv(sig_results_treat_day2,file="filename.csv")
#246 DE 

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
write.csv(results,file="filename.csv") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)

#78 DE genes

write.csv(sig_results,file="filename.csv")

## do all the same as above but for LOP

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

#### First for day0

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
write.csv(results_day0,file="filename.csv") #total results

sig_results_day0 <- results_day0[results_day0$padj<0.05, ]
sig_results_day0=subset(sig_results_day0,abs(sig_results_day0$log2FoldChange)>0.3)

write.csv(sig_results_day0,file="filename.csv",row.names=T)

#20 DE genes - I think this makes sense, if for treat day0, they were sampled as soon as the temp is set 
#to 25, then they have had no time to respond to the increase in temp. i.e they are basically the same 
#as control
# This also means that any difference you see in DE genes numbers while doing this comparison but keeping 
#all the samples from both populations is actually due to the inherent difference between the two 
#populaitons - something that might be interesting to look at too

### Next treat vs ctrl for day 1 - so we need to change the reference level to ctrl day1 first

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
write.csv(results_day1,file="filename.csv") #total results

sig_results_day1 <- results_day1[results_day1$padj<0.05, ]
sig_results_day1=subset(sig_results_day1,abs(sig_results_day1$log2FoldChange)>0.3)

write.csv(sig_results_day1,file="filename.csv",row.names=T)
#1874 DE

### Next treat vs ctrl for day 2 - so we need to change the reference level to ctrl day2 first

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
write.csv(results_day2,file="filename.csv") #total results

sig_results_day2 <- results_day2[results_day2$padj<0.05, ]
sig_results_day2=subset(sig_results_day2,abs(sig_results_day2$log2FoldChange)>0.3)

#686 DE genes - need to filter by basemean and log2fc

write.csv(sig_results_day2,file="filename.csv",row.names=T)
#709 DE

# Next we can look at how the gene expression changes by day within each treatment
# So, we would do this only for the samples at 25 (so the treatment), not for the control. 
#Cuz biologically speaking it doesn't make sense why control day 1 or day 2 should be different 
#from control day 0. Unless there is a specific reason you expect them to be different

# So, doing for the treatment
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
write.csv(results_treat_day1,file="filename.csv") #total results

sig_results_treat_day1 <- results_treat_day1[results_treat_day1$padj<0.05, ]
sig_results_treat_day1=subset(sig_results_treat_day1,abs(sig_results_treat_day1$log2FoldChange)>0.3)

write.csv(sig_results_treat_day1,file="filename.csv")
#1662 DE

#pairwise for day 2 vs 0 (in the Treatment)
res_treat_day2 <- results(dds, alpha = 0.05, contrast = c("factor","treat_day2","treat_day0"), test ="Wald")
#So the genes are up/ down regulated in treatment day2 vs treatment day0
res_treat_day2 <- res_treat_day2[res_treat_day2$baseMean > 10, ]

summary(res_treat_day2)
res_treat_day2

results_treat_day2 <- data.frame(res_treat_day2)
results_treat_day2 <- na.omit(results_treat_day2)
write.csv(results_treat_day2,file="filename.csv") #total results

sig_results_treat_day2 <- results_treat_day2[results_treat_day2$padj<0.05, ]
sig_results_treat_day2=subset(sig_results_treat_day2,abs(sig_results_treat_day2$log2FoldChange)>0.3)
# 700 DE

write.csv(sig_results_treat_day2,file="filename.csv")

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
write.csv(results,file="filename.csv") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)
#105 DE

write.csv(sig_results,file="filename.csv")

#* TO see what is the differnece between the two populations at time 0
#* we would do this only for control - cuz the idea is to see the baseline difference between the two 
#* populations

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

res <- results(dds, alpha = 0.05, contrast = c("population","LOP","CAS"), test ="Wald")
res <- res[res$baseMean > 10, ]

summary(res)
res

results <- data.frame(res)
results <- na.omit(results)
write.csv(results,file="filename.csv") #total results

sig_results <- results[results$padj<0.05, ]
sig_results=subset(sig_results,abs(sig_results$log2FoldChange)>0.3)

write.csv(sig_results,file="filename.csv")
#189 DE

#Now, proceed to merge and functional enrichment :)
