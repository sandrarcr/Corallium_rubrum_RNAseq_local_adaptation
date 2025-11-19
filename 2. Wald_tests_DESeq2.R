
######################Code accompanying:###################################################
#
#   Ramirez-Calero, S. et al. 2025. Does local adaptation influence thermal responses in 
#   red coral populations across depth gradients? Transcriptomic insights for effective conservation.
#
# Script written by: Sandra Ramirez
#
# This script contains all WALD TESTS used to quantify pairwise contrasts in gene expression.
#
# Wald tests are used to:
#   - Estimate direction and magnitude of differential expression with log2FoldChanges
#   - Identify specific up- or down-regulated genes under biological contrasts
#   - Compare individual conditions (e.g. Treatment vs Control, across populations and time points)
#
# Biological design:
#   2 Populations: CAS (shallow) and LOP (mesophotic)
#   2 Treatments: Control (18°C) vs Treatment (25°C)
#   3 Time points: T0, T5 and T10
#
# For this script we use the data sets already available in the repository

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

##########################Q1#########################################
# BASELINE DIFFERENCES BETWEEN POPULATIONS (CONTROL, T0) - FRONTLOADING
#
# Rationale:
#   - Before thermal stress, identify constitutive (baseline) differences between
#     shallow and mesophotic colonies under control conditions. Considered as 
#     frontloaded or pre-emptive expression (under normal temperature conditions)
#   - These genes may underlie local adaptation and differential stress tolerance.
#   - Since treatment × population interaction was not significant, we can analyze 
#     population differences within a treatment/time subset
#   - There is a strong treatment main effect and moderate treatment × time effect,
#     so it justifies comparing baseline expression separately, because time effects 
#     can be removed by restricting to T0.
#
#* Include only  control samples to have an real isolate effect to see the baseline 
#* difference between the two populations

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
nrow(sig_results)
#189

###############################Q2####################################
# GLOBAL TEMPERATURE EFFECT (all populations and days together)
#
# Rationale:
#   - The LRT showed a strong effect due to treatment
#   - Due to the limited sample size, we use the Wald test to identify which genes are up/down-regulated
#     due to temperature (all control samples vs all treatment samples) using a simpler additive model (), 
#   - This ignores small interaction effects for simplicity, as justified by the LRT.

#load data
data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)

# Design: population + treatment + day (no interactions)
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~population + treatment + day)

#check the levels of each factor
dds$population
dds$treatment
dds$day

# Keep only genes that have non-zero reads in total
keep <- rowSums(counts(dds)>= 10) >= 3
dds <- dds[keep,]

#obtain normalized counts
dds_norm <- estimateSizeFactors(dds)
sizeFactors(dds_norm)
norm_cts_C_vs_T <- counts(dds_norm, normalized=TRUE)
norm_cts_log <- log2(norm_cts_C_vs_T) #log transform normalised counts
write.csv(norm_cts_log, file="", row.names=T)

# Explicitly set the factor levels (control - reference)
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
nrow(sig_results)
#1801

################################Q3#######################################
# TEMPERATURE EFFECT WITHIN EACH POPULATION
#
# Rationale:
#   - The unified LRT confirmed a strong main effect of treatment and a negligible
#     treatment × population interaction, indicating that both populations were 
#     overall affected by temperature stress,
#   - We now analyze each population separately to see how each responds to heat.
#   - By applying a additive design (~ treatment + day) within each population,
#     we estimate the average treatment response while accounting for temporal 
#     variation, without overparameterizing the model due to small sampling size
#   - This allows us to assess local adaptation: whether shallow vs mesophotic colonies
#     activate different and specific ranscriptional pathways under the same stress.

#load data

### shallow ####
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
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

#obtain normalization counts if desired

# Explicitly set the factor levels - set control as baseline again
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
nrow(sig_results_cas)
#468

###mesophotic ####

#load data
lop_data <- read.csv("LOP_counts.csv",row.names=1)
lop_meta <- read.csv("LOP_metadata.csv",row.names=1)

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
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

#obtain normalized counts if desired

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
nrow(sig_results_lop)
#1008

################################Q4#######################################
# TIME EFFECTS OF TEMPERATURE PER POPULATION
#
# Rationale:
#   - Temporal contrasts assess how the effect of thermal stress evolves across exposure times.
#   - We compared treatment vs control samples at each day (T0, T5, T10) within each population.
#   - This approach identifies genes whose expression changes due to heat stress at early (T0), intermediate (T5),  and late time points (T10),

## load data
data <- read.csv("crubrum.gene.counts.matrix.csv",row.names=1)
meta <- read.csv("crubrum.coldata.csv",row.names=1)

data <- round(data)

(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

#defined factors
meta$treatment  <- factor(meta$treatment, levels = c("Control","Treatment"))
meta$population <- factor(meta$population, levels = c("CAS","LOP"))
meta$day        <- factor(meta$day, levels = c("0","1","2"), labels = c("T0","T5","T10"))

# Create a single factor that encodes - dummy variable for population_treatment_day
meta$factor <- factor(paste0(meta$population, "_", meta$treatment, "_", meta$day))

# Keep only genes that have non-zero reads in total
keep_genes <- rowSums(data >= 10) >= 3
data <- data[keep_genes, ]

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~ factor)

#obtain normalized counts if desired

#check the levels for each factor
dds$factor

# Set one reference
dds$factor <- relevel(dds$factor, ref = "CAS_Control_T0")

# Run DESeq
dds <- DESeq(dds)
resultsNames(dds)

#### shallow ####
#contrasts:
#### T0 (treatment T0 vs Control T0) for Shallow population
res_CAS_day0 <- results(dds, alpha = 0.05, contrast = c("factor", "CAS_Treatment_T0", "CAS_Control_T0")) #So the genes are up/ down regulated in treat day0 vs ctrl day0
res_CAS_day0 <- res_CAS_day0[ res_CAS_day0$baseMean > 10, ]
summary(res_CAS_day0)
res_CAS_day0

res_CAS_day0 <- data.frame(res_CAS_day0)
res_CAS_day0 <- na.omit(res_CAS_day0)
write.csv(res_CAS_day0, file = "")

sig_results_day0_CAS <- res_CAS_day0[res_CAS_day0$padj<0.05, ]
sig_results_day0_CAS=subset(sig_results_day0_CAS,abs(sig_results_day0_CAS$log2FoldChange)>0.3)
write.csv(sig_results_day0_CAS,file="",row.names=T)
nrow(sig_results_day0_CAS) 
#18

### T5 (treatment T5 vs Control T5) for shallow population
res_CAS_day5 <- results(dds, alpha = 0.05, contrast = c("factor", "CAS_Treatment_T5", "CAS_Control_T5")) #So the genes are up/ down regulated in treat day0 vs ctrl day0
res_CAS_day5 <- res_CAS_day5[ res_CAS_day5$baseMean > 10, ]
summary(res_CAS_day5)
res_CAS_day5

res_CAS_day5 <- data.frame(res_CAS_day5)
res_CAS_day5 <- na.omit(res_CAS_day5)
write.csv(res_CAS_day5, file = "")

sig_results_day5_CAS <- res_CAS_day5[res_CAS_day5$padj<0.05, ]
sig_results_day5_CAS=subset(sig_results_day5_CAS,abs(sig_results_day5_CAS$log2FoldChange)>0.3)
write.csv(sig_results_day5_CAS,file="",row.names=T)
nrow(sig_results_day5_CAS) 
#457

### T10 (treatment T10 vs Control T10)
res_CAS_day10 <- results(dds, alpha = 0.05, contrast = c("factor", "CAS_Treatment_T10", "CAS_Control_T10")) #So the genes are up/ down regulated in treat day0 vs ctrl day0
res_CAS_day10 <- res_CAS_day10[ res_CAS_day10$baseMean > 10, ]
summary(res_CAS_day10)
res_CAS_day10

res_CAS_day10 <- data.frame(res_CAS_day10)
res_CAS_day10 <- na.omit(res_CAS_day10)
write.csv(res_CAS_day10, file = "")

sig_results_day10_CAS <- res_CAS_day10[res_CAS_day10$padj<0.05, ]
sig_results_day10_CAS=subset(sig_results_day10_CAS,abs(sig_results_day10_CAS$log2FoldChange)>0.3)
write.csv(sig_results_day10_CAS,file="",row.names=T)
nrow(sig_results_day10_CAS) 
#490

#### mesophotic ####

#### T0 (treatment T0 vs Control T0)

res_LOP_day0 <- results(dds, alpha = 0.05, contrast = c("factor", "LOP_Treatment_T0", "LOP_Control_T0")) #So the genes are up/ down regulated in treat day0 vs ctrl day0
res_LOP_day0 <- res_LOP_day0[ res_LOP_day0$baseMean > 10, ]
summary(res_LOP_day0)
res_LOP_day0

res_LOP_day0 <- data.frame(res_LOP_day0)
res_LOP_day0 <- na.omit(res_LOP_day0)
write.csv(res_LOP_day0, file = "")

sig_results_day0_LOP <- res_LOP_day0[res_LOP_day0$padj<0.05, ]
sig_results_day0_LOP=subset(sig_results_day0_LOP,abs(sig_results_day0_LOP$log2FoldChange)>0.3)
write.csv(sig_results_day0_LOP,file="",row.names=T)
nrow(sig_results_day0_LOP) 
#13

### T5 (treatment T5 vs Control T5)
res_LOP_day5 <- results(dds, alpha = 0.05, contrast = c("factor", "LOP_Treatment_T5", "LOP_Control_T5")) #So the genes are up/ down regulated in treat day0 vs ctrl day0
res_LOP_day5 <- res_LOP_day5[ res_LOP_day5$baseMean > 10, ]
summary(res_LOP_day5)
res_LOP_day5

res_LOP_day5 <- data.frame(res_LOP_day5)
res_LOP_day5 <- na.omit(res_LOP_day5)
write.csv(res_LOP_day5, file = "")

sig_results_day5_LOP <- res_LOP_day5[res_LOP_day5$padj<0.05, ]
sig_results_day5_LOP=subset(sig_results_day5_LOP,abs(sig_results_day5_LOP$log2FoldChange)>0.3)
write.csv(sig_results_day5_LOP,file="",row.names=T)
nrow(sig_results_day5_LOP) 
#1730

### T10 (treatment T10 vs control T10)
res_LOP_day10 <- results(dds, alpha = 0.05, contrast = c("factor", "LOP_Treatment_T10", "LOP_Control_T10")) #So the genes are up/ down regulated in treat day0 vs ctrl day0
res_LOP_day10 <- res_LOP_day10[ res_LOP_day10$baseMean > 10, ]
summary(res_LOP_day10)
res_LOP_day10

res_LOP_day10 <- data.frame(res_LOP_day10)
res_LOP_day10 <- na.omit(res_LOP_day10)
write.csv(res_LOP_day10, file = "")

sig_results_day10_LOP <- res_LOP_day10[res_LOP_day10$padj<0.05, ]
sig_results_day10_LOP=subset(sig_results_day10_LOP,abs(sig_results_day10_LOP$log2FoldChange)>0.3)
write.csv(sig_results_day10_LOP,file="",row.names=T)
nrow(sig_results_day10_LOP) 
#562

########################## frontloading #########################################
# FRONTLOADING
#
# Rationale:
#   - Assess transcriptional frontloading (Barshis et al., 2013) between shallow (CAS) and
#   mesophotic (LOP) populations.
#
# Steps:
#   1) Start with genes significantly upregulated in mesophotic population (LOP) under heat stress.
#   2) Merge LOP and CAS DESeq2 results + baseline (Control T0) normalized counts.
#   3) Classify each gene into one of the four categories:
#
#        a) FRONTLOADED:
#            - LOP: significant upregulation under stress
#            - CAS: higher baseline expression (Control CAS > Control LOP)
#
#        b) REDUCED REACTION:
#            - Significant in both populations
#            - CAS log2FC < LOP log2FC
#
#        c) GREATER FOLD CHANGE:
#            - Significant in both populations
#            - CAS log2FC > LOP log2FC
#
#        d) NA:
#            - Significant in LOP but does not meet any criteria above
###

# Step 1

#load results from mesophotic under stress
res_lop <- read.csv(file = "", row.names = NULL)
rownames(res_lop) <- res_lop$SeqName
res_lop$SeqName <- NULL

#load results from shallow under stress
res_cas <- read.csv(file = "", row.names = NULL)
rownames(res_cas) <- res_cas$SeqName
res_cas$SeqName <- NULL

#select matching genes
common_genes <- intersect(rownames(res_lop), rownames(res_cas))

#step 2: generate conjunct table
df <- data.frame(
  SeqName = common_genes,
  
  # Mesophotic (LOP) Columns
  LOP_baseMean      = res_lop[common_genes, "baseMean"],
  LOP_log2FoldChange = res_lop[common_genes, "log2FoldChange"],
  LOP_padj           = res_lop[common_genes, "padj"],
  LOP_sig            = ifelse(res_lop[common_genes, "padj"] < 0.05, "sig", ""),
  
  # Shallow (CAS) Columns
  CAS_baseMean      = res_cas[common_genes, "baseMean"],
  CAS_log2FoldChange = res_cas[common_genes, "log2FoldChange"],
  CAS_padj           = res_cas[common_genes, "padj"],
  CAS_sig            = ifelse(res_cas[common_genes, "padj"] < 0.05, "sig", "")
)

#add baseline comparison — Load shallow vs mesophotic baseline (Control, T0) normalized counts you just created
data <- read.csv("",row.names=1)
meta <- read.csv("",row.names=1)

#recheck data
data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)
meta$population <- factor(meta$population, levels = c("CAS","LOP"))

dds_base <- DESeqDataSetFromMatrix(countData = data, colData = meta, 
                              design = ~population)

#extract normalized counts for baseline expression
dds_base <- estimateSizeFactors(dds_base)
norm_base <- counts(dds_base, normalized = TRUE)

### mean baseline expression in shallow and mesophotic (Control, T0)
cas_samples <- rownames(meta)[ meta$population == "CAS" ]
lop_samples <- rownames(meta)[ meta$population == "LOP" ]

mean_cas <- rowMeans(norm_base[, cas_samples])
mean_lop <- rowMeans(norm_base[, lop_samples])

#mean baseline exprssion
baseline_df <- data.frame(
  gene = rownames(norm_base),
  mean_CAS_ctrl = mean_cas,
  mean_LOP_ctrl = mean_lop
)

### Add baseline comparison (shallow vs Mesophotic at control)
df$Control_baseMean_CAS <- baseline_df$mean_CAS_ctrl[
  match(df$SeqName, baseline_df$gene)
]

df$Control_baseMean_LOP <- baseline_df$mean_LOP_ctrl[
  match(df$SeqName, baseline_df$gene)
]

df$Control_log2FoldChange <- log2( (df$Control_baseMean_CAS + 1) /
                                     (df$Control_baseMean_LOP + 1) )
#add log fold change difference
df$Control_sig <- ifelse(abs(df$Control_log2FoldChange) > 0.3, "sig", "")

#Step 3: define our gene classifications:

## Conditions
sig_cutoff <- 0.05
lfc_cutoff <- 0.3

lop_sig <- df$LOP_padj < sig_cutoff & abs(df$LOP_log2FoldChange) > lfc_cutoff
cas_sig <- df$CAS_padj < sig_cutoff & abs(df$CAS_log2FoldChange) > lfc_cutoff

## Frontloaded genes are:
frontloaded <- (df$Control_baseMean_CAS > df$Control_baseMean_LOP) & lop_sig

## Reduced reaction genes are (shallow reacts less)
reduced_reaction <- lop_sig & cas_sig & 
  (df$LOP_log2FoldChange > df$CAS_log2FoldChange)

## Greater logfold change (shallow reacts more)
greater_reaction <- lop_sig & cas_sig &
  (df$CAS_log2FoldChange > df$LOP_log2FoldChange)

## label left genes as NAs
df$Classification <- "NA"
df$Classification[greater_reaction]  <- "greater fold change"
df$Classification[reduced_reaction]  <- "reduced reaction"
df$Classification[frontloaded]       <- "frontloaded"

# Save full table:
write.csv(df,"write/full/path/and/filename.csv",
          row.names = FALSE) #filter table with upregulated mesophotic DEGs relevant for the study - add annotation and run enrichment.
