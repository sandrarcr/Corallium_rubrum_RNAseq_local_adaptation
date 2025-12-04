###################### Code accompanying:###################################################
#
#   Ramirez-Calero, S. et al. 2025. Does local adaptation influence thermal responses in
#   red coral populations across depth gradients? Transcriptomic insights for effective conservation.
#
# Script written by: Sandra Ramirez
#
# This script must be run within a Rproject containg the data folder and the renv.lock file.
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

# If you didn't run LRT*.R script before,
# then run this to retrieve the project library:
renv::restore()

# load libraries

library(digest)
library(XML)
library(RSQLite)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

########################## Q1#########################################
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

data <- read.csv("data/counts_cas_lop_ctrl_day0.csv", row.names = 1)
meta <- read.csv("data/coldata_cas_lop_ctrl_day0.csv", row.names = 1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)

dds <- DESeqDataSetFromMatrix(
  countData = data, colData = meta,
  design = ~population
)

# check the levels of each factor
dds$population

# Keep only genes that have non-zero reads in total
cat("Raw genes", nrow(counts(dds)), "genes before filtering")
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep, ]
cat("Kept genes", nrow(counts(dds)), "genes after filtering")

# Explicitly set the factor levels
dds$population <- relevel(dds$population, ref = "CAS")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, contrast = c("population", "CAS", "LOP"), test = "Wald") # L2FC > 0 upregulated in CAS
res <- res[res$baseMean > 10, ]

summary(res)

results <- data.frame(res)
results <- na.omit(results)
write.csv(results, file = "") # total results

sig_results <- results[results$padj < 0.05, ]
sig_results <- subset(sig_results, abs(sig_results$log2FoldChange) > 0.3)
write.csv(sig_results, file = "")
nrow(sig_results)
# 292

############################### Q2####################################
# GLOBAL TEMPERATURE EFFECT (all populations and days together)
#
# Rationale:
#   - The LRT showed a strong effect due to treatment
#   - Due to the limited sample size, we use the Wald test to identify which genes are up/down-regulated
#     due to temperature (all control samples vs all treatment samples) using a simpler additive model(),
#   - This ignores small interaction effects for simplicity, as justified by the LRT.

# load data
data1 <- read.csv("data/crubrum.gene.counts.matrix.csv", row.names = 1) #add data/
meta1 <- read.csv("data/crubrum.gene.coldata.csv", row.names = 1)

data1 <- round(data1)
(all(rownames(meta1) %in% colnames(data1)) || all(colnames(data1) %in% rownames(meta1)))
(all(colnames(data1) == rownames(meta1)))

meta1$day <- factor(meta1$day)

# Design: population + treatment + day (no interactions)
dds1 <- DESeqDataSetFromMatrix(
  countData = data1, colData = meta1,
  design = ~ population + treatment + day
)

# check the levels of each factor
dds1$population
dds1$treatment
dds1$day

# obtain normalized counts
dds1 <- estimateSizeFactors(dds1)
norm_cts_C_vs_T <- counts(dds1, normalized = TRUE)

# Keep only genes that have at least 10 reads in at least 3 samples
keep1 <- rowSums(norm_cts_C_vs_T >= 10) >= 3
cat("Genes before filter (Q2):", nrow(norm_cts_C_vs_T), "\n")
cat("Genes kept after normalized-count filter (Q2):", sum(keep1), "\n")

#recreate Deseq dataset
raw_counts1 <- counts(dds1, normalized = FALSE)
raw_counts1_filtered <- raw_counts1[keep1, , drop = FALSE]

dds1 <- DESeqDataSetFromMatrix(countData = raw_counts1_filtered,
                               colData = meta1,
                               design = ~ population + treatment + day)

# re-estimate size factors on the filtered dataset
dds1 <- estimateSizeFactors(dds1)

# Explicitly set the factor levels (control - reference)
dds1$treatment <- relevel(dds1$treatment, ref = "Control")

dds1 <- DESeq(dds1) # run DEseq2
resultsNames(dds1)

res1 <- results(dds1, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"), test = "Wald")
res1 <- res1[res1$baseMean > 10, ]

summary(res1)

results1 <- data.frame(res1)
results1 <- na.omit(results1)
write.csv(results1, file = "") # total results

# filter results
sig_results1 <- results1[results1$padj < 0.05, ]
sig_results1 <- subset(sig_results1, abs(sig_results1$log2FoldChange) > 0.3)
write.csv(sig_results1, file = "")
nrow(sig_results1)
# 2670

################################ Q3#######################################
# TEMPERATURE EFFECT WITHIN EACH POPULATION
#
# Rationale:
#   - The unified LRT confirmed a strong main effect of treatment and a negligible
#     treatment × population interaction, indicating that both populations were
#     overall affected by temperature stress,
#   - We now analyze each population separately to see how each responds to heat.
#   - By applying an additive design (~ treatment + day) within each population,
#     we estimate the average treatment response while accounting for temporal
#     variation, without overparameterizing the model due to small sampling size
#   - This allows us to assess local adaptation: whether shallow vs mesophotic colonies
#     activate different and specific transcriptional pathways under the same stress.

# load data

### shallow ####
cas_data <- read.csv("data/CAS_counts.csv", row.names = 1)
cas_meta <- read.csv("data/CAS_metadata.csv", row.names = 1)

cas_data <- round(cas_data)

(all(rownames(cas_meta) %in% colnames(cas_data)) || all(colnames(cas_data) %in% rownames(cas_meta)))
(all(colnames(cas_data) == rownames(cas_meta)))

cas_meta$day <- factor(cas_meta$day)

dds2 <- DESeqDataSetFromMatrix(
  countData = cas_data, colData = cas_meta,
  design = ~ treatment + day
)

# check the levels for each factor
dds2$treatment
dds2$day

# estimate size factors (so normalized counts are meaningful)
dds2 <- estimateSizeFactors(dds2)
norm_counts_cas <- counts(dds2, normalized = TRUE)

# Keep only genes that have at least 10 reads in at least 3 samples
keep2 <- rowSums(norm_counts_cas >= 10) >= 3
cat("Genes before filter (Q3):", nrow(norm_counts_cas), "\n")
cat("Genes kept after normalized-count filter (Q3):", sum(keep2), "\n")

#recreate Deseq dataset
raw_counts2 <- counts(dds2, normalized = FALSE)
raw_counts2_filtered <- raw_counts2[keep2, , drop = FALSE]

dds2 <- DESeqDataSetFromMatrix(countData = raw_counts2_filtered,
                               colData = cas_meta,
                               design = ~ treatment + day)

# re-estimate size factors on filtered dataset
dds2 <- estimateSizeFactors(dds2)

# Explicitly set the factor levels - set control as baseline again
dds2$treatment <- relevel(dds2$treatment, ref = "Control")

dds2 <- DESeq(dds2)
resultsNames(dds2)

res_cas <- results(dds2, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"), test = "Wald")
res_cas <- res_cas[res_cas$baseMean > 10, ]

summary(res_cas)
res_cas # This will be used for frontloading - keep it.

results_cas <- data.frame(res_cas)
results_cas <- na.omit(results_cas)
write.csv(results_cas, file = "") # total results

sig_results_cas <- results_cas[results_cas$padj < 0.05, ]
sig_results_cas <- subset(sig_results_cas, abs(sig_results_cas$log2FoldChange) > 0.3)
write.csv(sig_results_cas, file = "", row.names = T)
nrow(sig_results_cas)
# 692

### mesophotic ####

# load data
lop_data <- read.csv("data/LOP_counts.csv", row.names = 1)
lop_meta <- read.csv("data/LOP_metadata.csv", row.names = 1)

lop_data <- round(lop_data)

(all(rownames(lop_meta) %in% colnames(lop_data)) || all(colnames(lop_data) %in% rownames(lop_meta)))
(all(colnames(lop_data) == rownames(lop_meta)))

lop_meta$day <- factor(lop_meta$day)

dds3 <- DESeqDataSetFromMatrix(
  countData = lop_data, colData = lop_meta,
  design = ~ treatment + day
)

# check the levels for each factor
dds3$treatment
dds3$day

dds3 <- estimateSizeFactors(dds3)
norm_counts_lop <- counts(dds3, normalized = TRUE)

# Keep only genes that have at least 10 reads in at least 3 samples
keep3 <- rowSums(norm_counts_lop >= 10) >= 3
cat("Genes before filter (Q3):", nrow(norm_counts_lop), "\n")
cat("Genes kept after normalized-count filter (Q3):", sum(keep3), "\n")

#recreate Deseq dataset
raw_counts3 <- counts(dds3, normalized = FALSE)
raw_counts3_filtered <- raw_counts3[keep3, , drop = FALSE]

dds3 <- DESeqDataSetFromMatrix(countData = raw_counts3_filtered,
                               colData = lop_meta,
                               design = ~ treatment + day)

# re-estimate size factors
dds3 <- estimateSizeFactors(dds3)

# Explicitly set the factor levels
dds3$treatment <- relevel(dds3$treatment, ref = "Control")

dds3 <- DESeq(dds3)
resultsNames(dds3)

res_lop <- results(dds3, alpha = 0.05, contrast = c("treatment", "Treatment", "Control"), test = "Wald")
res_lop <- res_lop[res_lop$baseMean > 10, ]
summary(res_lop)
res_lop # this will be use for frontloading - keep it.

results_lop <- data.frame(res_lop)
results_lop <- na.omit(results_lop)
write.csv(results_lop, file = "") # total results

sig_results_lop <- results_lop[results_lop$padj < 0.05, ]
sig_results_lop <- subset(sig_results_lop, abs(sig_results_lop$log2FoldChange) > 0.3)
write.csv(sig_results_lop, file = "", row.names = T)
nrow(sig_results_lop)
# 1645

################################ Q4 #######################################
# TIME EFFECTS OF TEMPERATURE PER POPULATION
#
# Rationale:
#   - Temporal contrasts assess how the effect of thermal stress evolves across exposure times.
#   - We compared treatment vs control samples at each day (T0, T5, T10) within each population.
#   - This approach identifies genes whose expression changes due to heat stress at early (T0), intermediate (T5),  and late time points (T10),

## load data
data4 <- read.csv("data/crubrum.gene.counts.matrix.csv", row.names = 1)
meta4 <- read.csv("data/crubrum.gene.coldata.csv", row.names = 1)

data4 <- round(data4)

(all(rownames(meta4) %in% colnames(data4)) || all(colnames(data4) %in% rownames(meta4)))
(all(colnames(data4) == rownames(meta4)))

# defined factors
# Converting to factors.
meta4$treatment  <- factor(meta4$treatment, levels = c("Control", "Treatment"))
meta4$population <- factor(meta4$population, levels = c("Shallow", "Mesophotic"))
meta4$day        <- factor(meta4$day, levels = c("T0", "T5", "T10"))

# Create a single factor that encodes - dummy variable for population_treatment_day
meta4$factor <- factor(paste0(meta4$population, "_", meta4$treatment, "_", meta4$day))

dds4 <- DESeqDataSetFromMatrix(countData = data4, colData = meta4, design = ~ factor)
dds4 <- estimateSizeFactors(dds4)
norm_counts4 <- counts(dds4, normalized = TRUE)

# Keep genes with >= 10 normalized reads in >= 3 samples
keep_genes4 <- rowSums(norm_counts4 >= 10) >= 3
cat("Q4 - genes before normalized-count filter:", nrow(norm_counts4), "\n")
cat("Q4 - genes kept after normalized-count filter:", sum(keep_genes4), "\n")

# recreate Deseq dataset
raw_counts4 <- counts(dds4, normalized = FALSE)
raw_counts4_filtered <- raw_counts4[keep_genes4, , drop = FALSE]

dds4 <- DESeqDataSetFromMatrix(
  countData = raw_counts4_filtered, colData = meta4,
  design = ~factor
)

# re-estimate size factors on the filtered dataset
dds4 <- estimateSizeFactors(dds4)

# check the levels for each factor
dds4$factor

# Set one reference
dds4$factor <- relevel(dds4$factor, ref = "Shallow_Control_T0")

# Run DESeq
dds4 <- DESeq(dds4)
resultsNames(dds4)

#### shallow ####
# contrasts:
#### T0 (treatment T0 vs Control T0) for Shallow population
res_CAS_day0 <- results(dds4, alpha = 0.05, contrast = c("factor", "Shallow_Treatment_T0", "Shallow_Control_T0")) # So the genes are up/ down regulated in treat day0 vs ctrl day0
res_CAS_day0 <- res_CAS_day0[res_CAS_day0$baseMean > 10, ]
summary(res_CAS_day0)
res_CAS_day0

res_CAS_day0 <- data.frame(res_CAS_day0)
res_CAS_day0 <- na.omit(res_CAS_day0)
write.csv(res_CAS_day0, file = "")

sig_results_day0_CAS <- res_CAS_day0[res_CAS_day0$padj < 0.05, ]
sig_results_day0_CAS <- subset(sig_results_day0_CAS, abs(sig_results_day0_CAS$log2FoldChange) > 0.3)
write.csv(sig_results_day0_CAS, file = "", row.names = T)
nrow(sig_results_day0_CAS)
# 1

### T5 (treatment T5 vs Control T5) for shallow population
res_CAS_day5 <- results(dds4, alpha = 0.05, contrast = c("factor", "Shallow_Treatment_T5", "Shallow_Control_T5")) # So the genes are up/ down regulated in treat day0 vs ctrl day0
res_CAS_day5 <- res_CAS_day5[res_CAS_day5$baseMean > 10, ]
summary(res_CAS_day5)
res_CAS_day5

res_CAS_day5 <- data.frame(res_CAS_day5)
res_CAS_day5 <- na.omit(res_CAS_day5)
write.csv(res_CAS_day5, file = "")

sig_results_day5_CAS <- res_CAS_day5[res_CAS_day5$padj < 0.05, ]
sig_results_day5_CAS <- subset(sig_results_day5_CAS, abs(sig_results_day5_CAS$log2FoldChange) > 0.3)
write.csv(sig_results_day5_CAS, file = "", row.names = T)
nrow(sig_results_day5_CAS)
# 697

### T10 (treatment T10 vs Control T10)
res_CAS_day10 <- results(dds4, alpha = 0.05, contrast = c("factor", "Shallow_Treatment_T10", "Shallow_Control_T10")) # So the genes are up/ down regulated in treat day0 vs ctrl day0
res_CAS_day10 <- res_CAS_day10[res_CAS_day10$baseMean > 10, ]
summary(res_CAS_day10)
res_CAS_day10

res_CAS_day10 <- data.frame(res_CAS_day10)
res_CAS_day10 <- na.omit(res_CAS_day10)
write.csv(res_CAS_day10, file = "")

sig_results_day10_CAS <- res_CAS_day10[res_CAS_day10$padj < 0.05, ]
sig_results_day10_CAS <- subset(sig_results_day10_CAS, abs(sig_results_day10_CAS$log2FoldChange) > 0.3)
write.csv(sig_results_day10_CAS, file = "", row.names = T)
nrow(sig_results_day10_CAS)
# 514

#### mesophotic ####

#### T0 (treatment T0 vs Control T0)

res_LOP_day0 <- results(dds4, alpha = 0.05, contrast = c("factor", "Mesophotic_Treatment_T0", "Mesophotic_Control_T0")) # So the genes are up/ down regulated in treat day0 vs ctrl day0
res_LOP_day0 <- res_LOP_day0[res_LOP_day0$baseMean > 10, ]
summary(res_LOP_day0)
res_LOP_day0

res_LOP_day0 <- data.frame(res_LOP_day0)
res_LOP_day0 <- na.omit(res_LOP_day0)
write.csv(res_LOP_day0, file = "")

sig_results_day0_LOP <- res_LOP_day0[res_LOP_day0$padj < 0.05, ]
sig_results_day0_LOP <- subset(sig_results_day0_LOP, abs(sig_results_day0_LOP$log2FoldChange) > 0.3)
write.csv(sig_results_day0_LOP, file = "", row.names = T)
nrow(sig_results_day0_LOP)
# 6

### T5 (treatment T5 vs Control T5)
res_LOP_day5 <- results(dds4, alpha = 0.05, contrast = c("factor", "Mesophotic_Treatment_T5", "Mesophotic_Control_T5")) # So the genes are up/ down regulated in treat day0 vs ctrl day0
res_LOP_day5 <- res_LOP_day5[res_LOP_day5$baseMean > 10, ]
summary(res_LOP_day5)
res_LOP_day5

res_LOP_day5 <- data.frame(res_LOP_day5)
res_LOP_day5 <- na.omit(res_LOP_day5)
write.csv(res_LOP_day5, file = "")

sig_results_day5_LOP <- res_LOP_day5[res_LOP_day5$padj < 0.05, ]
sig_results_day5_LOP <- subset(sig_results_day5_LOP, abs(sig_results_day5_LOP$log2FoldChange) > 0.3)
write.csv(sig_results_day5_LOP, file = "", row.names = T)
nrow(sig_results_day5_LOP)
# 2470

### T10 (treatment T10 vs control T10)
res_LOP_day10 <- results(dds4, alpha = 0.05, contrast = c("factor", "Mesophotic_Treatment_T10", "Mesophotic_Control_T10")) # So the genes are up/ down regulated in treat day0 vs ctrl day0
res_LOP_day10 <- res_LOP_day10[res_LOP_day10$baseMean > 10, ]
summary(res_LOP_day10)
res_LOP_day10

res_LOP_day10 <- data.frame(res_LOP_day10)
res_LOP_day10 <- na.omit(res_LOP_day10)
write.csv(res_LOP_day10, file = "")

sig_results_day10_LOP <- res_LOP_day10[res_LOP_day10$padj < 0.05, ]
sig_results_day10_LOP <- subset(sig_results_day10_LOP, abs(sig_results_day10_LOP$log2FoldChange) > 0.3)
write.csv(sig_results_day10_LOP, file = "", row.names = T)
nrow(sig_results_day10_LOP)
# 765

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

# Step 1 - use res_cas and res_lop objects results already generated, otherwise:

# Load results from mesophotic under stress if you saved them
#res_lop <- read.csv(file = "res_lop.csv", row.names = NULL) # Unremark and complete this if you load the file
#rownames(res_lop) <- res_lop$SeqName
#res_lop$SeqName <- NULL

# load results from shallow under stress if you saved them
#res_cas <- read.csv(file = "res_cas.csv", row.names = NULL) # Unremark and complete this if you load the file
#rownames(res_cas) <- res_cas$SeqName
#res_cas$SeqName <- NULL

# select matching genes
common_genes <- intersect(rownames(data.frame(res_lop)), rownames(data.frame(res_cas)))

# step 2: generate conjunct table
df <- data.frame(
  SeqName = common_genes,

  # Mesophotic (LOP) Columns
  LOP_baseMean = res_lop[common_genes, "baseMean"],
  LOP_log2FoldChange = res_lop[common_genes, "log2FoldChange"],
  LOP_padj = res_lop[common_genes, "padj"],
  LOP_sig = ifelse(res_lop[common_genes, "padj"] < 0.05, "sig", ""),

  # Shallow (CAS) Columns
  CAS_baseMean = res_cas[common_genes, "baseMean"],
  CAS_log2FoldChange = res_cas[common_genes, "log2FoldChange"],
  CAS_padj = res_cas[common_genes, "padj"],
  CAS_sig = ifelse(res_cas[common_genes, "padj"] < 0.05, "sig", "")
)

# step 3 - add baseline comparison — Load shallow vs mesophotic baseline (Control, T0). Use data and meta objects already loaded, otherwise do:
# load data
data <- read.csv("data/counts_cas_lop_ctrl_day0.csv", row.names = 1) 
meta <- read.csv("data/coldata_cas_lop_ctrl_day0.csv", row.names = 1) 

# recheck data
data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

meta$day <- factor(meta$day)
meta$population <- factor(meta$population, levels = c("CAS", "LOP"))

dds_base <- DESeqDataSetFromMatrix(
  countData = data, colData = meta,
  design = ~population
)

# extract normalized counts for baseline expression
dds_base <- estimateSizeFactors(dds_base)
norm_base <- counts(dds_base, normalized = TRUE)

# Apply filtering
keep_base <- rowSums(norm_base >= 10) >= 3
norm_base <- norm_base[keep_base, , drop = FALSE]

### mean baseline expression in shallow and mesophotic (Control, T0)
cas_samples <- rownames(meta)[meta$population == "CAS"]
lop_samples <- rownames(meta)[meta$population == "LOP"]

mean_cas <- rowMeans(norm_base[, cas_samples])
mean_lop <- rowMeans(norm_base[, lop_samples])

# mean baseline expression
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

df$Control_log2FoldChange <- log2((df$Control_baseMean_CAS + 1) /
  (df$Control_baseMean_LOP + 1))
# add log fold change difference
df$Control_sig <- ifelse(abs(df$Control_log2FoldChange) > 0.3, "sig", "")

# Step 4: define our gene classifications:

## Conditions
sig_cutoff <- 0.05
lfc_cutoff <- 0.3

lop_sig <- df$LOP_padj < sig_cutoff & abs(df$LOP_log2FoldChange) > lfc_cutoff
cas_sig <- df$CAS_padj < sig_cutoff & abs(df$CAS_log2FoldChange) > lfc_cutoff

## Frontloaded genes are:
frontloaded <- (df$Control_baseMean_CAS > df$Control_baseMean_LOP) & lop_sig

## Reduced reaction genes are (shallow reacts less):
reduced_reaction <- lop_sig & cas_sig &
  (df$LOP_log2FoldChange > df$CAS_log2FoldChange)

## Greater logfold change (shallow reacts more):
greater_reaction <- lop_sig & cas_sig &
  (df$CAS_log2FoldChange > df$LOP_log2FoldChange)

## label left genes as NAs
df$Classification <- "NA"
df$Classification[greater_reaction] <- "greater fold change"
df$Classification[reduced_reaction] <- "reduced reaction"
df$Classification[frontloaded] <- "frontloaded"

# Save full table:
write.csv(df, "",
  row.names = FALSE
) # filter table with upregulated mesophotic DEGs relevant for the study - add annotation and run enrichment.
