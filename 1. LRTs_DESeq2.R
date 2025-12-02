################ Code accompanying: ########################################
#
#   Ramirez-Calero, S. et al. 2025. Does local adaptation influence thermal responses
#   in red coral populations across depth gradients? Transcriptomic insights for effective conservation.
#
# Script written by: Sandra Ramirez
#
# This script must be run within a Rproject containg the data folder and the renv.lock file.
# This script was run to obtain significant differential expressed genes due to a
# thermal treatment stress (24°C) vs a control (18°C) for two contrasting populations
# of octocoral Corallium rubrum across three different time points (T0, T5 and T10).
# Data sets used under this script are located in the data folder.
#
# Biological design:
#   2 Populations: CAS (shallow) and LOP (mesophotic)
#   2 Treatments: Control (18°C) vs Treatment (25°C)
#   3 Time points: T0, T5 and T10
#
#####################################################################

# Before running the script for the first time on your machine:
renv::restore() # To retrieve the exact library used for this project.

# load libraries
library(digest)
library(XML)
library(RSQLite)
library(DESeq2)
library(pcaExplorer)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("")

# load data
data <- read.csv("data/crubrum.gene.counts.matrix.csv", row.names = 1)
meta <- read.csv("data/crubrum.gene.coldata.csv", row.names = 1)

data <- round(data)
(all(rownames(meta) %in% colnames(data)) || all(colnames(data) %in% rownames(meta)))
(all(colnames(data) == rownames(meta)))

# Converting to factors.
meta$treatment  <- factor(meta$treatment, levels = c("Control", "Treatment"))
meta$population <- factor(meta$population, levels = c("Shallow", "Mesophotic"))
meta$day        <- factor(meta$day, levels = c("T0", "T5", "T10"))

# filter low expressed genes
cat("Raw genes", nrow(data), "genes before filtering")
keep <- rowSums(data >= 10) >= 3 # keep genes with at least 10 reads in >= 3 samples
data <- data[keep, ]
cat("Kept genes:", nrow(data), "genes after filtering")

# Building one unique full model to use for subsequent analysis
dds_full <- DESeqDataSetFromMatrix(countData = data,
                                   colData = meta,
                                   design = ~ treatment + day + population + treatment:day + treatment:population)

# exploratory PCA for all samples
vsd <- vst(dds_full, blind = FALSE)

rv <- rowVars(assay(vsd), useNames = TRUE)
topVarGenes <- head(order(rv, decreasing = TRUE), 5000)

pcaData <- plotPCA(vsd[topVarGenes, ], intgroup = c("population", "treatment", "day"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pop_shapes <- c("Shallow" = 21, "Mesophotic" = 24)
treat_colors <- c("Control" = "#6E8B3D", "Treatment" = "#FF7F50")
day_fill <- c("T0" = "#97FFFF", "T5" = "#528B8B", "T10" = "#2F4F4F")

ggplot(pcaData, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = treatment, group = interaction(population, treatment)),
    type = "norm", linetype = 1, linewidth = 1, alpha = 1
  ) +
  geom_point(aes(shape = population, fill = day, color = treatment),
    size = 3.8, stroke = 1.2, alpha = 0.95
  ) +
  scale_shape_manual(values = pop_shapes) +
  scale_color_manual(values = treat_colors) +
  scale_fill_manual(values = day_fill) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), legend.position = "right") +
  guides(
    shape = guide_legend(title = "Population"),
    color = guide_legend(title = "Treatment"),
    fill = guide_legend(title = "Day")
  )

# Save plot
ggsave("~")

###############################################################################
# 1) Test for interactions with LRT (global test). Rationale:
#    - If treatment effects differ across days or populations, main-effect tests
#      are not sufficient. We need to examine effects within each level.
#    - LRT compares the full model (with interaction) vs reduced model (without).
###############################################################################

# treatment x time interaction

# reduced model - remove the interaction term
dds_reduced <- DESeq(dds_full,
                     test = "LRT",
                     reduced = ~ treatment + population + day + treatment:population)

res_interaction_effect <- results(dds_reduced)
res_interaction_effect
summary(res_interaction_effect)
results_interaction_effect <- data.frame(res_interaction_effect)
results_interaction_effect <- na.omit(results_interaction_effect)

res_sig_interaction_effect <- results_interaction_effect[results_interaction_effect$padj < 0.05, ]
### 999 genes DE due to interaction effect

# treatment x population interaction

# reduced model - remove the interaction term
dds_reduced <- DESeq(dds_full,
                     test = "LRT",
                     reduced = ~ treatment + population + day + treatment:day)

res_pop_int <- results(dds_reduced)
summary(res_pop_int)
results_pop_int <- data.frame(res_pop_int)
results_pop_int <- na.omit(results_pop_int)

res_sig_pop_int <- results_pop_int[results_pop_int$padj < 0.05, ]
### 1 genes DE due to pop_int effect

# Interpretation:
# - treatment:day produced 616 genes indicating a moderate differential response across time points
# - treatment:population produced none, so there is little evidence that treatment effect *depends* on population
# => it is reasonable to proceed by modelling main effects without interactions for global summaries

#### treatment effect #####

###############################################################################
# 2) Main-factor global tests (LRTs): treatment, population, day
#    Rationale: LRT here deals with "does adding this factor to the model improve fit?"
#    This is a global test across all levels (for multiple factors) or overall effect.
###############################################################################

# reduced model - remove the Treatment term and its interactions
dds_reduced_treat <- DESeq(dds_full, test = "LRT", reduced = ~ day + population)

res_treatment <- results(dds_reduced_treat)
res_treatment
summary(res_treatment)
results_treatment <- data.frame(res_treatment)
results_treatment <- na.omit(results_treatment)

res_sig_treatment <- results_treatment[results_treatment$padj < 0.05, ]
### 3018 genes DE due to treatment effect

# Summary:
# The unified LRTs confirmed a strong treatment effect and negligible treatment × population interaction.
# We therefore proceed with Wald tests to quantify:
# 1) baseline population differences under control conditions (frontloading),
# 2) overall treatment effects in C. rubrum,
# 3) population-specific stress responses,
# 4) treatment effects across time points.

# Proceed with script Wald_tests_DESEq2.R
