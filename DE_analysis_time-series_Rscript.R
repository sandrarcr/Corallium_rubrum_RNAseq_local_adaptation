################################### Time-series ####################################

#set working directory
setwd("C:/Users/Sandra/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2") #mypc
setwd("C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/2. DESeq2") #Quim's PC

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

#To fix this, if happens, you can convert these variables to factors using the factor() 
#function:

coldata$day <- factor(coldata$day)
coldata$population <- factor(coldata$population)
coldata$Individual <- factor(coldata$Individual)

dds_time <- DESeqDataSetFromMatrix(countData = cts1, colData = coldata1, design = ~day + treatment + day:treatment + treatment:population) #full model

#day: This variable represents the time points or days at which the experiment was conducted. 
#It suggests that you have measured gene expression at different time points and want to account 
#for any time-dependent effects on gene expression.
#treatment: This variable represents the treatment conditions or experimental groups to which your samples belong. 
#It suggests that you want to compare the gene expression between different treatment groups.
#day:treatment: This component represents the interaction between the day and treatment variables. 
#It allows you to assess whether the effect of the treatment on gene expression varies depending on the time point. 
#In other words, it allows you to identify genes that show different responses to treatment at different time points.

#factors
dds_time$day <- factor(dds_time$day, levels = c("0","1","2")) #we set these factors as we do want to see the diff exp genes across the time of treatment
#taking into account all factors playing
dds_time$population <- factor(dds_time$population, levels = c("CAS", "LOP"))

# Remove genes which have at least 10 reads in all samples
keep <- rowSums(counts(dds_time)) >= 10
dds_time <- dds_time[keep,]

######## Run Differential expression Analysis

dds_LRT <- DESeq(dds_time, test="LRT", reduced = ~ day + treatment + day:treatment) #reduced --The LRT examines two 
#models for the counts, a full model with a certain number of terms and a reduced model,
#in which some of the terms of the full model are removed. The test determines if the 
#increased likelihood of the data using the extra terms in the full model is more than 
#expected if those extra terms are truly zero. The LRT is therefore useful for testing 
#multiple terms at once, for example testing 3 or more levels of a factor at once, or
#all interactions between two variables.
# you would need to decide which terms or interactions you want to test for significance 
#by removing them from the full model. In this case, let's say you want to test the 
#significance of the interaction term day:treatment while keeping all other terms.
#so you do: ~day:treatment + treatment:population

save(dds_LRT, file = "dds_lrt.RData")
#load(file = "C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/DESeq2/15. dds~pop+day+pop-day(timeseries)/dds_lrt_time.RData")

#Sandra- these are the new comparisons that appear when I include your suggestion in the design:

resultsNames(dds_LRT)
#[1] "Intercept"                        "day_1_vs_0"                       "day_2_vs_0"                       "treatment_Treatment_vs_Control"  
#[5] "day1.treatmentTreatment"          "day2.treatmentTreatment"          "treatmentControl.populationLOP"   "treatmentTreatment.populationLOP"

#"Intercept": This represents the baseline or reference level of gene expression when all variables in the model are at 
#their reference levels.
#"day_1_vs_0": This contrast compares the gene expression at day 1 to day 0, while keeping the treatment condition constant.
#"day_2_vs_0": This contrast compares the gene expression at day 2 to day 0, while keeping the treatment condition constant.
#"treatment_Treatment_vs_Control": This contrast compares the treatment condition (Treatment) to the control condition 
#(Control), while keeping the day variable constant.
#"day1.treatmentTreatment": This represents the interaction effect between day 1 and the treatment condition. It examines 
#whether the effect of treatment on gene expression is different at day 1 compared to the reference level (day 0).
#"day2.treatmentTreatment": This represents the interaction effect between day 2 and the treatment condition. It examines 
#whether the effect of treatment on gene expression is different at day 2 compared to the reference level (day 0).
#"treatmentControl.populationLOP": This represents the interaction between the treatment group and the population level "LOP" in the "Control" condition.
#"treatmentTreatment.populationLOP": This represents the interaction between the treatment group and the population level "LOP" in the "Treatment" condition.

#Tutorial: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#Two factors with two and three levels and interactions
#Day: 0, 1 and 2
#condition: control and treatment

# get the model matrix
mod_mat <- model.matrix(design(dds_LRT), colData(dds_LRT))

# Define coefficient vectors for each condition
day1_control <- colMeans(mod_mat[dds_LRT$day == "1" & dds_LRT$treatment == "Control", ])
day1_treat <- colMeans(mod_mat[dds_LRT$day == "1" & dds_LRT$treatment == "Treatment", ])
day0_control <- colMeans(mod_mat[dds_LRT$day == "0" & dds_LRT$treatment == "Control", ])
day0_treat <- colMeans(mod_mat[dds_LRT$day == "0" & dds_LRT$treatment == "Treatment", ])
day2_control <- colMeans(mod_mat[dds_LRT$day == "2" & dds_LRT$treatment == "Control", ])
day2_treat <- colMeans(mod_mat[dds_LRT$day == "2" & dds_LRT$treatment == "Treatment", ])
LOP_control <- colMeans(mod_mat[dds_LRT$population == "LOP" & dds_LRT$treatment == "Control", ])
LOP_treat <- colMeans(mod_mat[dds_LRT$population == "LOP" & dds_LRT$treatment == "Treatment", ])
LOP_0 <- colMeans(mod_mat[dds_LRT$population == "LOP" & dds_LRT$day == "0", ])
LOP_1 <- colMeans(mod_mat[dds_LRT$population == "LOP" & dds_LRT$day == "1", ])
LOP_2 <- colMeans(mod_mat[dds_LRT$population == "LOP" & dds_LRT$day == "2", ])
CAS_control <- colMeans(mod_mat[dds_LRT$population == "CAS" & dds_LRT$treatment == "Control"])
CAS_treat <- colMeans(mod_mat[dds_LRT$population == "CAS" & dds_LRT$treatment == "Treatment", ])
CAS_0 <- colMeans(mod_mat[dds_LRT$population == "CAS" & dds_LRT$day == "0", ])
CAS_1 <- colMeans(mod_mat[dds_LRT$population == "CAS" & dds_LRT$day == "1", ])
CAS_2 <- colMeans(mod_mat[dds_LRT$population == "CAS" & dds_LRT$day == "2", ])

#obtain results from dds object

#Sandra - here i'm a bit confused about what i'm obtaining.... I'm using contrasts for each as DESeq2 tutorial suggests, but i'm still unsure... :'(

###Overall results
res <-results(dds_LRT, alpha = 0.05, test = "Wald") #general results including all ind and counts. differential expression results for all comparisons or contrasts included in your model. This includes the fold changes, p-values, and adjusted p-values (if applicable) for each gene in the dataset.
res
#filter baseMean genes > 10
res <- res[res$baseMean > 10, ]
summary(res) #day_1_vs_0": This represents the comparison between day 1 and day 0, 

###day 1 vs 0 in the control
res1 <- results(dds_LRT, alpha = 0.05, contrast = day1_control - day0_control, test = "Wald")
# or equivalently
res1 <- results(dds_LRT, alpha = 0.05, contrast = list("day_1_vs_0"), test = "Wald")
res1
#filter baseMean genes > 10
res1 <- res1[res1$baseMean > 10, ]
summary(res1) #day_1_vs_0": This represents the comparison between day 1 and day 0, 
#specifically testing the difference in gene expression between these two time points.

###day 1 vs 0 (in the Treatment):
res2 <- results(dds_LRT, alpha = 0.05, contrast = day1_treat - day0_treat, test = "Wald")
# or equivalently
res2 <- results(dds_LRT, alpha = 0.05, contrast = list(c("day_1_vs_0",
                                       "day1.treatmentTreatment")), test = "Wald")
res2
#filter baseMean genes > 10
res2 <- res2[res2$baseMean > 10, ]
summary(res2) 

###day 2 vs 0 (in control)
res3 <- results(dds_LRT, alpha = 0.05, contrast = day2_control - day0_control, test = "Wald")
# or equivalently
res3 <- results(dds_LRT, alpha = 0.05,  contrast = list("day_2_vs_0"), test = "Wald")
res3
#filter baseMean genes > 10
res3 <- res3[res3$baseMean > 10, ]
summary(res3) #"day_2_vs_0": This represents the comparison between day 2 and day 0, 
#testing the difference in gene expression between these two time points.

###day 2 vs 0 (in the Treatment):
res4 <- results(dds_LRT, alpha = 0.05, contrast = day2_treat - day0_treat, test = "Wald")
# or equivalently
res4 <- results(dds_LRT, alpha = 0.05, contrast = list(c("day_2_vs_0",
                                           "day2.treatmentTreatment")), test = "Wald")
res4
#filter baseMean genes > 10
res4 <- res4[res4$baseMean > 10, ]
summary(res4) 

###Treatment vs Control (for day0):
res5 <- results(dds_LRT, alpha = 0.05, contrast = day0_treat - day0_control, test = "Wald")
# or equivalently
res5 <- results(dds_LRT, contrast = list(c("treatment_Treatment_vs_Control")), test = "Wald")
res5
#filter baseMean genes > 10
res5 <- res5[res5$baseMean > 10, ]
summary(res5) #treatment_Treatment_vs_Control": This represents the comparison between the 
#treatment group and the control group, specifically examining the difference in gene expression 
#between these two conditions.

###Treatment vs Control (for day 1):
res6 <- results(dds_LRT, alpha = 0.05, contrast = day1_treat - day1_control, test = "Wald")
# or equivalently
res6 <- results(dds_LRT, alpha = 0.05, contrast = list(c("treatment_Treatment_vs_Control", 
                                       "day1.treatmentTreatment")), test = "Wald")
res6
#filter baseMean genes > 10
res6 <- res6[res6$baseMean > 10, ]
summary(res6) 

###Treatment vs Control (for day 2):
res7 <- results(dds_LRT, alpha = 0.05, contrast = day2_treat - day2_control, test = "Wald")
# or equivalently
res7 <- results(dds_LRT, alpha = 0.05, contrast = list(c("treatment_Treatment_vs_Control", 
                                       "day2.treatmentTreatment")), test = "Wald")
res7
#filter baseMean genes > 10
res7 <- res7[res7$baseMean > 10, ]
summary(res7) 

###Interaction between day and condition (i.e. do day0, 1 and 2 respond differently to the treatment?):

#day 1 and 0 respond differently to treatment?
res8 <- results(dds_LRT, alpha = 0.05, contrast = (day1_treat - day1_control) - (day0_treat - day0_control), test = "Wald")
# or equivalently
res8 <- results(dds_LRT, alpha = 0.05, contrast = list("day1.treatmentTreatment"), test = "Wald")
res8
#filter baseMean genes > 10
res8 <- res8[res8$baseMean > 10, ]
summary(res8) #"day1.treatmentTreatment": This represents the interaction term between "day" and "treatment" for day 1, 
#specifically testing the difference in gene expression between the treatment group and the control group at day 1.

#day 2 and 0 respond differently to treatment?
res9 <- results(dds_LRT, alpha = 0.05, contrast = (day2_treat - day2_control) - (day0_treat - day0_control), test = "Wald")
# or equivalently
res9 <- results(dds_LRT, alpha = 0.05, contrast = list("day2.treatmentTreatment"), test = "Wald")
res9
#filter baseMean genes > 10
res9 <- res9[res9$baseMean > 10, ]
summary(res9) #"day2.treatmentTreatment": This represents the interaction term between "day" and "treatment" for day 2, 
#testing the difference in gene expression between the treatment group and the control group at day 2.

### LOP
res10 <- results(dds_LRT, alpha = 0.05, contrast = LOP_treat - LOP_control, test = "Wald")
res10
#filter baseMean genes > 10
res10 <- res10[res10$baseMean > 10, ]
summary(res10)

## LOP_treatment_day0 

res11 <- results(dds_LRT, alpha = 0.05, contrast = LOP_0 - LOP_treat, test = "Wald")
res11
#filer baseMean genes > 10
res11 <- res11[res11$baseMean > 10, ]
summary(res11)

### LOP_treatment_day1 

res12 <- results(dds_LRT, alpha = 0.05, contrast = LOP_1 - LOP_treat, test = "Wald")
res12
#filer baseMean genes > 10
res12 <- res12[res12$baseMean > 10, ]
summary(res12)

### LOP_treatment_day2 
  
res13 <- results(dds_LRT, alpha = 0.05, contrast = LOP_1 - LOP_treat, test = "Wald")
res13
#filer baseMean genes > 10
res13 <- res13[res13$baseMean > 10, ]
summary(res13)
  
### CAS
res14 <- results(dds_LRT, alpha = 0.05, contrast = CAS_treat - CAS_control, test = "Wald")
res14
#filter baseMean genes > 10
res14 <- res14[res14$baseMean > 10, ]
summary(res14)

### CAS_treatment_day0 

res15 <- results(dds_LRT, alpha = 0.05, contrast = CAS_0 - CAS_treat, test = "Wald")
res15
#filer baseMean genes > 10
res15 <- res15[res15$baseMean > 10, ]
summary(res15)

### CAS_treatment_day1 

res16 <- results(dds_LRT, alpha = 0.05, contrast = CAS_1 - CAS_treat, test = "Wald")
res16
#filer baseMean genes > 10
res16 <- res16[res16$baseMean > 10, ]
summary(res16)
  
### CAS_treatment_day2 

res17 <- results(dds_LRT, alpha = 0.05, contrast = CAS_2 - CAS_treat, test = "Wald")
res17
#filer baseMean genes > 10
res17 <- res17[res17$baseMean > 10, ]
summary(res17)

#whether CAS and LOP respond differently to Control-treatment

res18 <- results(dds_LRT, alpha = 0.05, contrast = (LOP_treat - LOP_control) - (CAS_treat - CAS_control), test = "Wald")
res18
#filer baseMean genes > 10
res18 <- res18[res18$baseMean > 10, ]
summary(res18)

#Sandra - In essence, I would like to obtain diff exp genes between cont-treat for LOP and CAS independently... I think I did it. Can you confirm?
#Sandra - I would also like to obtain diff exp genes for the treatment in LOP at day 0, 1 and 2... please let me know i did it correctly - same for CAS
#Sandra - If this is correct... then the rest should also be ok :) 


#save results:
save(res, file = "res.RData")  
save(res1, file = "res1vs0_control.RData")
save(res2, file = "res1cs0_treat.RData")
save(res3, file = "res2vs0_control.RData")
save(res4, file = "res2vs0_treat.RData")
save(res5, file = "res_treatvscontrol_0.RData")
save(res6, file = "res_treatvscontrol_1.RData")
save(res7, file = "res_treatvscontrol_2.RData")
save(res8, file = "res_int_1vs0_treat.RData")
save(res9, file = "res_int_2vs0_treat.RData")
save(res10, file = "res_LOP_treatvscontrol.RData")
save(res11, file = "res_LOP_0_treat.RData")
save(res12, file = "res_LOP_1_treat.RData")
save(res13, file = "res_LOP_2_treat.RData")
save(res14, file = "res_CAS_treatvscontrol.RData")
save(res15, file = "res_CAS_0_treat.RData")
save(res16, file = "res_CAS_1_treat.RData")
save(res17, file = "res_CAS_2_treat.RData")
save(res18, file = "res_int_LOPvsCAS_treat.RData")

#load(file = "C:/Users/Joaquim Garrabou/OneDrive - Universitat de Barcelona/Doctorado/ICM/Tesis/Chap2/RNA-Seq/DESeq2/15. dds~pop+day+pop-day(timeseries)/1) treatment + day + day-treatment/resX.RData")

plotDispEsts(dds_LRT)

#order res based on adjusted smallest p-values
resOrdered <- res[order(res$padj),]
resOrdered1 <- res1[order(res1$padj),]
resOrdered2 <- res2[order(res2$padj),]
resOrdered3 <- res2[order(res3$padj),]
resOrdered4 <- res2[order(res4$padj),]
resOrdered5 <- res2[order(res5$padj),]
resOrdered6 <- res2[order(res6$padj),]
resOrdered7 <- res2[order(res7$padj),]
resOrdered8 <- res2[order(res8$padj),]
resOrdered9 <- res2[order(res9$padj),]
resOrdered10 <- res1[order(res10$padj),]
resOrdered11 <- res1[order(res11$padj),]
resOrdered12 <- res2[order(res12$padj),]
resOrdered13 <- res2[order(res13$padj),]
resOrdered14 <- res2[order(res14$padj),]
resOrdered15 <- res2[order(res15$padj),]
resOrdered16 <- res2[order(res16$padj),]
resOrdered17 <- res2[order(res17$padj),]
resOrdered18 <- res2[order(res18$padj),]

summary(resX)

######
plotMA(res16, ylim=c(-12,12)) #shows the log2 fold changes attributable to a given variable over the mean of normalized 
#counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less 
#than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

plotCounts(dds_LRT, gene = which.min(res14$padj), intgroup = "treatment")
#plotCounts(dds_LRT, gene = which.min(res2$padj), intgroup = "treatment")

######## extract data from res objects
re <- data.frame(resOrdered)
re <- na.omit(re)
write.csv(re,file="deseq2_general_day-C-T_crubrum.csv")

#extract data from res object
re1 <- data.frame(resOrdered1)
re1 <- na.omit(re1)
write.csv(re1,file="deseq2_day1vs0_C_crubrum.csv")

re2 <- data.frame(resOrdered2)
re2 <- na.omit(re2)
write.csv(re2,file="deseq2_day1vs0_T_crubrum.csv")

#extract data from res object
re3 <- data.frame(resOrdered3)
re3 <- na.omit(re3)
write.csv(re3,file="deseq2_day2vs0_C_crubrum.csv")

re4 <- data.frame(resOrdered4)
re4 <- na.omit(re4)
write.csv(re4,file="deseq2_day2vs0_T_crubrum.csv")

#extract data from res object
re5 <- data.frame(resOrdered5)
re5 <- na.omit(re5)
write.csv(re5,file="deseq2_TvsC_day0_crubrum.csv")

re6 <- data.frame(resOrdered6)
re6 <- na.omit(re6)
write.csv(re6,file="deseq2_TvsC_day1_crubrum.csv")

#extract data from res object
re7 <- data.frame(resOrdered7)
re7 <- na.omit(re7)
write.csv(re7,file="deseq2_TvsC_day2_crubrum.csv")

re8 <- data.frame(resOrdered8)
re8 <- na.omit(re8)
write.csv(re8,file="deseq2_int_1vs0_T_crubrum.csv")

#extract data from res object
re9 <- data.frame(resOrdered9)
re9 <- na.omit(re9)
write.csv(re9,file="deseq2_int_2vs0_T_crubrum.csv")

#extract data from res object
re10 <- data.frame(resOrdered10)
re10 <- na.omit(re10)
write.csv(re10,file="deseq2_LOP_TvsC_C_crubrum.csv")

#extract data from res object
re11 <- data.frame(resOrdered11)
re11 <- na.omit(re11)
write.csv(re11,file="deseq2_LOP_0_T_crubrum.csv")

re12 <- data.frame(resOrdered12)
re12 <- na.omit(re12)
write.csv(re12,file="deseq2_LOP_1_T_crubrum.csv")

#extract data from res object
re13 <- data.frame(resOrdered13)
re13 <- na.omit(re13)
write.csv(re13,file="deseq2_LOP_2_T_crubrum.csv")

re14 <- data.frame(resOrdered14)
re14 <- na.omit(re14)
write.csv(re14,file="deseq2_CAS_TvsC_crubrum.csv")

#extract data from res object
re15 <- data.frame(resOrdered15)
re15 <- na.omit(re15)
write.csv(re15,file="deseq2_CAS_0_T_crubrum.csv")

re16 <- data.frame(resOrdered16)
re16 <- na.omit(re16)
write.csv(re16,file="deseq2_CAS_1_T_crubrum.csv")

#extract data from res object
re17 <- data.frame(resOrdered17)
re17 <- na.omit(re17)
write.csv(re17,file="deseq2_CAS_2_T_crubrum.csv")

re18 <- data.frame(resOrdered18)
re18 <- na.omit(re18)
write.csv(re18,file="deseq2_int_LOPvsCAS_T_crubrum.csv")

#### extract data which pass an adjusted p value threshold (0.05)

sig=subset(re,re$padj<0.05)
sig1=subset(sig,abs(sig$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig1, file = "sig2_deseq2_general_day-C-T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig1,file="diff_genes_general_day-C-T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig1=subset(re1,re1$padj<0.05)
sig2=subset(sig1,abs(sig1$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig2, file = "sig2_deseq2_day1vs0_C_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig2,file="diff_genes_day1vs0_C_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig2=subset(re2,re2$padj<0.05)
sig3=subset(sig2,abs(sig2$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig3, file = "sig2_deseq2_day1vs0_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig3,file="diff_genes_day1vs0_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig3=subset(re3,re3$padj<0.05)
sig4=subset(sig3,abs(sig3$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig4, file = "sig2_deseq2_day2vs0_C_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig4,file="diff_genes_day2vs0_C_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig4=subset(re4,re4$padj<0.05)
sig5=subset(sig4,abs(sig4$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig5, file = "sig2_deseq2_day2vs0_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig5,file="diff_genes_day2vs0_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig5=subset(re5,re5$padj<0.05)
sig6=subset(sig5,abs(sig5$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig6, file = "sig2_deseq2_TvsC_day0_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig6,file="diff_genes_TvsC_day0_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig6=subset(re6,re6$padj<0.05)
sig7=subset(sig6,abs(sig6$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig7, file = "sig2_deseq2_TvsC_day1_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig7,file="diff_genes_TvsC_day1_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig7=subset(re7,re7$padj<0.05)
sig8=subset(sig7,abs(sig7$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig8, file = "sig2_deseq2_TvsC_day2_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig8,file="diff_genes_TvsC_day2_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig8=subset(re8,re8$padj<0.05)
sig9=subset(sig8,abs(sig8$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig9, file = "sig2_deseq2_int_1vs0_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig9,file="diff_genes_int_1vs0_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig9=subset(re9,re9$padj<0.05)
sig10=subset(sig9,abs(sig9$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig10, file = "sig2_deseq2_int_2vs0_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig10,file="diff_genes_int_2vs0_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig10=subset(re10,re10$padj<0.05)
sig11=subset(sig10,abs(sig10$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig11, file = "sig2_deseq2_LOP_TvsC_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig11,file="diff_genes_LOP_TvsC_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig11=subset(re11,re11$padj<0.05)
sig12=subset(sig11,abs(sig11$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig12, file = "sig2_deseq2_LOP_0_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig12,file="diff_genes_LOP_0_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig12=subset(re12,re12$padj<0.05)
sig113=subset(sig12,abs(sig12$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig13, file = "sig2_deseq2_LOP_1_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig13,file="diff_genes_LOP_1_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig13=subset(re13,re13$padj<0.05)
sig14=subset(sig13,abs(sig13$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig14, file = "sig2_deseq2_LOP_2_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig14,file="diff_genes_LOP_2_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig14=subset(re14,re14$padj<0.05)
sig15=subset(sig14,abs(sig14$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig15, file = "sig2_deseq2_CAS_TvsC_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig15,file="diff_genes_CAS_TvsC_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig15=subset(re15,re15$padj<0.05)
sig16=subset(sig15,abs(sig15$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig16, file = "sig2_deseq2_CAS_0_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig16,file="diff_genes_CAS_0_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig16=subset(re16,re16$padj<0.05)
sig17=subset(sig16,abs(sig16$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig17, file = "sig2_deseq2_CAS_1_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig17,file="diff_genes_CAS_1_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig17=subset(re17,re17$padj<0.05)
sig18=subset(sig17,abs(sig17$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig18, file = "sig2_deseq2_CAS_2_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig18,file="diff_genes_CAS_2_T_crubrum.csv") #save as cvs file

#extract data which pass an adjusted p value threshold (0.05)
sig18=subset(re18,re18$padj<0.05)
sig19=subset(sig18,abs(sig18$log2FoldChange)>0.3) #Log2 --> normalization of the data to minimize differences between samples due to small counts.

save(sig19, file = "sig2_deseq2_int_LOPvsCAS_T_crubrum.RData")#save as object
#load(file = "sig2_C-T.RData")
write.csv(sig19,file="diff_genes_int_LOPvsCAS_T_crubrum.csv") #save as cvs file

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
sig
topVarGenes <- head(order(sig$padj ), decreasing=TRUE , 100)
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

#Normalized counts for a gene with condition-specific changes over time.
cru <- plotCounts(dds_LRT, which.min(res2$padj), 
                  intgroup = c("day","treatment"), returnData = TRUE)
cru$day <- as.numeric(as.character(cru$day))
ggplot(cru,
       aes(x = day, y = count, color = treatment, group = treatment)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()+
  theme_classic()

betas <- coef(dds_LRT)
colnames(betas)

##Heatmap of log2 fold changes for genes with smallest adjusted p value.

topGenes <- head(order(res$padj),50)
mat <- betas[topGenes, -c(1)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)

dev.off()
