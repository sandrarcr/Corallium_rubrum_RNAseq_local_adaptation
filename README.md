## README

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

This repository contains the analysis for the paper:
**"Does local adaptation influence thermal responses in red coral populations across depth gradients? Transcriptomic insights for effective conservation."**

---

## **Overview**
This project investigates the transcriptomic responses of the precious coral *Corallium rubrum* across depth gradients. It provides insights into local adaptation by analyzing differential gene expression between populations from contrasting environments (shallow vs mesophotic zones) subjected to thermal stress. In this reposity we included the scripts to conduct Differential Gene Expression Analysis. Associated data to this study can be found in ENA accession number: 

---

## **Methodology**
1. **Sample Collection**:
    - Shallow Population (15–18m): CAS (N=3) sampled across 3 time points (days 5, 10, 15) – 9 replicates.
    - Deep Population (48m): LOP (N=3) sampled across 3 time points (days 5, 10, 15) – 9 replicates.
      
<img width="669" height="540" alt="Figure 1 - paper 2" src="https://github.com/user-attachments/assets/45d13429-cfa6-4373-a7c4-83e4a0870106" />


2. **RNA-seq Workflow**:
    - RNAseq quality control and pipeline details were conducted using FASTQC, Trimmomatic, Kraken, Trinity, Bowtie2, RSEM, BLAST and BUSCO. Details can be found in the paper.
    - Differential expression analysis using DESeq2.
    - Statistical tests:
        - Likelihood Ratio Tests (LRT) for group-level differences.
        - Wald tests for pairwise comparisons.

3. **Visualization**:
    - We followed the current DESeq2 manual vignette for plots.
  



