
# Differential Expression Analysis (DEA) Using edgeR


### What is DEA?

Differential expression analysis is used to discover quantitative changes in expression levels between experimental groups by taking normalized read count data to identify genes that show statistically significant differences. 

### Why edgeR?

edgeR is a Bioconductor software package for examining differential expressed genes of RNAseq count data between different groups. edgeR is based on negative binomial (NB) distributions. It is important to consider the experimental design when choosing an analysis method. edgeR can perform multiple comparisons while some of the differential expression tools can only perform pairwise comparison.


## Objective


A differential expression analysis of RNA sequencing data by using edgeR performed  to compare samples from Psoriatic Arthritis patients to Healthy controls in order to identify deregulated genes whose expression levels are significantly different that helps to understand molecular mechanism of Psoriatic Arthritis, can lead to identify potential biomarkers for therapeutic targets. 


##  Sample

There are total 27 samples that contains read counts of RNA sequencing data. Each sample's classify into one of the three categories:

- Healthy
- PsA_les (Psoriatic Arthritis Lesional)
- PsA_uninv (Psoriatic Arthritis Uninvolved - areas that are not affected by Psoriatic Arthritis)



Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205748


## Main Workflow
![Snip20240814_28](https://github.com/user-attachments/assets/20ca21b2-ee87-4634-9cd8-8207be6c6838)



## Overview

Part 1: Pre-processing

- Installation and Setup
- Loading packages
- Data import
- Knowing the data structure

Part 2: Data cleaning and preparation

- Rearrange data
- Handle missing values

Part 3: Normalization

- Factor
- Model Matrix
- Normalize counts
- Barplot
- Filter

Part 4: Exploratory data analysis
- Boxplot
- PCA 

Part 5: Differential expression analysis

- Dispersion and model filter
- Differential test
- Multiple comparison
- Filtering
- Export result

Part 6: Results annotation

- Biomart
- Final data frame

Part 7: Result visualization

- Volcano plot
- Heatmap
- Save final workspace


## References

1. Johnsson, H., Cole, J., Siebert, S., McInnes, I. B., & Graham, G. (2023). Cutaneous lesions in psoriatic arthritis are enriched in chemokine transcriptomic pathways. Arthritis Research & Therapy, 25(1), 73. https://doi.org/10.1186/s13075-023-03034-6

2. Johnsson, H., Cole, J., McInnes, I. B., Graham, G., & Siebert, S. (2024). Differences in transcriptional changes in psoriasis and psoriatic arthritis skin with immunoglobulin gene enrichment in psoriatic arthritis. Rheumatology, 63(1), 218-225.
https://doi.org/10.1093/rheumatology/kead195

