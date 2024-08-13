#Part 1 : Pre-processing




## Installation and set-up


### Install RStudio
*Running RStudio locally?*
  [**Download RStudio**](http://rstudio.com/download)

*Want to try the latest 'Preview' version of RStudio?*
  [**RStudio Preview version**]((https://www.rstudio.com/products/rstudio/download/preview/))


### Install R
*... if you haven't already!*
The RStudio startup message should specify your current local version of R.
For _e.g.,_ `R v4.4.4.1

- [Download R](https://cran.r-project.org)
- [*Update to the latest version of R?*](https://www.linkedin.com/pulse/3-methods-update-r-rstudio-windows-mac-woratana-ngarmtrakulchol)
- [A Windows-specific solution to updating R](https://www.r-statistics.com/2015/06/a-step-by-step-screenshots-tutorial-for-upgrading-r-on-windows/)


### Install packages

**Trouble installing BiocManager?**

- Check your R/RStudio versions and readjust it
- Install individual constituent packages 


```{r}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install('biomaRt')
BiocManager::install("ComplexHeatmap")
install.packages('circlize')
install.packages('openxlsx')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('FactoMineR')
install.packages('plotly')
```


###Load required libraries

```{r}
library(edgeR)                  # Normalization, dispersion estimation
library(openxlsx)               # Importng and exporting data to excel file
library(ggplot2)                # Data visualization
library(dplyr)                  # Data manipulation
library(FactoMineR)             # Performing PCA
library(circlize)               # Visualizing genomic data
library(ComplexHeatmap)         # Creating complex heat maps
library(plotly)                 #Creating interactive graphs
library(biomaRt)                #Retrieve genomic information
```


## Data import


```{r}

# Downloading sample transcriptomics dataset from NCBI's  site
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205748/suppl/GSE205748%5Fread%5Fcounts.csv.gz"


```



###Set directory

```{r}
#setwd()
```
###Check your current directory

```{r}
getwd()
```

## Knowing your data
*Dataset details*:
  
  There are total 27 samples. It contains read counts of RNA sequencing data. Each sample's classified into three categories:

- HC      : Healthy control
- PSA_L   :Psoriatic arthritis skin leison
- PSA_U   :Psoriatic arthritis skin leison uninvolved


##Reading raw expression

##Read count from RNA sequencing data of normal and Psoriatic Arthritis tissue samples.
```{r}
expr0 <- read.table(file = 'GSE205748_read_counts_PsA.csv', 
                    header = TRUE,
                    row.names = 1)
```


##Reading metadata (sample information)

##Metadata containing other information such as sample type, sample title, geo accession  and tissue type
```{r}
meta1 <- read.table('GSE205748_series_matrix_edit.txt',
                    header = TRUE)
```



## Data structure
```{r}
str(expr0)              # Structure of the dataframe
head(expr0)             # Shows the top few observations (rows) of your dataframe
glimpse(expr0)          # Info-dense summary of the data
View(head(expr0, 100))  # View data in a visual GUI-based spreadsheet-like format

colnames(expr0)         # Column names
nrow(expr0)             # No. of rows
ncol(expr0)             # No. of columns


```


##selecting rows and columns

#dataframe[row,column]
```{r}
meta1[1,] #first row
meta1[,1] #first column
```

# Part 2: Data cleaning and preparation


##Checking if column order of expr table is the same as row order of meta1

```{r}
identical(colnames(expr0), meta1$Sample_code)
```


##Since exist difference, so need to rearrange
##Rearranging expression table column order to match sample information order

```{r}
expr1 <- expr0[,meta1$Sample_code]
#rm(expr0) #removing expr0
expr1[ , 'HC_8'] # choose specific one and check with expr table
```


##Checking again if column order of expr table is the same as row order of meta1

```{r}
identical(colnames(expr1), meta1$Sample_code)

```


##Check for missing values

```{r}
table(is.na(expr1))
```


###since there are no missing values and we matched the columns of the expression data and metadata



#Part 3: Normalization

**What is normalization?**

 -Normalization methods minimize the systematic variation of differences between samples   or within samples, gene length and GC content in a raw read counts that affect    differential expression analysis. 




##Categorical variables to be factors that ensure the models treat these variables as discrete groups rather than continuous numeric values

```{r}
factor(meta1$Tissue_type)
```


## Data must be grouped based on experimental conditions to build matrix

```{r}
group1 <- as.factor(meta1$Tissue_type)
group1
```


##Creating the model matrix - With an intercept term (design 1) or without an intercept term (design2) 

```{r}
design1 <- model.matrix(~group1)
design1

design2 <- model.matrix(~group1+0)
design2
```

##Compare matrices

```{r}
head(design1,3)
head(design2,3)
```


##Matrix choose group with no intercept (this only changes the way you define your comparisons later at the contrasts)


```{r}
design <- design2
head(design,3)
rm(design1,design2)



##Create differential gene expression object contain gene expression data as a count matrix and groups as a biological condition we are comparing

```{r}
d <- DGEList(counts=expr1, group=group1 )
```

##Calculate normalization factors. Essential step for normalization!
##Normalizes to adjust the differences in sequencing depth and RNA composition between samples

```{r}
d <- calcNormFactors(d)
```



##Normalize counts

 **what is read count?** 
 - A read count for a  transcript is the total number of reads that align to the exonic regions of that  transcript
 - Read counts are directly proportional to the abundance of the corresponding RNA molecules in the sample
 
 

 **Normalize method**
 
 ***Counts Per Million (CPM)***

-Calculating Counts Per Million (CPM) is a normalization method used in RNA-seq data analysis to account for differences in library sizes (total read counts) across samples -CPM normalizes the read counts by the total number of reads in the library and then scales the counts to a per-million reads basis


###Visualizing the normalized counts with boxplots


###Check current expression, las = 2 turns x axis labels vertically

```{r}
boxplot(cpm(d), las=2)   #cpm(d) function transforms the raw counts in the d object into counts per million
                         #the las=2 parameter rotates the axis labels to be perpendicular to the axis for better readability

```


###Check with log values to eliminate the effect of outliers

```{r}
boxplot(cpm(d, log=TRUE), las=2)

```

**Boxplot**
-Most genes wihtin each sample are expressed very low


##Filtering lowly expressed genes 

-Keep genes that have more than one count per million (cpm) for at least 5% of samples

# 1) cpm function normalizes counts
# 2) cpm(d) >1 checks if the gene expression is more than 1 cpm
# 3) rowSums adds the number of TRUE per Row (aka the number of samples with more than 1 cpm)
# 4) dim(d) return the dimensions of the dataset; dim(d)[1] is the number of genes, dim(d)[2] is the number of samples
# 5) ceiling rounds the given number to the closest higher integer



```{r}
keep <- rowSums( cpm( d ) > 1 ) >= ceiling(0.05*dim(d)[2])
head(keep)
```

###Keep only rows with TRUE
```{r}
d1 <- d[keep, ]
```



###Check how many genes were kept


```{r}
dim(d1)[1] #Number
(dim(d1)[1]/dim(d)[1])*100 #Percentage
```



##Create a new DGElist object with the new counts

```{r}
d2 <- DGEList(counts = d1, group = group1)

```


###Recalculate normalization factors

```{r}
d2 <- calcNormFactors(d2)
```


###Get the normalized counts

```{r}
norm_d <- cpm(d2)

```


### Get the log values of the data

```{r}
logCPM <- cpm(d2, log=TRUE)

```


###Export your normalized values

```{r}
write.table(norm_d, file = "GSE205748_cpm.txt", sep='\t')
write.table(logCPM, file = "GSE205748_logcpm.txt", sep = '\t')

```



***

#Part 4 - Exploratory data analysis

##Create boxplot

###Create a directory to save your plots

```{r}
dir.create("Plots")

```


###Export boxplots to pdf file
#open file

```{r}
pdf(file = "Plots/GSE205748_boxplots.pdf", width = 16, height = 8)
```


###Convert graphics window to 2 columns
```{r}
par(mfrow=c(1,2))

```


#Check sample count distribution again with boxplot

```{r}
boxplot(norm_d, las=2)

```

##With log values 


```{r}
boxplot(logCPM, las=2)

```


##Close graphics device to save pdf

###To check the current active device. If it returns 1, it means there is no active plotting device

```{r}
dev.cur()

```


###If there is no active device, open a new one using windows() on Windows, x11() on Linux, or quartz() on macOS.

```{r}
quartz()
```

###Create a plot to ensure the device is active.

```{r}
plot(1:10)
```

###Use a loop to close all active devices.

```{r}
while(dev.cur() > 1) {
    dev.off()
}

```


###safely close the device using dev.off()

```{r}
dev.off()

```


###Convert graphics window back to 1 column

```{r}
par(mfrow=c(1,1))

```


***


##Principal Component Analysis

**What is PCA?**

-Principal component analysis (PCA), is a statistical technique that reduces the number of dimensions in large datasets into a smaller one that contains most of the information called principal components. The first principal component (PC1) have the highest or most variance then the second principal component (PC2)


###To perform the PCA we need to use the logCPM values and transpose the table

```{r}
tlogcpm <- t(logCPM)    #logCPM has normalized values in log form 

```


###Performing the PCA

```{r}
pcacpm <- PCA(tlogcpm, 
              scale.unit = T, #We scale for better PCA results - TRUE ensures that each gene contributes equally to the principal components, regardless of its variance
              graph = F) #False as we will create our own graphs
```


###Isolate the data for plotting 

```{r}
data_pca <- as.data.frame(pcacpm$ind$coord)    #PCA co-oridinate of each sample from index 

```


###Isolate variance explained by each PC

```{r}
head(pcacpm$eig) #we need the second column
class(pcacpm$eig) # this is a matrix, not a data.frame
pca_perc <- pcacpm$eig[,2] #choosing the second column
```


###Generalize the plotting function

```{r}
PCx=1 #Change accordingly
PCy=2 #Change accordingly

```



### Plotting PCA

```{r}
g <- ggplot(data_pca, 
       aes(data_pca[, 1], data_pca[,2],
           color = meta1$Tissue_type)) + 
  geom_point(key_glyph = "point") + 
  labs(title = 'PCA - GSE205748', 
       x = paste0("PC1(", round(pca_perc[1], 2), ")"), 
       y = paste0("PC2(", round(pca_perc[2], 2), ")"),
       color = 'Tissue type') + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5))

```



###Export to pdf
###Open file

### If R is unable to create or open the specified PDF file for saving the plot
###create the directory if it doesn't exist using the following command before generating the plot
###Ensure that the file name and directory names do not contain special characters that might cause issues with file creation

```{r}
# Ensure the directory exists
if (!dir.exists("Plots")) {
  dir.create("Plots")
}

# Create the PDF plot
pdf(file = "Plots/PCA_GSE205748.pdf", width = 8, height = 6)

# Generate the plot
plot(1:10)

# Close the PDF device
dev.off()


```


###Export to pdf

```{r}
pdf(file = 'Plots/PCA_GSE205748.pdf', width = 8, height = 6)
#print the plot
print(g)
#close the file
dev.off()

```

#Part 5 - Differential expression analysis



##Part 5.1 - Dispersion and model fit

-Dispersion measure the variability or diversity in gene expression levels across different samples

### Using a common estimate across all genes.
-d3 <- estimateGLMCommonDisp( d2, design, verbose=TRUE) 

### Fitting an estimate based on the mean-variance trend across the dataset, such that genes similar abundances have similar variance estimates (trended dispersion)

-d3 <- estimateGLMTrendedDisp(d2, design) 

###Computing a genewise dispersion (tagwise dispersion). Needs one of the former two as prerequisite. Best for multifactorial analysis

-d3 <- estimateGLMTagwiseDisp(d2, design) 


###Performs all 3 of the above dispersion estimations

```{r}
d3 <- estimateDisp( d2, design, verbose=TRUE)

```


###Fit your model

```{r}
fit <- glmQLFit(d3, design)

```


##Part 5.2 - Contrast define and DEA


###Check design matrix categories

```{r}
head(design,3)

```


##Create Contrast parameter
#PsA Lesional vs Healthy would be PsA_les - Healthy, therefore based on the column order of the design matrix:

```{r}
contr <- c(-1,1,0) # this corresponds to group1Healthy*-1 + group1PsA_les*1 + group1PsA_uninv*0, therefore PsA_les - Healthy

```


### We can also make the contrast parameter using makeContrasts if that is convenient.

```{r}
contr2 <- makeContrasts(group1PsA_les-group1Healthy, levels = design) 

```


##Run differential expression test

```{r}
lrt <- glmQLFTest(fit, contrast = contr )

```


###Check the resulting object

```{r}
str(lrt)

```
###Check the gene results

```{r}
head(lrt$table).  # visualize only table

```

** Table Summary**
  
  -logFC(log Fold Change) indicates the magnitude of change in gene expression between the two conditions  
-Minus value indicates the gene is downregulated in the condition of interest compared to the control, 
-Pluss value is indicate upregulated.

-logCPM(log Counts Per Million) - expression level of the gene across all samples, expressed as log2 of the counts per million reads
-The larger the F value means against the null hypothesis (no differential expression)


###Multiple comparisons adjustment

**Why need**
  -Multiple comparison adjustment control the rate of false positive in gene expression analysis, where thousands of genes are tested for differential expression between conditions (like diseased and healthy),  the likelihood of identifying genes as significant purely by chance increases as the number of comparisons increases. Here's why multiple comparisons adjustment is necessary:



###False Discovery Rate Correction (Benjamini-Hochberg methods)

```{r}
lrt$table$fdr <- p.adjust(lrt$table$PValue, method="BH")

```


**Results Filtering**

  -FDR (false discovery rate) is expected proportion of false positives among the declared significant results
  -FDR threshold of 0.01 is very strict and means that we are accepting only a 1% chance of false positives among the significant results
  
  -logFC of 1 indicates a two-fold change in expression, either up or down (since it is in logarithmic scale, 2^1 = 2).
  -Setting a higher threshold (|logFC| >= 1) is stricter and captures genes with more substantial changes in expression.

#Filtering for fdr <= 0.01 (strict) and absolute logFC >=1 (strict)
top = lrt$table[c(lrt$table$fdr <= 0.01 & abs(lrt$table$logFC)>=1),]

  -FDR threshold of 0.05 is less strict, allowing for a 5% chance of false positives
  -logFC of 0.58 corresponds to approximately a 1.5-fold change (2^0.58 ≈ 1.5)
  -Lower threshold (|logFC| >= 0.58) is less strict and captures genes with more modest changes in expression



###Filtering for fdr <= 0.01 (strict) and absolute logFC >=1 (strict)

```{r}
top <- lrt$table[c(lrt$table$fdr <= 0.01 & abs(lrt$table$logFC)>=1),]

```


###Filtering for fdr <= 0.05 (normal) and absolute logFC >=0.58 (normal)

```{r}
top2 <- lrt$table[c(lrt$table$fdr <= 0.05 & abs(lrt$table$logFC)>=0.58),]

```


###Export results

```{r}
write.xlsx(top, file = "DE_Results_GSE205748_FDR_0_01_logFC1.xlsx", rowNames=TRUE, colNames=TRUE)
write.xlsx(top2, file = "DE_Results_GSE205748_FDR_0_05_logFC0_58.xlsx", rowNames=TRUE, colNames=TRUE)


```


###Save workspace

```{r}
save.image("GSE205748_DE_results.RData")

```



***

#Part 6 - Result annotation 




**BioMart**

  -BioMart is a powerful, web-based data management and analysis tool that integrates data from multiple biological databases, including Ensembl (genome data) and others
  -This integration allows  to query across different types of biological data.
  
  
##Choose biomart version

```{r}
mart <- useEnsembl(biomart = 'ensembl', 
                   dataset = 'hsapiens_gene_ensembl', 
                   version = 105) 
```



##Search for annotation based on ensembl_gene_id identifiers

```{r}
annot <- getBM(filters= "ensembl_gene_id", #which identifier you are using
               attributes= c("ensembl_gene_id", 
                             "description", 
                             "start_position", 
                             "end_position", 
                             "strand", 
                             "hgnc_symbol"), #which attributes you want to collect
               values=rownames(lrt$table), #the names of your genes
               mart= mart) #the mart you defined previously


#getBM(get BioMart) allows to retrieve specific data based on your query
#filters are used to find genes based on their Ensembl gene IDs.
#values means we are requesting data for the genes listed in the differential expression table (lrt$table)

```



##Creating final dataframe

```{r}
lrt$table$ensembl_gene_id <- rownames(lrt$table)  ##add ensenmle_gene_id to lrt$table and add rownames to that new column

Final_version <-  lrt$table %>%
  left_join(annot, by = 'ensembl_gene_id')
head(Final_version)

write.xlsx(Final_version, file="GSE205748_DE_results_PsA_les_vs_healthy.xlsx", colNames=TRUE, rowNames=TRUE)

```



***


#Part 7 - Results Visualization

## 7.1 - Volcano plot

**Volcano Plot**

  -A volcano plot is a type of scatter plot that is used in the analysis of gene expression studies and represents the relationship between the statistical significance and the magnitude of change (effect size) for each data point for a gene
  
  -X-axis represents  log2 fold change in gene expression  between two conditions. A positive value indicates upregulation, and a negative value indicates downregulation

  -Y-axis- represents the negative logarithm (base 10) of the false discovery rate (FDR) p-value.
  - The higher the point on the y-axis, the more statistically significant the change is. Using the negative log10 of the p-value ensures that small p-values (indicating high significance) are represented as large positive values, making them easy to identify

*Thresholds and Cutoffs*

  -Significance Threshold (Y-Axis): A horizontal line is often drawn to indicate the threshold for statistical significance, such as a p-value cutoff of 0.05. Points above this line are considered statistically significant.
  -Fold Change Threshold (X-Axis): Vertical lines may be drawn to indicate thresholds for biologically meaningful fold changes (e.g., log2 fold change > 1 or < -1).




```{r}
ggplot(data=Final_version, 
       aes(x=logFC, y=-log10(fdr))) + 
  geom_point() + 
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") + 
  theme_minimal()



```

###Adding color to our genes according to their significance level

```{r}
DE_res <- Final_version

DE_res <- DE_res %>%
  mutate(significance = case_when(
    logFC >=1 & fdr <= 0.01 ~ 'Upregulated',
    logFC <=-1 & fdr <= 0.01 ~ 'Downregulated',
    abs(logFC) < 1 | fdr > 0.01 ~ 'Not significant'
  ))

```



###Plotting again, this time with color
###Opening file

```{r}
pdf("Plots/VolcanoPlot_GSE205748.pdf", height = 8, width = 8)
#preparing the plot
g <- ggplot(DE_res,
       aes(x = logFC,
           y = -log10(fdr),
           color = significance,
           label = hgnc_symbol) # we will use rownames with plotly below, it is not required here for the volcano
       ) + 
  geom_point() + 
  geom_vline(xintercept = c(-1,1), color = 'red') + 
  geom_hline(yintercept = -log10(0.01), color = 'red') + 
  labs(x = 'LogFC', y = '-Log10(FDR)',color = 'Significance') + 
  scale_color_manual(values = c('blue','grey90','red'))+
  theme_minimal() 
```


###Printing

```{r}
print(g)

```

###Closing the file


```{r}
dev.off()

```



###Interactive volcano plot for better identification of DE significant genes
```{r}
library(plotly)
print(g)
ggplotly()
```

##Part 7.2 - Heatmap


**Heatmap**


  -A heatmap is a graphical representation of gene expression data that is commonly used to visualize the expression levels of a large number of genes across multiple samples and  conditions, and allows to quickly identify patterns, such as groups of genes with similar expression profiles 
  
  - Rows: Typically represent individual genes or transcripts
  
  -Columns: Represent different samples, conditions, or time points
  
  -Each cell in the heatmap shows the expression level of a particular gene in a specific sample

  -Genes and/or samples are often clustered based on their expression patterns, so that genes with similar expression profiles or samples with similar overall gene expression patterns are grouped together
  


###Select the top 50 upregulated genes

```{r}
Up50 <- DE_res %>%
  slice_max(order_by = logFC, n =50)

```


###Select the top 50 downregulated genes

```{r}
Down50 <- DE_res %>%
  slice_min(order_by = logFC, n =50)
```


###Bind dataframes

```{r}
Top100 <- bind_rows(Up50, Down50)

```


###Isolate the expression of top 100

```{r}
logCPM_100 <- logCPM[Top100$ensembl_gene_id,] 

```


###Save heatmap to a file

```{r}

pdf('Plots/GSE205748_heatmap_top100_DE.pdf', width = 8, height =8)

```


#Create the heatmap

```{r}
Heatmap(logCPM_100,
        row_labels = Top100$hgnc_symbol, 
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7), 
        heatmap_legend_param = list(title = "LogCPM\nexpression"), 
        top_annotation = HeatmapAnnotation(Condition = meta1$Tissue_type, 
                                           which = 'column', 
                                           col = list(Condition = c(PsA_les = 'turquoise4',
                                                                    PsA_uninv = 'red3',
                                                                    Healthy = 'green2')
                                                      )
                                           )
        )

```


###Close the file

```{r}
dev.off()

```


###Save final workspace
```{r}
date <- Sys.Date()
save.image(file = paste0(date,"_GSE205748_DE_analysis_complete.RData"))

```



###Save session info to a text file for reproducibility purposes

```{r}
sink(file = 'session_info_GSE205748_DE_analysis_complete.txt')
sessionInfo()
sink()



```








