## Tutorial 3

**Objective:** Perform a differential expression analysis with DESeq2 in R.


## Download and install R on your computer 

If you don't have R on your computer, please install it: [https://www.r-project.org/](https://www.r-project.org/)

R is freely available software for statistical analyses and figure generation.
R packages are bundles of code that perform myriad operations and analyses and can be installed from a CRAN mirror.  
Bioconductor is essentially a project for developing and maintaining code, written in R, for analysis of biological data.
We will install [Bioconductor](https://www.bioconductor.org/) as well as three of its packages: 
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
[tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

Open an R window on your computer and install Bioconductor, DESeq2, tximport, and clusterProfiler.

	install.packages("BiocManager")
	BiocManager::install("DESeq2")
	BiocManager::install("tximport")
	BiocManager::install("clusterProfiler")
	
Please also install these other useful packages:

	install.packages("ggplot2")
	install.packages("pheatmap")

Now load the packages in your R session

	library(DEseq2)
	library(tximport)
	library(clusterProfiler)
	library(ggplot2)
	library(pheatmap)

<br>

## Import your Salmon sf files and summarize transcript read counts for a gene-level analysis.

To perform a DESeq2 analysis we will need our Salmon sf files, our `tx2gene.txt` file, and a table that specifies the conditions associated with our samples (the experimental design).
You can manually create the table in R, Excel (exporting as a tab-delimited text file), or in any text editor, such as `nano`. 
We'll create this file after we import our sf and tx2gene files.

If you haven't transferred your files from the VM to your computer, you can obtain copies of them here: [data-files](https://github.com/elasekness/RNASeq_workshop/tree/main/data-files)

Import the `sf` and `tx2gene.txt` files with the `tximport()` function. 

	files = file.path(c("wt-1.sf", "wt-2.sf", "pqse-1.sf", "pqse-2.sf"))
	names(files) <- c("wt1", "wt2", "pqse1", "pqse2")
	tx2gene = read.table("tx2gene.txt", header=T, sep='\t')
	txi = tximport(files, type="salmon", tx2gene=tx2gene)

> The first command specifies the path to our files. <br>
> The second command names the files.  **These names should be in the same order as listed in our experimental design table** <br>
> We import the `tx2gene.txt` file as a table, specifying that it has a header and columns/fields separated by tabs. <br>
> Finally, we create the txi object, which will summarize transcript-level abundances to the gene-level. <br>
> In our case, the abundances will not change because the transcripts are equivalent to the genes. <br>
> **Note** including the option `countsFromAbundance=“lengthScaledTPM”` in the tximport command removes read count correlation with length as needed for the DESeq2 analysis. However, the `DESeqDataSetFromTximport()` function will do this for us when we create our DESeqDataSet object.

The `txi` object contains the matrices for abundance (TPM), counts, and length from our sf files. You can get the names of the matrices and access the information in the matrices with the following commands.

	names(txi)
	head(txi$counts)
	
> `head` returns the first few lines of the counts table, just as it does on our VMs. <br>


## Create a table that specifies the experimental design.

Now create a design table that specifies the condition(s) associated with each file.  Again, ensure that the file names are listed in the same order as they were imported.
This file will tell DESeq2 how the data should be analyzed.

	ColData = data.frame(row.names=1, c("wt1", "wt2", "pqse1", "pqse2"), Condition=c("wt", "wt", "pqse", "pqse"))

> `pqse` is our treatment that will be compared against our wild-type, `wt` control. <br>
> `row.names=1` assigns that the first column of data (gene names) as the row names.  Without this, R would create numerical row names. <br>


## Create the DESeqDataSet object.

We will load the data from the txi object and the `ColData` table into the DESEqDataSet object.  
We will also specify that the design formula can be found in the column entitled `Condition` in  our `ColData` table.

	dds = DESeqDataSetFromTximport(txi, ColData, ~Condition)
	
> **Note** raw count data can be imported from a table as well: `dds <- DESeqDataSetFromMatrix(countData = counts, colData = ColData, design = ~ Condition)`.


Our DESeqDataSet object has specific slots for the different data we provide, including the original count data and the design of the study.

	head(counts(dds))
	design(dds)

## Normalize the read counts.

DESeq2 normalizes read counts with the median of ratios method, which accounts for variable sequencing depth and RNA composition. Briefly, the geometric mean for each gene is calculated across all samples and used to generate a (gene read count)/(geometric mean) ratio for each gene in each sample.
The normalization factor or size factor for each sample is the median of that sample's ratios, which should eliminate the influence of differentially expressed genes.  Normalized counts are generated by dividing raw counts of a sample
by its normalization or size factor.


We need only one function call to calculate the size factors and perform the normalization.

	dds <- estimateSizeFactors(dds)

> Normalized read counts are now saved back to the dds object. <br>

We can retrieve the size factors with the following command:

	sizeFactors(dds)
	
We can also retrieve the normalized read counts and save them to a table.

	head(counts(dds, normalized=T))
	write.table(counts(dds, normalize=T), file="normalized_readcounts.txt", sep='\t', quote=F, col.names=NA)
	
> **Note** DESeq2 uses the raw read count data for DE analysis and models the normalization inside the Generalized Linear Model (GLM). <br>
> The normalized read counts we generated are useful for visualization purposes.

## Visually evaluate the quality of your data.

Ideally, we would like to see that our biological replicates look very similar to each other and that the variation across treatments is greater than within treatments.
PCA (Principal Component Analysis) reduces multi-dimensional data into the principal components that explain the most variation within our dataset.  The first
two principal components are plotted against each other for each sample. Samples will appear close together in the plot if they share similar expression levels 
for genes that contribute most to the variation in PC1 and PC2.

	rld<- rlogTransformation(dds, blind=TRUE)
	
> This function log-transforms our data for better separation/clustering of our data. <br>
> Setting the `blind` option to `TRUE` ensures that the calculation is not influenced by sample information. <br>

	plotPCA(rld, intgroup="Condition")
	
> We are coloring by the values in `Condition`.  If we had additional variables, such as sex, time, or stain, we could include these in our design and color the samples accordingly to determine which variables best explain the variation we see in PC1 and PC2.

* Based on your PCA plot, do you think we have good replicates?
* How do you think this will influence the number of differentially expressed genes we are able to detect?


You could also plot the normalized read counts of replicates against each other.
Ideally, all points would be close to or on the 45 degree line, indicating a 1:1 ratio of expression between replicates.

	normalized_counts = counts(dds, normalize=T)
	plot(normalized_counts$wt1, normalized_counts$wt2)
	
Or we can make a heatmap showing the correlation of expression values across samples.

	rld_mat <- assay(rld)
	rld_cor <- cor(rld_mat)
	pheatmap(rld_cor, annotation = ColData)
	
> The first command retrieves the matrix of log-transformed values from the rld object. <br>
> The second command generates the pairwise correlations among samples. <br>
> The third command uses the pheatmap() function from the pheatmap package to make a heatmap, adding annotations from our ColData table. <br>

## Perform the differential expression analysis.

This requires a single function that is really performing multiple steps, including estimating size factors and gene-wise dispersion, and fitting the data to the model.

	dds <- DESeq(dds)
	
> DESeq2 uses a negative binomial distribution to model the read count data, which accounts for the fact that variance > mean (overdispersion).

DESeq2 uses the Wald test to determine whether the expression changes observed across treatments are statistically significant 
(the null hypothesis being that there is no change in expression between two conditions as measured by the log-fold change).

	res = results(dds, contrast=c("Condition", "pqse", "wt"))
	write.table(res, file="pqsE-vs-wt-results.txt", sep='\t', quote=F, col.names=NA)
	
> The `results()` function performs the hypothesis testing for us, comparing `pqse` to `wt`.
> **Note** the order of the sample names in the contrast is important. 
> The `res` object can be saved as a table for further inspection.

The output includes:
* `baseMean` mean of normalized counts for all samples
*`log2FoldChange` log2fold change (in our case, a negative number indicates a decrease in expression in pqsE with respect to the wild-type condition)
* `lfcSE` standard error of baseMean
* `stat` Wald statistic
* `p-value` p-value
* `padj` FDR adjusted p-value to take into account the increase in false positives due to multiple testing

Please refer to this excellent tutorial for additional explanation of the differential expression analysis with DESeq2: [hbc_training_DGE_analysis_workshop](https://hbctraining.github.io/DGE_workshop_salmon/schedule/)
