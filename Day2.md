## Tutorial 2

**Objective:** Quantify transcript abundance from our RNASeq data.

RNA-Seq analyses are often employed to determine the differential gene expression that occurs between two conditions.
Two common strategies to quantify gene expression are 1) to map reads to a reference genome and count the reads that overlap genic coordinates or 
2) to quantify the reads that map to transcripts. When no reference assembly is available, a de novo transcriptome assembly is generated and reads are mapped back to it for quantification purposes. 

Here we are studying the effects of deleting the pqsE gene on overall gene expression in _Pseudomonas aeruginosa_ (data generously provided by the Paczkowski lab).
PqsE directly interacts with the RhlR transcription factor, which is involved in the quorum sensing cascade and in activating genes required for
pathogenesis and biofilm formation. The Paczkowski lab demonstrated that PqsE complexed with RhlR enhanced binding affinity of RhlR target promoters and
increased transcription of genes encoding virulence factors associated with pathogenesis.

Our questions are:

* What genes are significantly up or down-regulated by the deletion of the pqsE gene?
* What pathways are significantly enriched for differentially expressed genes?

Our read data are in 150x150 bp paired-end format (each replicate as an R1 and R2 read file) and were sequenced on an Illumina NextSeq machine.
We have RNASeq results from two biological replicates of a pqsE deletion strain and two biological replicates of the wild type strain.

The number of replicates included is an important consideration of your experimental design - increasing the number of replicates increases the power to detect significant changes.  
A good reference regarding the number of biological replicates to consider is [Schurch et al 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4878611/).


  
<br>


## Download the fastq files for our RNA-Seq analysis from the SRA

This step is unnecessary since we already have the fastq files.  However, if these data had been submitted to NCBI's SRA (Sequencing Read Archive) database,
we could download the files directly from the SRA to our VM using the **`prefetch`** and **`fasterq-dump`** commands.

<br>

## Get some basic read statistics

You can generate a summmary of the quality of your data with [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), which provides nice graphics or we can use [SeqKit](https://bioinf.shenwei.me/seqkit/) for a more pared-down summary of our data.

	cd fastq/
	seqkit stats wt-1-R1.fastq.gz
 	seqkit stats wt-1_R2.fastq.gz
 	
Or you can run stats on both files at the same time

	seqkit stats wt-1*gz -T | csvtk pretty -t

<br>

## Clean your reads with [TrimGalore](https://github.com/FelixKrueger/TrimGalore).

You can process your PE fastqs one-by-one or you could execute a BASH for-loop to do the job for you.

The long way:

	trim_galore -q 30 --length 100 --trim-n --paired wt-1_R1.fastq.gz wt-1_R2.fastq.gz

> TrimGalore will automatically detect the sequencing adapter and trim it from our reads. <br>
> We also trim reads to a PHRED quality score of 30 (1/1000 chance of being a miscalled base), remove ambiguous bases, and only retain reads with a minimum length of 100 bp. <br>
> The **`--paired`** option keeps R1 and R2 reads together, which is necessary for most (all?) mapping software. <br>
> Other popular read trimming programs include [fastp](https://github.com/OpenGene/fastp) and [Trimmomatic](https://github.com/usadellab/Trimmomatic) <br>

The short way:

	ls *R1.fastq.gz | cut -d "_" -f 1 > seqlist
	for filn in `cat seqlist`; do trim_galore -q 30 --length 100 --trim-n --paired $filn"_R1.fastq.gz" $filn"_R2.fastq.gz"; done

> We list all R1 fastq files and cut them on the underscore delimiter to take the first field, which is the base name for each PE fastq file (For example, `wt-1`). <br>
> We save that output to a file called `seqlist`.  You can **`cat`** the seqlist file to see the results. <br>
> We loop through the `seqlist` file line-by-line to obtain the basename of each PE library (instead of looping through a series of files). <br>
> The back ticks represent a subroutine. The output of the subroutine command `cat` is being passed to the for-loop.  <br>
> Thus, `$filn` is a variable that will represent each basename in our `seqlist` file. We use this basename to specify the R1 and R2 fastq files by filling in the reminder of the unique part of each PE file's name. <br>
> For example, `$filn` will be interpreted as `wt-1` and the file endings `"\_R1.fastq.gz"` and `"\_R2.fastq.gz"` will be interpreted literally because of the quotation marks, which gives us the full file names: `wt-1_R1.fastq.gz` and `wt-1_R2.fastq.gz`.

<br>


## Quantify transcripts with [Salmon](https://combine-lab.github.io/salmon/)

[Salmon](https://combine-lab.github.io/salmon/) is a 'wicked fast' program that allows the direct quantification of reads against a transcriptome (no need for an initial read mapping step). 
We'll use the coding sequences that we downloaded as our transcripts.

To make typing downstream commands easier, let's move our reference coding sequence fasta file and our trimmed read files to a new directory.

	cd
	cd rnaseq_workshop/
	mkdir salmon_analyses
	mv fastq/*fq.gz salmon_analyses
	mv reference_db/pa14_cds.fna salmon_analyses

> Typing **`cd`** by itself returns you to your home directory (as you might have learned in tutorial 1). <br>

Index the coding sequence file.

	salmon index -t pa14_cds.fna -i salmon_index
	
> If you **`ls`** your directory, you'll see that salmon has put the index files to your 'transcriptome' in a sub-directory called `salmon_index`. <br>

Quantify transcript abundance.

	for filn in `cat seqlist`; do salmon quant -i salmon_index -l A -1 $filn"_R1_val_1.fq.gz" -2 $filn"_R2_val_2.fq.gz" --validateMappings -o "$filn; done
	
> The **`quant`** function of **`salmon`** allows direct quantification of reads agains the 'transcriptome' index and will output the results into a directory called `salmon_quant`. <br>
> The tab-delimited output files (sf files) for each library will be in a directory named for the library basename.  <br>
> The **`-l`** option specifies the automatic detection of the library type.  The **`--validateMappings`** option is the recommended default.  It essentially checks that the mappings are plausible enough to be quantified. <br>


The `Salmon` output files required for downstream analyses are the `quant.sf` files. 
The sf file is a tab delimited text file containing the length and effective length of each transcript (effective length relating to the expectation of sampling more or less reads from a transcript), the normalized TPM (transcripts per million) values, and read counts. The TPM values are referred to as pseudocounts and need to be non-normalized for DESeq2 analyses.

All of our sf files are currently named `quant.sf`. Let's use a for-loop to rename them according to the basename of the fastq files.

	for filn in `cat seqlist`; do mv $filn"/quant.sf" $filn".sf"; done

<br>

## Generate a table that associates transcripts with genes.

If we had true transcript data we would need to collapse our transcript-level abundances to gene-level abundances.
The first step is to create a tab-delimited table that associates transcripts to genes (with transcript ID being placed in the first column). 

In our case, transcripts and genes are the same so we will create a table with the locus tag listed twice:

	TXNAME	GENEID
	PA14_00010	PA14_00010
	
> The header names can be anything. <br>

We can generate this table with various methods, including a combination of **`grep`** and **`sed`** commands.
Navigate to the location of your reference coding sequence file, which should now be in the `salmon_analyses` directory.

	cd ~/rnaseq_workshop/salmon_analyses
	grep ">" pa14_cds.fna | sed "s/>//" | sed "s/\(.*\)/\1\t\1/" > tx2gene.txt
	
> The tilde **`~`** in the **`cd`** command is a shorthand way of specifying your home directory. <br>
> The **`grep`** command grabs all the definition lines from our reference transcript file. <br>
> The first **`sed`** command removes the **`>`** symbol, as it is not part of our transcript name. <br>
> The second **`sed`** command searches for any number of characters any number of times with **`.*`**. This represents all locus tags, which are saved via the special parenthetical notation. <br>
> We recall the locus tag twice in our replacement command with **`\1`** and separate the transcript name from the gene name with a tab - **`\t`**. <br>
> You can manually add the headers `TXNAME` and `GENEID` in **`nano`** and save the output.

<br>

## Transfer data from VM to poject's storage bucket.

The `sf` files and the `tx2gene.txt` file are necessary for our DEG analyses.  While our files are currently located on the training VMs, we will perform our DEG analyses in **`R`** on our personal computers. We need to transfer our data from the VMs to the project's associated storage bucket.  We can then download the files directly to our computers from the storage bucket.  We will make use of GCP's [gsutil](https://docs.cloud.google.com/storage/docs/gsutil) commands (note that **`gsutil`** is a legacy package and Google now recommends using **`gcloud storage`** commands).

**`gsutil`** is Python application that allows you to interact with cloud storage (the bucket) from the command line. Functions include, copying files between bucket and VM, and deleting, moving, and listing files in storage.

Before we can copy anything to the bucket, we need to authenticate our accounts using a function from the **`gcloud`** application.

	gcloud auth login

> Now follow the instructions provided.

<br>

You can copy data to the bucket once you've successfully authenticated but before we copy anything, let's create a location in the bucket to which you can transfer your files.

Navigate to `Cloud Storage` from the GCP window in your web browser. Click `Google Cloud` in the top left of the window and then the `Cloud Storage` option under `Quick access`. 

Click the `Buckets` option from the menu on the left side of the page, which will show one bucket associated with the project: `wc-bms-bi-training-bucket.` 

Click the link to this bucket to see its contents. Click on the folder entitled, `rnaseq_workshop0126` and then click `create folder` option to create a folder for your data.  The name of the folder should be unique so use your initials.

Return to your terminal window and navigate to the location of your `sf` and `tx2gene.txt` files, if not already there (**`cd ~/rnaseq_workshop/salmon_analyses`**.  Now copy these files to your folder in the bucket.

	gsutil cp *sf gs://wc-bms-bi-training-bucket/rnaseq_workshop0126/your_initials
	gsutil cp tx2gene.txt gs://wc-bms-bi-training-bucket/rnaseq_workshop0126/your_initials

> Notice that the **`gsutil cp`** command is similar in structure to a BASH **`cp`** command. <br>
> Also notice that we can use metacharacters to specify multiple files. <br>
> The `gs://` prefix must be added to the path in your bucket to indicate a resoure in Cloud Storage.

<br>

Now you can click the link to each file in the bucket and choose the `download` option to copy the data to your personal computer.
