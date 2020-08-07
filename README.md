# RNA-seq-Data-Analysis-Work-Flows
This is a hands-on, step-by-step tutorial on how to conduct a basic RNA sequence analysis for absolute beginners!  
![workflow](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows/blob/master/figures/rna_workflow.png)
The picture above is the RNA seq workflow. The following steps in this document will go through this workflow.  

#### Contents  
1. [Download SRA sequences from NCBI](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#1-download-sra-sequences-from-ncbi)  
  a. [Obtain accession numbers from search results](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#a-obtain-accession-numbers-from-search-results)  
  b. [Obtain run accessions](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#b-obtain-run-accessions)  
  c. [Download sequence data files using SRA Toolkit](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#c-download-sequence-data-files-using-sra-toolkit)  
2. [RNA Seq Pipeline 1](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#2-rna-seq-pipeline-1)  
  a. [Setting up](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#a-setting-up)  
  b. [Alignment with STAR](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#b-alignment-with-star)  
3. [RNA seq Analysis: Read count and DESeqDataSet](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#3-rna-seq-analysis-read-count-and-deseqdataset)  
  a. [Start an interactive R Studio](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#a-start-an-interactive-r-studio)  
  b. [Convert RNA-seq reads into counts](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#b-convert-rna-seq-reads-into-counts)  
4. Analysis  
  4.1. [RNA seq Analysis: Differential Expression (DE) Analysis](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#41-rna-seq-analysis-differential-expression-de-analysis)  
  4.2. [RNA seq Analysis: Down Stream Plots](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#42-rna-seq-analysis-down-stream-plots)  
    - a. [Normalized count data for down-stream plots](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#a-normalized-count-data-for-down-stream-plots)  
    - b. [Heatmap](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#b-heatmap)  
    - c. [Volcano plot](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows#c-volcano-plot)  

## 1. Download SRA sequences from NCBI

### a. Obtain accession numbers from search results
- Search by keywords or accessions. https://www.ncbi.nlm.nih.gov/sra  
NCBI is a major database for DNA sequences.  
![pic_01](https://github.com/xinyi15/ChIP-seq-Data-Analysis-Work-Flows/raw/master/figures/pic_01.jpg)  
- Click the checkboxes next to records (experiments) to select data of interest. Leave all checkboxes unchecked to select all records (experiments) from your search.

### b. Obtain run accessions
- Click on _Send to:_, _File_, and select _RunInfo_ to create your file.

### c. Download sequence data files using SRA Toolkit
- Connect to Longleaf.
- Change your current directory to Your /pine directory.  
```
cd /pine/scr/<o>/<n>/<onyen> (the “o/n/” are the first two letters of your ONYEN)
```
- Use fastq-dump command to download the fastq data via bash (changing run accessions to whatever the accessions are) [code](https://github.com/xinyi15/RNA-seq-Data-Analysis-Work-Flows/blob/master/scr/Download_sequence_data_RNA.sh)  
**Remember to replace the code within <>.** To execute the code file type `bash <filename>.sh`. Alternatively, copy and paste the following:
```
#!/bin/bash

#SBATCH --job-name=run_5981
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --ntasks=1
#SBATCH --mem=500g
#SBATCH --mail-type=END   
#SBATCH --mail-user=<Email Address>@email.unc.edu 


module load sratoolkit/2.9.6

for run in <Run Accessions> #eg: SRR5385266 SRR5385267 SRR5385268 SRR5385269 SRR5385270 SRR5385271 SRR5385272 SRR5385273 
do
        fastq-dump $run 
done
```
_SBATCH_ lines are the task execution properties.  
The SRA Toolkit will convert data in SRA format to fastaq format.  

- If the it is paired-end sequencing, we can use the following code to download the data.
```
fastq-dump --split-files
```

&nbsp;
&nbsp;
&nbsp;

## 2. RNA Seq Pipeline 1

### a. Setting up
- Change your current directory to Your /pine directory and create the output folder.  
```
mkdir <Output>
```
- Load packages (The bash code is here: [code file](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows/blob/master/scr/RNA_seq_pipeline_1.sh)  
These packages will respectively convert .fastq files to .sam (human-readable files) and then to .bam (compressed, binary equivalent) files. 
```
module load star
module load samtools
```
- Set downloaded fstq file name as f. e.g:
```
f=SRR5385266 
```
### b. Alignment with STAR
Align reads to the genome using STAR  
- STAR performs a two-step process:
- Seed searching
- Clustering, stitching, and scoring.  
To get a better understanding of the alignment method STAR utilizes, please click [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)
```
 STAR --genomeDir /pine/scr/x/i/xinyi7/SRAdata/mm10_UCSC_index # for this tutorial we will retrieve data from this location
      --runThreadN 12 --readFilesIn ${f}.fastq  --outFileNamePrefix <Output>/${f}.
```
Change file format from SAM to BAM  
- -b: Output in the BAM format.
- -S: Input in the SAM format
```
samtools view -bS <Output>/${f}.Aligned.out.sam -o <Output>/${f}.bam
```
Sort BAM files by genomic coordinates  
```
samtools sort -o <Output>/${f}.sorted.bam <Input>/${f}.bam
```
&nbsp;
&nbsp;
&nbsp;


## 3. RNA seq Analysis: Read count and DESeqDataSet

### a. Start an interactive R Studio
This will open up an R Studio window
```
module add r/3.5.0
module add rstudio
srun --mem=200g -t 5:00:00 -p interact -N 1 -n 1 --x11=first rstudio
```
Setting memory as 200g may be too large, you can adjust that if there are errors.  
In your interactive R Studio window, create a R file and install the following packages  
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("NMF")

library(Rsubread)
library(org.Mm.eg.db)
library ( DESeq2 )
library ( NMF)
```
### b. Convert RNA-seq reads into counts
After letting our reads aligned to the genome, the next step is to count how many reads have mapped to each gene.
- Save the path of each bam file into a vector
```
all_bam_file <-list.files(path = <Path of folder>) #eg: "C:/Users/Dropbox/"
all_bam_file= sapply(all_bam_file,function(x)  {paste0(<Path of folder>,x)})  #eg: "C:/Users/Dropbox/"
```
- Check proportion of mapping
```
## this step takes a lot of time
propmapped(all_bam_file)
```
- We use Featurecounts to report the "raw" counts of reads that map to a single location (uniquely mapping). Essentially, total read count associated with a gene (meta-feature) = the sum of reads associated with each of the exons (feature) that "belong" to that gene. Please click [here](https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/05_counting_reads.md) to see more details.
```
## counts read from all .bam files, which are all using mm10 annotation (same as the reference genome when alignment)
fc_mm10 <- featureCounts(all_bam_file, annot.inbuilt="mm10", isPairedEnd=FALSE) #you can change it to other reference gemone
save(fc_mm10,file=<Path of file>)#eg:"C:/Users/Dropbox/read_count_mm10.RData"
```
- Output The read count table(fc_mm10$counts) should look like this:
![read_count](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows/blob/master/figures/read_count.png?raw=true)
Each row is a geneId and each column is a sample.
- Convert geneId into geneName  
```
data=fc_mm10$counts
a=mapIds(org.Mm.eg.db, rownames(data), keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rownames(data)=a
```
- Use sample names as column names
```
my_counts<-data
#### eg: examples of sample names

#RNA-seq of WT Tconv replicate 1
#RNA-seq of WT Tconv replicate 2

#RNA-seq of Cd4CreSatb1CKO peripheral Treg replicate 1
#RNA-seq of Cd4CreSatb1CKO peripheral Treg replicate 2

#RNA-seq of WT peripheral Treg replicate 1
#RNA-seq of WT peripheral Treg replicate 2

#RNA-seq of Cd4CreSatb1CKO tTreg precursor and tTreg replicate 1
#RNA-seq of Cd4CreSatb1CKO tTreg precursor and tTreg replicate 2

#RNA-seq of WT tTreg precursor and tTreg replicate 1
#RNA-seq of WT tTreg precursor and tTreg replicate 2

####

my_counts[is.na(my_counts)] <- 0

colnames(my_counts) <-  c( paste0("WT_Tconv", 1:2), paste0("Cd4CreSatb1CKO_peripheral_Treg", 1:2), 
paste0("WT_peripheral_Treg", 1:2),  paste0("  Cd4CreSatb1CKO_tTreg_precursor_and_tTreg", 1:2),
paste0(" WT_tTreg_precursor_and_tTreg", 1:2))

sample_names <- colnames(data)[1:10]  ## make a copy of colomn names before renaming them 


sample_info <- data.frame ( condition =sapply(colnames(my_counts),function(x) substr(x,1,(nchar(x)-1))),
                            row.names = colnames( my_counts) )
```

The sample_info table should look like this:  
![sample_info](https://github.com/colachocolat/RNA-seq-Data-Analysis-Work-Flows/blob/master/figures/sample_info.png?raw=true)
Each row is the sample name and the column is condition.  
- Create a DESeqDataSet object to proceed with other analysis
```
rownames(my_counts)=NULL
 #3.1 generate the DESeqDataSet: with Batch effect added 
DESeq.ds <- DESeqDataSetFromMatrix ( countData = my_counts,
                                     colData = sample_info ,
                                     design = ~condition)
a=a[rowSums( counts(DESeq.ds) ) > 0]
DESeq.ds <- DESeq.ds[ rowSums( counts(DESeq.ds) ) > 0, ]
```

&nbsp;
&nbsp;
&nbsp;

## 4.1 RNA seq Analysis: Differential Expression (DE) Analysis
We will be generating plots using DESeq2 in R Studio.  
We are interested in for each gene, whether the differences in expression (counts) between groups is significant given the amount of variation observed within groups (replicates).

- Calculate the size factor/Dispersions/tests and add it to the data set  
  - sizeFactors: Using median-ratio method to estimate scale factors for each sample using stably expressed genes (genes with no zeros). For each sample, the scale factor is equal to the median of each gene's ratio to it's geometric mean. The geometric mean for each gene is calculated over all samples. This can be used to effectively normalize bulk RNA-seq data.  
  - Dispersion is the the variability between replicates. If you estimate dispersion = 0.19, then sqrt(dispersion) = BCV = 0.44. This means that the expression values vary up and down by 44% between replicates.  
  - Wald Test For The GLM Coefficients: To test the significance of coefficients in a Negative Binomial GLM, first we need to calculated sizeFactors (or normalizationFactors) and dispersion estimates.  
  - DESeq: Differential Expression Analysis Based On The Negative Binomial (A.K.A. Gamma-Poisson) Distribution. After the DESeq function returns a DESeqDataSet object, results tables (log2 fold changes and p-values) can be generated using the results function.  

```
DESeq.ds <- estimateSizeFactors(DESeq.ds)
DESeq.ds <- estimateDispersions(DESeq.ds)
DESeq.ds <- nbinomWaldTest(DESeq.ds)

## 3.5 deseq DE analysis restuls: DESeq.ds   
## You can run the code below to get SizeFactors, Dispersions, nbinomWaldTest at once.
DESeq.ds <- DESeq(DESeq.ds)
```

## 4.2 RNA seq Analysis: Down Stream Plots
### a. Normalized count data for down-stream plots
To generate a heatmap of RNA-seq results, we need a table of normalized counts.  
- Normalize count data for down-stream analysis like box-plot, volcano plots ...
counts() allows you to immediately retrieve the normalized read counts
```
DESeq.ds <- DESeq(DESeq.ds) # It can compute the sizeFactors Dispersions and tests automatically.
counts.sf_normalized <- counts(DESeq.ds , normalized = TRUE ) #normalization need sizeFactors
counts_not_normal <- counts(DESeq.ds)
names(a)=NULL
rownames(counts.sf_normalized )=a
rownames(counts_not_normal )=a
```
- Obtain Log tansformed values
```
# transform size - factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2( counts.sf_normalized + 1) 
log.not.norm.counts <- log2(counts_not_normal + 1)
```
- Obtain regularized log - transformed values
```
##  Note: blind = FALSE (use design information) / TRUE (do not use design information) ... 
DESeq.rlog <- rlog (DESeq.ds , blind = FALSE )  ## this step takes a while ... 
rlog.norm.counts <- assay (DESeq.rlog )
rownames(rlog.norm.counts )=a
```
### b. Heatmap
Heatmap is useful to for visualizing the expression of genes across the samples.
```
#plot2
 hm.mat_sig.gene <- subset(rlog.norm.counts, rownames(rlog.norm.counts) %in% <gene names>) 
#eg: c('Kdm6b','Kdm5b','Satb1','Tet1','Hdac2','Skp1a','Tdrkh','Cdk1','Aurka','Aurkb','Hmgb3','Ezh2','Uhrf1')
pdf("<Plot name>.pdf")
par(mar=c(1,1,1,1))
#Input the number of genes
aheatmap (hm.mat_sig.gene[,1:<numOfgenes>],
          Rowv = TRUE , Colv = T,
          distfun = "euclidean", hclustfun = "average",
          scale = "row") 
dev.off()
```
### c. Volcano plot
This tutorial [here](http://dputhier.github.io/jgb71e-polytech-bioinfo-app/practical/rna-seq_R/rnaseq_diff_Snf2.html#volcano-plot) is very helpful.

&nbsp;
&nbsp;
&nbsp;

## References
1. SRA: https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
2. Longleaf: https://its.unc.edu/research-computing/longleaf-cluster/
3. STAR: https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
4. Workflow: https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lectures/NGS_workflows.pdf
5. read count: https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/05_counting_reads.md
6. DE analysis: https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2/blob/master/lessons/DE_analysis.md
7. Volcano plot: http://dputhier.github.io/jgb71e-polytech-bioinfo-app/practical/rna-seq_R/rnaseq_diff_Snf2.html#volcano-plot
