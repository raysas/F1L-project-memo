# F1L Internship Emulator - personal memo

As introduced in [Dean Lee's linkedin newsletter](https://www.linkedin.com/pulse/introducing-f1l-internship-emulator-dean-lee-xckee/), the F1L Internship Emulator is an example of computational bio (intern) work in biotech/pharma. Following the recommendation of the author in [the first week's atricle](https://www.linkedin.com/pulse/week-1-f1l-internship-emulator-ksq-dean-lee-354ke/?trackingId=7%2BNj91jVSHm3z92bY3RYIg%3D%3D), this will be a documentation of my progress in the F1L project.

## Week 1: The KSQ

***The KSQ: Using available scRNA-seq data from cancer cell lines, how would you explore the use of the following FDA-approved antibody therapies in additional cancers?***

### Terms:  

* **scRNA-seq**: single-cell RNA sequencing, technique used to quantify gene expression levels of cells, where unlike bulk RNA-seq, the samples are _cells_, the feature are genes. The data is processed into a matrix of cells x genes. It allows for cell to cell comparison by highlighting the differences between a population of cells. This has won Nature method of the year in 2013.  
  * Workflow:
    * cell isolation
    * extraction of the genetic matrial + amplification
    * library prep
    * sequencing (NGS)
    * analysis (QC, alignment, quantification, normalization, clustering, etc.)
  * Uses:
    * Allow identification of cell types that are involved in a particular disease
    * The heterogenity of cells allows the identification of rare cell types 
    * Identification of expression patters of genes among cells -> gene co-regulated modules -> GRN between cells
    * PCA for cell subpopulations identification 
    * Cancer cell lines (?)
  * HCA: Human Cell Atlas, a global effort to map all cells in the human body. This would form a reference map - like the Human Genome Project did for the genome. Fully released biological networks for the lungs and brain are out on the [data portal](https://data.humancellatlas.org/), and efforts are still being done until the whole organism is mapped. Works are now on 18 biological networks. Each network is a group aiming to form the map for a specific organ/system/tissues as a step towards the ultimate goal of creating a human reference map. On the data portal, each project contains the raw sequencing data + associated files + metadata (cell type, seq details...).
  * Analysis: tools known for scRNA-seq analysis are: [Seurat R package](https://satijalab.org/seurat/), [Scanpy python toolkit](https://scanpy.readthedocs.io/en/stable/)
    * Preprocessing: as is the case in sequencing data in general and common with bulk RNA-seq we have teh following preprocessinf steps:
      * QC: quality control, to remove low quality cells. Known tools are FastQC and trimmomatic
      * Alignment: mapping the reads to the reference genome. e.g., STAR, HISAT2
      * Filtering (?)
      * Quantification: turn it into counts of reads per gene. e.g., featureCounts
    * Normalization: taking into consideration library size and sequencing depth
    * DEA: identify differential expression between cell groups - e.g., DESeq2
    * Clustering: grouping based on gene expression patterns, unsupervised machine learning can be performed using Seurat functions
    * Visualization: using dimentionality reduction techniques like PCA, t-SNE, UMAP
* **Cancer Cell Lines**:

### Data & code:

- [Smart-Seq2 dataset from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130)  
```bash 
# downloading the data - 65 mb
curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102130/suppl/GSE102130%5FK27Mproject.RSEM.vh20170621.txt.gz > data/Smart-Seq2.txt.gz
gunzip data/Smart-Seq2.txt.gz #makes it ~250mb
```


### References:  
* [Single Cell Sequencing in a Nutshell - TheScientist](https://www.the-scientist.com/single-cell-sequencing-in-a-nutshell-71048)    
* [Human Cell Atlas](https://www.humancellatlas.org/learn-more/#event-launch-of-the-human-cell-atlas)  
* [The Why and How of scRNA-Seq – A Guide for Beginners](https://www.parsebiosciences.com/blog/the-why-and-how-of-scrna-seq-a-guide-for-beginners/)  
* [Pre-processing of Single-Cell RNA Data - Galaxy training](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing/tutorial.html)  
* [Single-Cell RNA Sequencing (scRNA-seq) Data Analysis With R - medium](https://medium.com/@s.islambey/single-cell-rna-sequencing-scrna-seq-data-analysis-with-r-7e9fd4bc9e88)  
* [ANALYSIS OF SINGLE CELL RNA-SEQ DATA](https://broadinstitute.github.io/2019_scWorkshop/data-preprocessing.html)
* [Smart-Seq2 dataset from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130)