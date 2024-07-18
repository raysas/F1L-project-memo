# F1L Internship Emulator - personal memo 📝

As introduced in [Dean Lee's linkedin newsletter](https://www.linkedin.com/pulse/introducing-f1l-internship-emulator-dean-lee-xckee/), the F1L Internship Emulator is an example of computational bio (intern) work in biotech/pharma. Following the recommendation of the author in [the first week's atricle](https://www.linkedin.com/pulse/week-1-f1l-internship-emulator-ksq-dean-lee-354ke/?trackingId=7%2BNj91jVSHm3z92bY3RYIg%3D%3D), this will be a documentation of my progress in the F1L project.

## Week 1: The Key Scientific Question 🔎

***The KSQ: Using available scRNA-seq data from cancer cell lines, how would you explore the use of the following FDA-approved antibody therapies in additional cancers?***

### Terms:  

**scRNA-seq**:  
single-cell RNA sequencing, technique used to quantify gene expression levels of cells, where unlike bulk RNA-seq, the samples are _cells_, the feature are genes. The data is processed into a matrix of cells x genes. It allows for cell to cell comparison by highlighting the differences between a population of cells. This has won Nature method of the year in 2013.  
* Workflow:
  * cell isolation
  * extraction of the genetic matrial + amplification
  * library prep
  * sequencing (NGS)
  * analysis (QC, alignment, quantification, normalization, clustering, etc.)
Uses:
* Allow identification of cell types that are involved in a particular disease
* The heterogenity of cells allows the identification of rare cell types 
* Identification of expression patters of genes among cells -> gene co-regulated modules -> GRN between cells
* PCA for cell subpopulations identification 
* Cancer cell lines (?)  
  
HCA: Human Cell Atlas, a global effort to map all cells in the human body. This would form a reference map - like the Human Genome Project did for the genome. Fully released biological networks for the lungs and brain are out on the [data portal](https://data.humancellatlas.org/), and efforts are still being done until the whole organism is mapped. Works are now on 18 biological networks. Each network is a group aiming to form the map for a specific organ/system/tissues as a step towards the ultimate goal of creating a human reference map. On the data portal, each project contains the raw sequencing data + associated files + metadata (cell type, seq details...).

Analysis: tools known for scRNA-seq analysis are: [Seurat R package](https://satijalab.org/seurat/), [Scanpy python toolkit](https://scanpy.readthedocs.io/en/stable/)
* Preprocessing: as is the case in sequencing data in general and common with bulk RNA-seq we have teh following preprocessinf steps:
* QC: quality control, to remove low quality cells. Known tools are FastQC and trimmomatic
    * Alignment: mapping the reads to the reference genome. e.g., STAR, HISAT2
    * Filtering (?)
    * Quantification: turn it into counts of reads per gene. e.g., featureCounts
  * Normalization: taking into consideration library size and sequencing depth
  * DEA: identify differential expression between cell groups - e.g., DESeq2
  * Clustering: grouping based on gene expression patterns, unsupervised machine learning can be performed using Seurat functions
  * Visualization: using dimentionality reduction techniques like PCA, t-SNE, UMAP

**Cancer cells**: 
When things go wrong, normal cells transform into cancer cells that are able to proliferate rapidly. Usually it can be caused due to mutations in:  
* Tumor suppressor genes: genes that prevent cell division, allow cell death, repair DNA... When there is a mutation, the products of these genes (proteins) will lose their functionality => no more brakes on cell division.
The two hit hypothesis: a cell needs two mutations in a tumor suppressor gene to become cancerous. Why two? for both alleles to be turned off (one mutation is like having one out of 2 brakes broken).
_e.g., BRCA1, BRCA2. These are eventually translated into proteins that repair DNA_   
* Oncogenes: to define them would need to know about a _proto-oncogene_ which is a gene that is involved in cell growth and division. When a proto-oncogene mutates, it becomes an oncogene that can cause cells to grow and divide uncontrollably (like being turned on all the time). This will lead to cancer. First oncogene discovered is in Sarcoma (a type of cancer) and is called _SRC_, when mutated it gets over expressed.  

It's worth noting that not all mutations have effects, we have: silent mutations, missense mutations, nonsense mutations, frameshift mutations, etc.  
Common cancer characteristics that allow spreading in different cancer types are called cancer hallmarks (~10), which work for:  
* Mutation of its own genetic material: can change cellular makeup to resist and survive - change DNA (property that normal cells dont really have) (this one is called genome instability and mutation)
* Cell division: allow proliferation by being unchecked, sustain signals - cancer cells keep proliferating (3 hallmarks)
  * Evasion of growth suppressors: escape and keep on dividing 
  * Enable replicative immortality: keep on dividing without dying (immortal)
  * Dyregulating cellular energetics: change the way they get energy (aerobic glycolysis)
* Spreading disease: (2 hallmarks)
  * Inducing angiogenesis: make blood vessels grow towards them to get nutrients. Blood vessels spread all around body for the purpose of spreading nutrients, moving cells... Cancer cells use them to spread
  * Activation of invasion: invasion is moving from one tissue to another, usually normal cells stay in their place - not tha case in cancer cells. They can move to other tissues and start growing there.
* Escape normal body defense mechanisms:
  * Resistsing cell death: resist the mechanisms to grow
  * Avoid immune destruction: can not be recognized by immune systems (immune system is supposed to recognize and kill cancer cells)
* Tumor promoting inflamation: (inflamation is body repsonse to disease), cancer cells can use this inflamation to grow, they use the inflamation and take advantage of it  

Research aims to find ways to target these hallmarks to stop cancer growth.

Neoplasia is the abnormal growth of cells, it can be benign or malignant. In short benign is not cancerous and malignant is cancerous. Normal tissue have a particular order (epithelial => basement => blood/lumph vessels).  
* _Dysplasia_: losing this order and becoming disorganized.
* _Preinvasive carcinoma_: when cells start to invade and take the whole layer (p.s. carcinoma = cancer of the epithelium).  
* _Invasive carcinoma_: when cells start to bottom layer => basement membrane, get to vessels.  
* _Metatstatic carcinoma (metastasis)_: when cells start to spread to other tissues.

Now metastasis (last step of its progression) is the most dangerous part of cancer, it's when cells start to spread to other tissues. It's the most dangerous because it's hard to treat, detect, remove. Also found in important and vital organs. They can metastsize in different ways:  
- lymphatic spread: spread through lymphatic system
- hematogenous spread: spread through blood vessels

They can spread to cavities, bones, brain, liver, lungs, lymph nodes, peritoneum, pleura, skin, soft tissues, etc. difficult to know where its going through these routes. There are several steps to be achieved for a cell to metastasize. it can be local (same organ) or distant (other organ).  

_What's metastatic tropism?_  
It's the ability of cancer cells to spread to specific organs. It's not random, it's specific. It's not clear why it happens, but it's known that some cancers have a preference to spread to specific organs.


**Cancer Cell Lines**:

**FDA-approved antibody therapies**:


### Data & code:

- [Smart-Seq2 dataset from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130)  
```bash 
# downloading the data - 65 mb
curl https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102130/suppl/GSE102130%5FK27Mproject.RSEM.vh20170621.txt.gz > data/Smart-Seq2.txt.gz
gunzip data/Smart-Seq2.txt.gz #makes it ~250mb
```

- [R code for Seurat exploration](code/week1_Seurat_tutorial.R)

### References:  
* [Single Cell Sequencing in a Nutshell - TheScientist](https://www.the-scientist.com/single-cell-sequencing-in-a-nutshell-71048)    
* [Human Cell Atlas](https://www.humancellatlas.org/learn-more/#event-launch-of-the-human-cell-atlas)  
* [The Why and How of scRNA-Seq – A Guide for Beginners](https://www.parsebiosciences.com/blog/the-why-and-how-of-scrna-seq-a-guide-for-beginners/)  
* [Pre-processing of Single-Cell RNA Data - Galaxy training](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing/tutorial.html)  
* [Single-Cell RNA Sequencing (scRNA-seq) Data Analysis With R - medium](https://medium.com/@s.islambey/single-cell-rna-sequencing-scrna-seq-data-analysis-with-r-7e9fd4bc9e88)  
* [ANALYSIS OF SINGLE CELL RNA-SEQ DATA](https://broadinstitute.github.io/2019_scWorkshop/data-preprocessing.html)
* [Smart-Seq2 dataset from GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102130)  
* [Cancer Cell Biology for Beginners - youtube playlist](https://www.youtube.com/playlist?list=PLNrGtdJ6nPMh8JSe1JhYmBDrQ52whsZps)