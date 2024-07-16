# trying Seurat package as recommended by sezer islambey
# -- reference: medium post: Single-Cell RNA Sequencing (scRNA-seq) Data Analysis With R
#    https://medium.com/@s.islambey/single-cell-rna-sequencing-scrna-seq-data-analysis-with-r-7e9fd4bc9e88

install.packages("Seurat")
library(Seurat) 

getwd()
setwd('..')
data <- read.table('data/Smart-Seq2.txt', header = TRUE, row.names = 1)

# -------- explore --------



# -------- tutorial --------


# Create the Seurat object
glioma <- CreateSeuratObject(counts = mydata, min.cells = 300, min.features = 200)

# the "counts" slot
GetAssayData(object = glioma, slot = "counts")

# check the class of the object
class(glioma)

# Violin plots - based on nGenes and nUMIS
VlnPlot(object = glioma, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# make the GenePlot - number of genes vs number of UMIS I
FeatureScatter(object = glioma, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# check the dimension of the data
dim(glioma@assays$RNA)

# a QC step - filter cells having less than 500 and more than 8000 genes
glioma <- subset(x = glioma, subset = nFeature_RNA > 500 & nFeature_RNA < 8000)

# check again the dimension of the data
dim(glioma@assays$RNA)

# scale the data
all.genes <- rownames(x = glioma)
glioma <- ScaleData(object = glioma, features = all.genes)

# glioma@assays$RNA@counts
# glioma@assays$RNA@data

# do some checks
expr1 <- GetAssay(glioma, slot="counts")[,"MGH66. P08.BØ5"]
# first, let's see the raw expression value of the gene A2M of the cell with name MGH66.P08.B05
expr1["A2M",] # should be equal to 18.78
log1p(expr1["A2M" ,] * 10000 / sum(expr1)) 

# Normalization
glioma <- NormalizeData(object = glioma, normalization.method = "LogNormalize", scale.factor = 10000)

# Now, the normalized data is stored in the slot counts
# check what is the normalized expression of the same gene in the same cell
expr2 <- GetAssay(glioma, slot="counts")[, "MGH66.P08.B05"]
expr2["A2M",] # should be equal to 0.1721429
# that is how it was computed!
# should give the same result

# glioma[, "MGH66.P08.B05"]

# Finding variable genes

# Expmean function: Calculate mean of logged values in non-log space (return answer in log-space)
# LogVMR function: Calculate the variance to mean ratio (VMR) in non-logspace (return answer in log-space)
glioma <- FindVariableFeatures(object = glioma, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = glioma), 10)

# VariableFeatures(glioma)

plot1 <- VariableFeaturePlot(object = glioma)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
LabelPoints(plot = plot1, points = top10, repel = TRUE) 

