if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SimBenchData")
BiocManager::install("ExperimentHub")

library(SimBenchData)
??SimBenchData::SimBenchData

# -- following vignette
library(ExperimentHub)
eh <- ExperimentHub()
alldata <- query(eh, "SimBenchData")
alldata 

data_1 <- alldata[["EH5384"]]  

metadata <- showMetaData()
View(metadata)
