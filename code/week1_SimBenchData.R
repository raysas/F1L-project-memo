if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SimBenchData")

library(SimBenchData)
??SimBenchData::SimBenchData
