# install required packages
BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)


# load count matrix
dataFilt <- get(load(file = 'dataFilt.RData'))

# load clinical data
Clinical_data <- read.csv("trait_data.csv")
Clinical_data <- data.frame(Clinical_data, row.names = T)


all(rownames(Clinical_data) %in% colnames(dataFilt))
all(rownames(Clinical_data) == colnames(dataFilt))


dds <- DESeqDataSetFromMatrix(countData = dataFilt,
                              colData = primary,
                              design= ~braf_status)


# remove gene symbols with low read counts
Keep <- rowSums(counts(dds) >= 10)
dds <- dds[Keep,]

# set reference label
dds$braf_status <- relevel(dds$braf_status, ref="WT")


# Perform DESeq for our expression values
dds <- DESeq(dds)    
resultsNames(dds)
res <- results(dds)
summary(res)

# for p-value less than 0.05 and threshold 1 (>1 up-regulated and <1 down-regulated)
res.05 <- data.frame(results(dds,alpha = 0.05,lfcThreshold=1))
summary(res.05)

# save file
write.csv(res.05, file = "res.05.csv", row.names = T)
