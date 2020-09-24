## Install the package
#BiocManager::install("DESeq2")

## Define woring environment
Path = getwd()
setwd(Path)
library(DESeq2)



## Read input data
cts <- as.matrix(read.csv("gene_count_matrix.csv",sep=",",row.names="gene_id"))
coldata <- read.csv("sample_info.csv",head=T,row.names = 1)
coldata$condition <- factor(coldata$condition)



## Load data into DESeq required format with experimental design
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Tube")
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)



## Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds)


## Output results
summary(res)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered),  file="Gene_level_differential_expression_analysis.csv")
