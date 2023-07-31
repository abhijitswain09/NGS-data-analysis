library(limma)
library(DESeq2)
counts = read.table("/home/scbb/abhijit/genes_counts.txt", sep="\t", header=T,row.names=c(1))
x = counts[, c(5,6,7,8,14,15)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_1vsbtp_2_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_1vsbtp_2_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
counts = read.table("/home/scbb/abhijit/genes_counts.txt", sep="\t", header=T,row.names=c(1))
x = counts[, c(8,14,15,2,3,4)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_1vsbtp_3_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_1vsbtp_3_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(8,14,15,1,12,13)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_1vsbtp_4_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_1vsbtp_4_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(8,14,15,9,10,11)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_1vsbtp_5_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_1vsbtp_5_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(5,6,7,2,3,4)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_2vsbtp_3_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_2vsbtp_3_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(5,6,7,1,12,13)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_2vsbtp_4_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_2vsbtp_4_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(5,6,7,9,10,11)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_2vsbtp_5_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_2vsbtp_5_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(2,3,4,1,12,13)]
condition<-c("C","C","C","T","T","T")

coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_3vsbtp_4_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_3vsbtp_4_down_regulated.csv", quote=F, sep= "\t", col.names = NA)
x = counts[, c(2,3,4,9,10,11)]
condition<-c("C","C","C","T","T","T")
coldata <- data.frame(row.names=colnames(x), condition)
#The next step is to generate an object, called “dds”, containing the counts and conditions to be compared as entry for DESeq2
dds <- DESeqDataSetFromMatrix( countData=x[rowSums
(x)>0,], colData=coldata, design=~condition)
#Run the DESeq2 function to start the statistical analysis. It retrieves a normalized expression level in logarithmic base 2 scale that will be stored in “results” object
results <- results(DESeq(dds))
#Filter the data in “results” with the desired thresholds of fold change level and adjusted p-values. In this case fold change is log2FC >1 with adjusted p-values<0. Store the data in “filter” object
res <- na.exclude(as.data.frame(results))
#for counting up-regulated &down-regulated genes
resSig <- subset(res, padj < 0.1)
up_regulated <- subset(resSig, log2FoldChange > 0)
down_regulated <- subset(resSig, log2FoldChange < 0)
#Export the files
write.table(up_regulated,"btp_3vsbtp_5_up_regulated.csv", quote=F, sep= "\t", col.names = NA)
write.table(down_regulated,"btp_3vsbtp_5_down_regulated.csv", quote=F, sep= "\t", col.names = NA)




