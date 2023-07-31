#install packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
#import library
library("limma")
library("edgeR")
#import the data
count_matrix <- as.matrix(read.table("/home/scbb/abhijit/fc0.counts.txt",header = TRUE, sep = "\t" ,row.names = "X"))
#check if it's import properly
head(count_matrix, 2)
#get the desire column from the dataset
x = count_matrix[, c(3,4,5,6)]
#set the condition
sample_info <- c("ctr", "ctr",  "trt", "trt")
#create DGEList data class for count and sample information
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
et <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = et, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]  
#summary
summary(decideTests(object = et, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="ct4vsct8_condition_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="ct4vsct8_up_condition_dge.csv")
write.csv(as.data.frame(down_regulated), file="ct4vsct8_down_condition_dge.csv")
x = count_matrix[, c(3,4,7,8)]
#set the condition
sample_info <- c("ctr", "ctr",  "trt", "trt")
#create DGEList data class for count and sample information
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
et <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = et, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]  
#summary
summary(decideTests(object = et, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="ct4vsft4_condition_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="ct4vsft4_up_condition_dge.csv")
write.csv(as.data.frame(down_regulated), file="ct4vsft4_down_condition_dge.csv")
x = count_matrix[, c(7,8,5,6)]
#set the condition
sample_info <- c("ctr", "ctr",  "trt", "trt")
#create DGEList data class for count and sample information
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
et <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = et, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]  
#summary
summary(decideTests(object = et, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="ft4vsct8_condition_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="ft4vsct8_up_condition_dge.csv")
write.csv(as.data.frame(down_regulated), file="ft4vsct8_down_condition_dge.csv")
x = count_matrix[, c(7,8,9,10)]
#set the condition
sample_info <- c("ctr", "ctr",  "trt", "trt")
#create DGEList data class for count and sample information
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
et <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = et, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]  
#summary
summary(decideTests(object = et, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="ft4vsft8_condition_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="ft4vsft8_up_condition_dge.csv")
write.csv(as.data.frame(down_regulated), file="ft4vsft8_down_condition_dge.csv")
x = count_matrix[, c(9,10,5,6)]
#set the condition
sample_info <- c("ctr", "ctr",  "trt", "trt")
#create DGEList data class for count and sample information
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
et <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = et, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]  
#summary
summary(decideTests(object = et, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="ft8vsct8_condition_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="ft8vsct8_up_condition_dge.csv")
write.csv(as.data.frame(down_regulated), file="ft8vsct8_down_condition_dge.csv")
x = count_matrix[, c(3,4,9,10)]
#set the condition
sample_info <- c("ctr", "ctr",  "trt", "trt")
#create DGEList data class for count and sample information
dge <- DGEList(counts = x, group = factor(sample_info))
#Filter out the genes with low counts
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
#Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
#Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
#Testing for differential gene expression
et <- exactTest(object = dge)
#extract the table with adjusted p values (FDR)
top_degs = topTags(object = et, n = "Inf")
filter <- (abs(top_degs$table$logFC)>=1)
DEG <- top_degs$table[filter,]  
#summary
summary(decideTests(object = et, lfc = 1))
#export the data
write.csv(as.data.frame(DEG), file="ct4vsft8_condition_dge.csv")
resSig <- subset(DEG, FDR < 0.1)
up_regulated <- subset(resSig, logFC > 0)
down_regulated <- subset(resSig, logFC < 0)
write.csv(as.data.frame(up_regulated), file="ct4vsft8_up_condition_dge.csv")
write.csv(as.data.frame(down_regulated), file="ct4vsft8_down_condition_dge.csv")


