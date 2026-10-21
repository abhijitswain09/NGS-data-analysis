#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) stop("Usage: deseq2.R counts.tsv samplesheet.csv outdir reference target")
counts <- as.matrix(read.delim(args[1], row.names=1, check.names=FALSE))
ss <- read.csv(args[2], stringsAsFactors=FALSE)
ref <- args[4]; target <- args[5]
ss <- ss[ss$condition %in% c(ref,target), c("sample","condition")]
if (!all(c(ref,target) %in% ss$condition)) stop("Both requested conditions must be present in samplesheet")
if (!all(ss$sample %in% colnames(counts))) stop("Samplesheet sample IDs must match count-matrix columns")
counts <- counts[,ss$sample,drop=FALSE]
ss$condition <- relevel(factor(ss$condition), ref=ref)
rownames(ss) <- ss$sample
dds <- DESeqDataSetFromMatrix(countData=counts, colData=ss, design=~condition)
dds <- dds[rowSums(counts(dds) >= 10) >= min(table(ss$condition)),]
dds <- DESeq(dds)
coef <- resultsNames(dds)[2]
res <- lfcShrink(dds, coef=coef, res=results(dds, name=coef), type="apeglm")
tab <- data.frame(gene_id=rownames(res), as.data.frame(res))
padj_cut <- as.numeric(Sys.getenv("PADJ","0.05"))
lfc_cut <- as.numeric(Sys.getenv("LFC","1"))
sig <- subset(tab, !is.na(padj) & padj < padj_cut & abs(log2FoldChange) > lfc_cut)
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
write.table(tab, file.path(args[3],"all_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sig, file.path(args[3],"DEGs.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(subset(sig,log2FoldChange>0), file.path(args[3],"DEGs_up.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(subset(sig,log2FoldChange<0), file.path(args[3],"DEGs_down.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
