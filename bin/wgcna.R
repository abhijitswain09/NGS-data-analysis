#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: wgcna.R count_matrix.tsv samplesheet.csv outdir")
suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"
samples <- read.csv(args[2], stringsAsFactors = FALSE)
if (!all(c("sample", "condition") %in% names(samples))) stop("Samplesheet needs sample and condition columns")
if (nrow(samples) < 15) stop("WGCNA needs at least 15 samples; increase sample count or omit --run_wgcna")
if (anyDuplicated(samples$sample)) stop("Samplesheet sample IDs must be unique")
if (!all(samples$sample %in% colnames(counts))) stop("Samplesheet IDs must match count-matrix columns")
counts <- counts[, samples$sample, drop = FALSE]
keep <- rowSums(counts >= 10) >= 3
counts <- counts[keep, , drop = FALSE]
if (nrow(counts) < 100) stop("Fewer than 100 genes remain after low-count filtering")
condition <- factor(samples$condition)
coldata <- data.frame(condition = condition, row.names = samples$sample)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
vst_counts <- assay(vst(dds, blind = TRUE))
datExpr <- as.data.frame(t(vst_counts))
good <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[good$goodSamples, good$goodGenes, drop = FALSE]
if (nrow(datExpr) < 15) stop("Fewer than 15 valid samples remain for WGCNA")
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 0)
fit <- sft$fitIndices
eligible <- fit$Power[!is.na(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.80]
soft_power <- if (length(eligible)) eligible[1] else 12
net <- blockwiseModules(datExpr, power = soft_power, networkType = "unsigned",
  TOMType = "unsigned", minModuleSize = 20, mergeCutHeight = 0.25,
  numericLabels = FALSE, saveTOMs = FALSE, verbose = 2)
module_colors <- net$colors
names(module_colors) <- colnames(datExpr)
module_table <- data.frame(gene_id = names(module_colors), module = module_colors)
me <- net$MEs
trait <- model.matrix(~ condition - 1, data = data.frame(condition = condition))
rownames(trait) <- samples$sample
module_trait_cor <- cor(me, trait, use = "pairwise.complete.obs")
module_trait_p <- corPvalueStudent(module_trait_cor, nSamples = nrow(datExpr))
module_trait <- data.frame(module = rep(rownames(module_trait_cor), ncol(module_trait_cor)),
  trait = rep(colnames(module_trait_cor), each = nrow(module_trait_cor)),
  correlation = as.vector(module_trait_cor), pvalue = as.vector(module_trait_p))
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
write.csv(module_table, file.path(args[3], "module_membership.csv"), row.names = FALSE)
write.table(module_trait, file.path(args[3], "module_trait_correlations.tsv"),
  sep = "\\t", quote = FALSE, row.names = FALSE)
saveRDS(list(network = net, soft_threshold = sft, expression = datExpr, module_trait = module_trait),
  file.path(args[3], "wgcna_network.rds"))
pdf(file.path(args[3], "soft_threshold.pdf"), width = 9, height = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  xlab = "Soft threshold", ylab = "Scale-free topology fit (signed R^2)", type = "n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
  labels = sft$fitIndices[,1], col = "red")
abline(h = 0.80, col = "blue")
dev.off()
pdf(file.path(args[3], "module_dendrogram.pdf"), width = 12, height = 7)
plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
  "WGCNA modules", dendroLabels = FALSE, hang = 0.03)
dev.off()
