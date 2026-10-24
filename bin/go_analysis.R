#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) stop("Usage: go_analysis.R all_genes.tsv annotation.gtf annotations.gaf.gz go-basic.obo outdir padj lfc")

results <- read.delim(args[1], check.names = FALSE)
gtf_lines <- readLines(args[2], warn = FALSE)
gtf_lines <- gtf_lines[grepl("\\t(exon|gene)\\t", gtf_lines)]
get_attribute <- function(x, key) {
  hit <- regexec(paste0(key, ' "([^"]+)"'), x)
  vapply(regmatches(x, hit), function(m) if (length(m) >= 2) m[2] else NA_character_, character(1))
}
gene_map <- unique(data.frame(
  gene_id = get_attribute(gtf_lines, "gene_id"),
  gene_name = get_attribute(gtf_lines, "gene_name"),
  stringsAsFactors = FALSE
))
gene_map <- gene_map[!is.na(gene_map$gene_id) & !is.na(gene_map$gene_name), ]
if (!nrow(gene_map)) stop("GTF must include gene_id and gene_name attributes to map GO annotations")

obo <- readLines(args[4], warn = FALSE)
blocks <- split(obo, cumsum(grepl("^\\[Term\\]$", obo)))
terms <- lapply(blocks, function(x) {
  id <- sub("^id: ", "", grep("^id: GO:", x, value = TRUE)[1])
  name <- sub("^name: ", "", grep("^name: ", x, value = TRUE)[1])
  ns <- sub("^namespace: ", "", grep("^namespace: ", x, value = TRUE)[1])
  obsolete <- any(grepl("^is_obsolete: true", x))
  if (is.na(id) || is.na(name) || is.na(ns) || obsolete) return(NULL)
  data.frame(go_id = id, term = name, namespace = ns)
})
term_table <- do.call(rbind, terms)
if (is.null(term_table) || !nrow(term_table)) stop("Could not parse GO terms from the ontology file")

gaf <- read.delim(gzfile(args[3]), header = FALSE, sep = "\t", comment.char = "!", quote = "", fill = TRUE, stringsAsFactors = FALSE)
if (ncol(gaf) < 9) stop("GO annotation file is not a valid GAF with at least 9 columns")
gaf <- gaf[!grepl("(^|\\|)NOT(\\||$)", gaf[[4]]) & gaf[[5]] %in% term_table$go_id, c(3, 5)]
names(gaf) <- c("gene_name", "go_id")
annotations <- unique(merge(gaf, gene_map, by = "gene_name"))
if (!nrow(annotations)) stop("No GO symbols matched GTF gene_name values; supply matching GO annotations or a compatible GTF")

results$gene_id <- as.character(results$gene_id)
map <- gene_map[match(results$gene_id, gene_map$gene_id), c("gene_id", "gene_name")]
results$gene_name <- map$gene_name
results <- results[!is.na(results$gene_name), ]
if (!"stat" %in% names(results)) {
  if (!all(c("log2FoldChange", "pvalue") %in% names(results))) stop("DESeq2 results need stat or log2FoldChange and pvalue columns")
  results$stat <- sign(results$log2FoldChange) * -log10(pmax(results$pvalue, .Machine$double.xmin))
}
results <- results[is.finite(results$stat), ]
results <- results[order(abs(results$stat), decreasing = TRUE), ]
results <- results[!duplicated(results$gene_name), ]
ranked <- results$stat
names(ranked) <- results$gene_name
ranked <- sort(ranked, decreasing = TRUE)

annotations <- annotations[annotations$gene_name %in% unique(results$gene_name), ]
background <- unique(annotations$gene_name)
degs <- unique(results$gene_name[
  !is.na(results$padj) & results$padj <= as.numeric(args[6]) &
    abs(results$log2FoldChange) >= as.numeric(args[7]) &
    results$gene_name %in% background
])
pathways <- lapply(split(annotations$gene_name, annotations$go_id), unique)
pathways <- pathways[lengths(pathways) >= 10 & lengths(pathways) <= 500]
if (!length(pathways)) warning("No GO sets of size 10-500 overlap the tested genes; writing empty enrichment tables")

ora_rows <- lapply(names(pathways), function(go_id) {
  term_genes <- intersect(pathways[[go_id]], background)
  overlap <- intersect(term_genes, degs)
  if (length(overlap) == 0 || length(degs) == 0) return(NULL)
  a <- length(overlap); b <- length(degs) - a
  c <- length(term_genes) - a; d <- length(background) - length(degs) - c
  p <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")$p.value
  data.frame(go_id = go_id, overlap = length(overlap), set_size = length(term_genes),
             genes = paste(overlap, collapse = ","), pvalue = p)
})
ora <- do.call(rbind, ora_rows)
if (is.null(ora)) ora <- data.frame(go_id = character(), overlap = integer(), set_size = integer(), genes = character(), pvalue = numeric())
ora$padj <- p.adjust(ora$pvalue, method = "BH")
ora <- merge(ora, term_table, by = "go_id", all.x = TRUE)
ora <- ora[order(ora$padj), ]

suppressPackageStartupMessages(library(fgsea))
gsea <- if (length(pathways)) {
  as.data.frame(fgsea(pathways = pathways, stats = ranked, minSize = 10, maxSize = 500, scoreType = "std"))
} else {
  data.frame(go_id = character(), pathway = character(), pval = numeric(), padj = numeric(), NES = numeric(), size = integer())
}
if (nrow(gsea)) {
  gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ",", character(1))
  names(gsea)[names(gsea) == "pathway"] <- "go_id"
  gsea <- merge(gsea, term_table, by = "go_id", all.x = TRUE)
  gsea <- gsea[order(gsea$padj), ]
} else {
  gsea <- data.frame(go_id = character(), term = character(), NES = numeric(), pval = numeric(), padj = numeric())
}

dir.create(args[5], recursive = TRUE, showWarnings = FALSE)
write.table(ora, file.path(args[5], "GO_ORA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea, file.path(args[5], "GO_GSEA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(ora)) {
  pdf(file.path(args[5], "GO_ORA_top_terms.pdf"), width = 10, height = 7)
  top <- head(ora[!is.na(ora$padj), ], 20)
  if (nrow(top)) {
    top$term <- factor(top$term, levels = rev(top$term))
    print(ggplot2::ggplot(top, ggplot2::aes(term, -log10(padj))) +
      ggplot2::geom_col(fill = "steelblue") + ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "-log10 adjusted p-value", title = "GO over-representation") +
      ggplot2::theme_minimal())
  }
  dev.off()
}
if (nrow(gsea)) {
  pdf(file.path(args[5], "GO_GSEA_top_terms.pdf"), width = 10, height = 7)
  top <- head(gsea[!is.na(gsea$padj), ], 20)
  if (nrow(top)) {
    top$term <- factor(top$term, levels = top$term[order(top$NES)])
    print(ggplot2::ggplot(top, ggplot2::aes(term, NES, fill = NES > 0)) +
      ggplot2::geom_col() + ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Normalized enrichment score", title = "GO gene set enrichment") +
      ggplot2::theme_minimal() + ggplot2::guides(fill = "none"))
  }
  dev.off()
}
