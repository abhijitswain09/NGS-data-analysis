process DESEQ2 {
    label 'process_medium'
    conda 'bioconda::bioconductor-deseq2 bioconda::bioconductor-apeglm conda-forge::r-ggplot2 conda-forge::r-pheatmap'
    input:
    path counts
    path samplesheet
    path script
    output:
    path 'deseq2/*', emit: results
    script:
    """
    PADJ=${params.padj} LFC=${params.lfc} Rscript ${script} ${counts} ${samplesheet} deseq2 ${params.contrast_reference} ${params.contrast_target}
    """
    stub:
    """mkdir -p deseq2 && touch deseq2/DEGs.tsv deseq2/all_genes.tsv"""
}
