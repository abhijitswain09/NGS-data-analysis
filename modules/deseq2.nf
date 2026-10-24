process DESEQ2 {
    label 'process_medium'
    conda 'bioconda::bioconductor-deseq2 bioconda::bioconductor-apeglm conda-forge::r-ggplot2 conda-forge::r-pheatmap'
    input:
    path counts
    path samplesheet
    path script
    output:
    path 'all_genes.tsv', emit: all_genes
    path 'DEGs.tsv', emit: degs
    path 'DEGs_up.tsv'
    path 'DEGs_down.tsv'
    script:
    """
    PADJ=$params.padj LFC=$params.lfc Rscript "$script" "$counts" "$samplesheet" . "$params.contrast_reference" "$params.contrast_target"
    """
    stub:
    """
    touch DEGs.tsv DEGs_up.tsv DEGs_down.tsv all_genes.tsv
    """
}
