process EDGER {
    label 'process_medium'
    conda 'bioconda::bioconductor-edger'
    input:
    path counts
    path samplesheet
    path script
    output:
    path 'all_genes.tsv', emit: results
    path 'DEGs.tsv'
    path 'DEGs_up.tsv'
    path 'DEGs_down.tsv'
    path 'edger_plots.pdf', optional: true, emit: plots
    script:
    """
    FDR=$params.padj LFC=$params.lfc Rscript "$script" "$counts" "$samplesheet" . "$params.contrast_reference" "$params.contrast_target"
    """
    stub:
    """
    touch DEGs.tsv DEGs_up.tsv DEGs_down.tsv all_genes.tsv edger_plots.pdf
    """
}
