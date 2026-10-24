process WGCNA {
    label 'process_high'
    conda 'bioconda::bioconductor-deseq2 bioconda::bioconductor-wgcna'
    input:
    path counts
    path samplesheet
    path script
    output:
    path 'module_membership.csv'
    path 'module_trait_correlations.tsv'
    path 'wgcna_network.rds'
    path 'soft_threshold.pdf', optional: true
    path 'module_dendrogram.pdf', optional: true
    script:
    """
    Rscript "$script" "$counts" "$samplesheet" .
    """
    stub:
    """
    touch module_membership.csv module_trait_correlations.tsv wgcna_network.rds
    """
}
