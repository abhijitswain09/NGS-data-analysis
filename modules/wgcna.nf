process WGCNA {
    label 'process_high'
    conda 'bioconda::bioconductor-deseq2 bioconda::bioconductor-wgcna'
    input:
    path counts
    path samplesheet
    path script
    output:
    path 'wgcna/*', emit: results
    script:
    """
    Rscript "$script" "$counts" "$samplesheet" wgcna
    """
    stub:
    """
    mkdir -p wgcna
    touch wgcna/module_membership.csv wgcna/wgcna_network.rds
    """
}
