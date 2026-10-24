process GO_ANALYSIS {
    label 'process_medium'
    conda 'bioconda::bioconductor-fgsea conda-forge::r-ggplot2'
    input:
    path all_genes
    path gtf
    path gaf
    path obo
    path script
    output:
    path 'GO_ORA.tsv', emit: ora
    path 'GO_GSEA.tsv', emit: gsea
    path 'GO_ORA_top_terms.pdf', optional: true
    path 'GO_GSEA_top_terms.pdf', optional: true
    script:
    """
    Rscript "$script" "$all_genes" "$gtf" "$gaf" "$obo" . "$params.padj" "$params.lfc"
    """
    stub:
    """
    touch GO_ORA.tsv GO_GSEA.tsv
    """
}
