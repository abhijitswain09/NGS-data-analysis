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
    path 'go_analysis/*', emit: results
    script:
    """
    Rscript "$script" "$all_genes" "$gtf" "$gaf" "$obo" go_analysis "$params.padj" "$params.lfc"
    """
    stub:
    """
    mkdir -p go_analysis
    touch go_analysis/GO_ORA.tsv go_analysis/GO_GSEA.tsv
    """
}
