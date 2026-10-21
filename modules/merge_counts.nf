process MERGE_COUNTS {
    label 'process_low'
    conda 'conda-forge::python=3.12'
    input:
    path counts
    output:
    path 'count_matrix.tsv', emit: matrix
    script:
    """python ${projectDir}/bin/merge_counts.py ${counts} > count_matrix.tsv"""
    stub:
    """python ${projectDir}/bin/merge_counts.py ${counts} > count_matrix.tsv"""
}
