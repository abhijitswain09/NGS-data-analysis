process HISAT2_BUILD {
    tag "${fasta.simpleName}"
    label 'process_high'
    conda 'bioconda::hisat2=2.2.1'
    input:
    path fasta
    output:
    path 'hisat2_index', emit: index
    script:
    """
    mkdir hisat2_index
    hisat2-build -p ${task.cpus} "${fasta}" hisat2_index/genome
    """
    stub:
    """mkdir hisat2_index && touch hisat2_index/genome.1.ht2"""
}
