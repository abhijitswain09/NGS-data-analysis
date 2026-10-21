process HTSEQ_COUNT {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::htseq=2.0.5'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    output:
    tuple val(meta), path("${meta.id}.counts.txt"), emit: counts
    script:
    """
    htseq-count -f bam -r pos -s ${params.strandedness} -t exon -i ${params.gtf_attribute} ${bam} ${gtf} > ${meta.id}.counts.txt
    """
    stub:
    """printf 'GENE_A\\t10\\nGENE_B\\t20\\n__no_feature\\t5\\n' > ${meta.id}.counts.txt"""
}
