process HISAT2_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::hisat2=2.2.1 bioconda::samtools=1.21'
    input:
    tuple val(meta), path(reads)
    path index
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "${meta.id}.hisat2.log", emit: log
    script:
    """
    hisat2 -p ${task.cpus} -x "${index}/genome" -1 "${reads[0]}" -2 "${reads[1]}" \
      --summary-file "${meta.id}.hisat2.log" | samtools sort -@ ${task.cpus} -o "${meta.id}.bam" -
    samtools index "${meta.id}.bam"
    """
    stub:
    """touch ${meta.id}.bam ${meta.id}.bam.bai ${meta.id}.hisat2.log"""
}
