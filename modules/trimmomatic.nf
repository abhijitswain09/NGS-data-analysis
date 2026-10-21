process TRIMMOMATIC {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::trimmomatic=0.39'
    input:
    tuple val(meta), path(reads)
    path adapters
    output:
    tuple val(meta), path("${meta.id}_R{1,2}.paired.fq.gz"), emit: reads
    path "${meta.id}.trimmomatic.log", emit: log
    script:
    """
    trimmomatic PE -threads ${task.cpus} -phred33 \
      ${reads[0]} ${reads[1]} \
      ${meta.id}_R1.paired.fq.gz ${meta.id}_R1.unpaired.fq.gz \
      ${meta.id}_R2.paired.fq.gz ${meta.id}_R2.unpaired.fq.gz \
      ILLUMINACLIP:${adapters}:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${params.min_len} \
      2> ${meta.id}.trimmomatic.log
    rm -f ${meta.id}_R1.unpaired.fq.gz ${meta.id}_R2.unpaired.fq.gz
    """
    stub:
    """touch ${meta.id}_R1.paired.fq.gz ${meta.id}_R2.paired.fq.gz ${meta.id}.trimmomatic.log"""
}
