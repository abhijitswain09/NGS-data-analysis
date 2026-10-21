nextflow.enable.dsl=2

params.samplesheet = "${projectDir}/samplesheet.csv"
params.reference = null
params.annotation = null
params.outdir = "results"

process TRIM_READS {
    tag "$sample"
    conda "${projectDir}/environment.yml"
    cpus 4
    input:
    tuple val(sample), path(read1), path(read2)
    output:
    tuple val(sample), path("${sample}_R1.paired.fq.gz"), path("${sample}_R2.paired.fq.gz"), emit: reads
    script:
    """
    trimmomatic PE -threads ${task.cpus} -phred33 "${read1}" "${read2}" \
      ${sample}_R1.paired.fq.gz ${sample}_R1.unpaired.fq.gz \
      ${sample}_R2.paired.fq.gz ${sample}_R2.unpaired.fq.gz \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process ALIGN_AND_COUNT {
    tag "$sample"
    conda "${projectDir}/environment.yml"
    cpus 8
    input:
    tuple val(sample), path(read1), path(read2)
    path reference
    path annotation
    output:
    path "${sample}.bam", emit: bam
    path "${sample}.counts.tsv", emit: counts
    script:
    """
    hisat2-build -p ${task.cpus} "${reference}" genome_index
    hisat2 -p ${task.cpus} -x genome_index -1 "${read1}" -2 "${read2}" 2> ${sample}.hisat2.log | samtools sort -@ ${task.cpus} -o ${sample}.bam
    htseq-count -f bam -r pos -s no -t exon -i gene_id ${sample}.bam "${annotation}" > ${sample}.counts.tsv
    """
}

workflow {
    if (!params.reference || !params.annotation) error "Provide --reference FASTA and --annotation GTF"
    channel.fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header: true)
      .map { row -> tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)) }
      .set { samples }
    ALIGN_AND_COUNT(TRIM_READS(samples), file(params.reference, checkIfExists: true), file(params.annotation, checkIfExists: true))
}