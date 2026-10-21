nextflow.enable.dsl=2
params.samplesheet = "${projectDir}/samplesheet.csv"
params.reference = null
params.annotation = null
params.adapters = null
params.outdir = 'results'
process BUILD_INDEX {
  conda "${projectDir}/environment.yml"
  input: path reference
  output: path 'genome_index.*', emit: index
  script: "hisat2-build -p ${task.cpus} '${reference}' genome_index"
}
process TRIM_READS {
  tag "$sample"
  conda "${projectDir}/environment.yml"
  cpus 4
  input: tuple val(sample), path(read1), path(read2), path(adapters)
  output: tuple val(sample), path("${sample}_R1.paired.fq.gz"), path("${sample}_R2.paired.fq.gz"), emit: reads
  script:
  """
  trimmomatic PE -threads ${task.cpus} -phred33 "${read1}" "${read2}" ${sample}_R1.paired.fq.gz ${sample}_R1.unpaired.fq.gz ${sample}_R2.paired.fq.gz ${sample}_R2.unpaired.fq.gz ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  """
}
process ALIGN_AND_COUNT {
  tag "$sample"
  conda "${projectDir}/environment.yml"
  cpus 8
  input: tuple val(sample), path(read1), path(read2)
  path index_files
  path annotation
  output: path "${sample}.bam", emit: bam
  output: path "${sample}.counts.tsv", emit: counts
  script:
  """
  hisat2 -p ${task.cpus} -x genome_index -1 "${read1}" -2 "${read2}" 2> ${sample}.hisat2.log | samtools sort -@ ${task.cpus} -o ${sample}.bam
  htseq-count -f bam -r pos -s no -t exon -i gene_id ${sample}.bam "${annotation}" > ${sample}.counts.tsv
  """
}
workflow {
  if (!params.reference || !params.annotation || !params.adapters) error 'Provide --reference FASTA, --annotation GTF, and --adapters FASTA'
  channel.fromPath(params.samplesheet, checkIfExists: true).splitCsv(header: true).map { row -> tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true), file(params.adapters, checkIfExists: true)) }.set { samples }
  index = BUILD_INDEX(file(params.reference, checkIfExists: true))
  trimmed = TRIM_READS(samples)
  ALIGN_AND_COUNT(trimmed.reads, index.index, file(params.annotation, checkIfExists: true))
}