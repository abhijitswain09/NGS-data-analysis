#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { TRIMMOMATIC } from './modules/trimmomatic'
include { HISAT2_BUILD } from './modules/hisat2_build'
include { HISAT2_ALIGN } from './modules/hisat2_align'
include { HTSEQ_COUNT } from './modules/htseq_count'
include { MERGE_COUNTS } from './modules/merge_counts'
include { DESEQ2 } from './modules/deseq2'
include { EDGER } from './modules/edger'
include { MULTIQC } from './modules/multiqc'
workflow {
    if (!params.input) error "Missing --input samplesheet CSV"
    if (!params.fasta) error "Missing --fasta genome FASTA"
    if (!params.gtf) error "Missing --gtf annotation GTF"
    samplesheet = file(params.input, checkIfExists: true)
    ch_reads = channel.fromPath(samplesheet).splitCsv(header: true, strip: true).map { row ->
        if (!row.sample || !row.fastq_1 || !row.fastq_2 || !row.condition)
            error "Samplesheet rows need sample, fastq_1, fastq_2, and condition"
        tuple([id: row.sample, condition: row.condition],
              [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)])
    }
    fasta = file(params.fasta, checkIfExists: true)
    gtf = file(params.gtf, checkIfExists: true)
    adapters = file(params.adapters, checkIfExists: true)
    FASTQC_RAW(ch_reads)
    TRIMMOMATIC(ch_reads, adapters)
    FASTQC_TRIMMED(TRIMMOMATIC.out.reads)
    HISAT2_BUILD(fasta)
    index = HISAT2_BUILD.out.index.first()
    HISAT2_ALIGN(TRIMMOMATIC.out.reads, index)
    HTSEQ_COUNT(HISAT2_ALIGN.out.bam, gtf)
    count_files = HTSEQ_COUNT.out.counts.map { meta, counts -> counts }.collect()
    MERGE_COUNTS(count_files)
    DESEQ2(MERGE_COUNTS.out.matrix, samplesheet, file("${projectDir}/bin/deseq2.R"))
    EDGER(MERGE_COUNTS.out.matrix, samplesheet, file("${projectDir}/bin/edger.R"))
    qc = FASTQC_RAW.out.reports.mix(FASTQC_TRIMMED.out.reports)
        .mix(TRIMMOMATIC.out.log).mix(HISAT2_ALIGN.out.log).collect()
    MULTIQC(qc)
}
