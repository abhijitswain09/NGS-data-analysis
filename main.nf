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
include { GO_ANNOTATIONS } from './modules/go_annotations'
include { GO_ANALYSIS } from './modules/go_analysis'
include { WGCNA } from './modules/wgcna'

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

    if (!params.skip_go) {
        GO_ANNOTATIONS()
        GO_ANALYSIS(DESEQ2.out.all_genes, gtf, GO_ANNOTATIONS.out.gaf, GO_ANNOTATIONS.out.obo,
                    file("${projectDir}/bin/go_analysis.R"))
    }
    if (params.run_wgcna) {
        WGCNA(MERGE_COUNTS.out.matrix, samplesheet, file("${projectDir}/bin/wgcna.R"))
    }
}
