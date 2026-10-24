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
    required_columns = ['sample', 'fastq_1', 'fastq_2', 'condition']
    if (!(params.strandedness in ['no', 'yes', 'reverse']))
        error "--strandedness must be one of: no, yes, reverse"
    if (!(params.gtf_attribute ==~ /[A-Za-z_][A-Za-z0-9_]*/))
        error "--gtf_attribute must be a valid GTF attribute name"
    if (!(params.contrast_reference ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) ||
        !(params.contrast_target ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) ||
        params.contrast_reference == params.contrast_target)
        error "Contrast labels must be distinct and contain only letters, numbers, dot, underscore, or hyphen"
    try {
        if (params.min_len.toInteger() < 1) error "--min_len must be a positive integer"
        if (!(params.padj.toBigDecimal() > 0 && params.padj.toBigDecimal() <= 1))
            error "--padj must be greater than 0 and at most 1"
        if (params.lfc.toBigDecimal() < 0) error "--lfc must be non-negative"
    } catch (NumberFormatException e) {
        error "--min_len, --padj, and --lfc must be numeric (min_len must be an integer)"
    }
    rows_ch = channel.fromPath(samplesheet).splitCsv(header: true, strip: true).collect()
    ch_reads = rows_ch.flatMap { rows ->
        if (!rows || !required_columns.every { rows[0].containsKey(it) })
            error "Samplesheet must contain columns: ${required_columns.join(', ')}"
        if (rows.any { row -> required_columns.any { !row[it] } })
            error "Every samplesheet row needs non-empty sample, fastq_1, fastq_2, and condition values"
        if (rows.any { !(it.sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) || !(it.condition ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) })
            error "Sample IDs and condition labels may contain only letters, numbers, dot, underscore, or hyphen, and must start with a letter or number"
        sample_ids = rows.collect { it.sample }
        if (sample_ids.toSet().size() != sample_ids.size()) error "Samplesheet sample IDs must be unique"
        conditions = rows.collect { it.condition }.toSet()
        if (!(params.contrast_reference in conditions) || !(params.contrast_target in conditions))
            error "Both contrast conditions must occur in the samplesheet: ${params.contrast_reference}, ${params.contrast_target}"
        ref_count = rows.count { it.condition == params.contrast_reference }
        target_count = rows.count { it.condition == params.contrast_target }
        if (ref_count < 2 || target_count < 2)
            error "At least two samples are required for each contrast condition"
        rows.collect { row ->
            tuple([id: row.sample, condition: row.condition],
                  [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)])
        }
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
