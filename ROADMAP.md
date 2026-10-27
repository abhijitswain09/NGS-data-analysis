# Pipeline Roadmap

This roadmap prioritizes improvements to the paired-end RNA-seq workflow described in the [README](README.md). It is a sequence of proposed work, not a release calendar; priorities may be revisited as requirements evolve.

## Current foundation

The workflow uses Nextflow DSL2 for paired-end read QC, Trimmomatic trimming, HISAT2 alignment, HTSeq counting, merged count matrices, DESeq2 and edgeR differential expression, GO over-representation and GSEA, optional WGCNA, and MultiQC reporting. Nextflow also writes timeline, report, and trace files.

## Prioritized roadmap

| Priority | Initiative | Why it matters | Depends on |
|---|---|---|---|
| 1 | Early input and design validation | Finds malformed samplesheets, missing reads, invalid contrasts, and inconsistent matrices before costly computation. | None |
| 2 | Run provenance manifest | Makes parameters, software, references, and changing GO resources identifiable for each run. | None; benefits from the validated effective configuration in priority 1. |
| 3 | Multi-factor designs and explicit contrasts | Supports covariates such as batch and makes DESeq2 and edgeR comparisons consistent. | Priority 1, plus agreed design and contrast semantics. |
| 4 | Single-end read support | Broadens the supported input type across QC, trimming, and alignment. | Priority 1, plus an explicit policy for mixed read layouts. |

### Priority 1 — Input and design validation

Validate required fields, non-empty values, unique sample IDs, read-file availability, condition labels, contrast levels, count-matrix/sample consistency, and annotation expectations before expensive processes begin. Messages should point users to the specific sample, field, or file that needs correction.

**Ready to complete when:** invalid cases fail early with actionable messages, while the existing valid paired-end workflow is unchanged.

### Priority 2 — Run provenance

Write a machine-readable record of effective parameters, workflow and software versions, reference and annotation identifiers/checksums, and the exact GO annotation and ontology inputs used. In particular, identify the downloaded files behind the current GO resource URLs.

**Ready to complete when:** each run's outputs include a parseable record that identifies the analysis configuration and external resources used.

### Priority 3 — Multi-factor differential expression

Allow design covariates and explicit contrasts, and use the same documented design intent in DESeq2 and edgeR. Define accepted covariate types and validation for confounded or otherwise unidentifiable designs before exposing the interface.

**Ready to complete when:** documented multi-factor examples produce consistent, intended contrasts across both tools, with clear errors for unsupported designs.

### Priority 4 — Single-end reads

Support rows with one FASTQ input and route them through suitable trimming, FastQC, and alignment paths. Preserve paired-end behavior and state whether single-end and paired-end samples may coexist in one run.

**Ready to complete when:** single-end runs produce the expected outputs and existing paired-end tests remain successful; mixed-layout behavior is documented and validated.

## Roadmap principles

- Preserve the existing paired-end defaults and output expectations unless a deliberate compatibility change is documented.
- Prefer actionable validation and traceable outputs over silent assumptions.
- Add focused test fixtures for each new behavior and retain regression coverage for the existing test profile.
- Revisit priorities based on user needs; no target dates or release commitments are implied.
