# RNA-seq Workflow Improvement Plan

This plan describes a proposed sequence for improving the Nextflow RNA-seq workflow. It is documentation for future work; it does not mean these changes have been implemented or scheduled.

## Goals

- Catch invalid inputs and unsupported analysis designs before expensive compute starts.
- Make each run's configuration, software, references, and external annotation resources traceable.
- Support common experimental designs beyond a single two-condition comparison.
- Extend the workflow to single-end data without regressing paired-end analyses.

## Scope

The plan covers the samplesheet and analysis interface, provenance outputs, differential-expression design configuration, and read-layout support. Existing paired-end behavior should remain the compatibility baseline. Each feature should include user documentation and focused automated tests.

Out of scope for this plan are new analysis methods, changes to biological defaults or statistical thresholds without separate justification, and release-date commitments.

## Implementation phases

### Phase 1 — Validate inputs and analysis design

1. Define and document samplesheet requirements, including required columns, non-empty values, and unique sample identifiers.
2. Validate read paths and paired-read expectations, condition labels, selected contrast levels, count-matrix/sample correspondence, and annotation requirements.
3. Run validation before resource-intensive alignment and downstream analysis wherever possible.
4. Return actionable errors that identify the relevant sample, field, file, or design issue.
5. Add negative tests for malformed inputs and regression coverage for the existing valid paired-end test dataset.

**Acceptance criteria:** supported valid inputs continue to run unchanged; representative invalid inputs stop early with specific, actionable messages.

### Phase 2 — Record run provenance

1. Identify the effective workflow parameters and available Nextflow, pipeline, and process-tool versions.
2. Record reference and annotation identifiers and checksums where they can be obtained reliably.
3. Record the exact GO annotation and ontology inputs used, including downloaded resource identity rather than only moving `current` URLs.
4. Publish a machine-readable provenance file with the run outputs and document its fields.
5. Test that the file is valid, parseable, and present in the test run.

**Acceptance criteria:** a user can inspect an output artifact to determine the effective analysis configuration and external resources used, with unavailable metadata represented explicitly rather than guessed.

### Phase 3 — Add multi-factor designs and explicit contrasts

1. Specify a samplesheet or parameter interface for covariates and contrast definitions.
2. Validate categorical/numeric covariates, missing values, sample alignment, and model identifiability; reject confounded or unsupported designs with actionable messages.
3. Implement the same documented design intent for DESeq2 and edgeR and define how effect direction and output columns correspond.
4. Add simple and multi-factor examples plus focused tests, including an invalid/confounded design.
5. Keep the current two-condition workflow available as the default-compatible path.

**Acceptance criteria:** documented designs fit successfully on representative fixtures, both methods test the requested comparison, and unsupported designs fail clearly.

### Phase 4 — Support single-end reads

1. Define whether samplesheet rows may mix single-end and paired-end reads; default to rejecting mixed layouts until explicitly designed and tested.
2. Make the second FASTQ input optional only for the supported single-end mode.
3. Route each layout through appropriate trimming, raw/trimmed QC, alignment, and output handling.
4. Add single-end test data or a small fixture and retain paired-end integration regression tests.
5. Update samplesheet examples, parameters, and output documentation.

**Acceptance criteria:** single-end input completes the supported workflow and emits expected outputs; paired-end behavior remains unchanged; mixed-layout behavior is explicit and validated.

## Cross-cutting validation

- Run the Nextflow stub profile after changes to workflow wiring.
- Run the existing integration test profile when available, plus targeted tests for each new validation or analysis behavior.
- Check documentation examples against the actual supported interface.
- Preserve and compare the existing paired-end test outputs and expected output paths.
- Review external-data and software-version assumptions before claiming reproducibility.

## Risks and decisions to resolve during implementation

- Checksums or tool-version metadata may not be available for every remote resource or execution environment; record limitations explicitly.
- Model formulas, factor reference levels, and contrast direction must be specified clearly to avoid misleading differential-expression results.
- Mixed single-end and paired-end runs require additional process and samplesheet semantics; do not infer a layout silently.
- Input validation must happen early without duplicating or contradicting checks inside the R analysis scripts.

## Relationship to the roadmap

See [ROADMAP.md](ROADMAP.md) for the ranked overview, rationale, and completion signals. This plan expands those initiatives into implementation phases and cross-cutting acceptance criteria.
