process GO_ANNOTATIONS {
    label 'process_low'
    conda 'conda-forge::curl'
    output:
    path 'go-annotations.gaf.gz', emit: gaf
    path 'go-basic.obo', emit: obo
    script:
    """
    curl --fail --location --retry 3 "$params.go_annotation_url" -o go-annotations.gaf.gz
    curl --fail --location --retry 3 "$params.go_obo_url" -o go-basic.obo
    """
    stub:
    """
    touch go-annotations.gaf.gz go-basic.obo
    """
}
