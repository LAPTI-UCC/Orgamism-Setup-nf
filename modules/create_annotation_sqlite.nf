

process CREATE_ANNOTATION_SQLITE {
    tag "$organism"
    label 'process_medium'

    publishDir "${params.base}/${organism}/${version}", mode: 'copy'

    conda "conda-forge::python=3.9 conda-forge::intervaltree=3.1.0 conda-forge::sqlitedict=2.1.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
    //     'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path gtf
    path transcriptome_fasta
    val organism
    val version

    output:
    path "*.sqlite", emit: sqlite_db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def pseudo_utr_len = task.ext.pseudo_utr_len ?: 0
    """
    create_annotation_sqlite.py \\
        ${organism} \\
        ${gtf} \\
        ${transcriptome_fasta} \\
        ${pseudo_utr_len}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        intervaltree: \$(python -c "import intervaltree; print(intervaltree.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = gtf.baseName
    """
    touch ${prefix}.sqlite

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffutils: 0.12
        python: 3.9.0
    END_VERSIONS
    """
}