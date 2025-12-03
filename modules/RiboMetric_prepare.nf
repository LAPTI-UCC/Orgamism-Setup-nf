

process RIBOMETRIC_PREPARE {
    label 'process_high'

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'

    conda "${projectDir}/conda/ribometric.yml"
    container "ghcr.io/lapti-ucc/riboseqorg-nf-ribometric:latest"

    input:
    path(gtf)
    path(fasta)
    val organism
    val ensembl_version

    output:
    path("*.tsv"), emit: ribometric_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    RiboMetric prepare -p ${task.cpus} -g $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    touch ${prefix}_ribometric.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: 0.12.7
    END_VERSIONS
    """
}