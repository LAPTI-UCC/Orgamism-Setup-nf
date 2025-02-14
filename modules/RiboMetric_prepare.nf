

process RIBOMETRIC_PREPARE {
    label 'process_high'

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'

    conda conda/ribometric.yml
    container "community.wave.seqera.io/library/pip_ribometric:6264b49edef91023"

    input:
    path(gtf)
    path(fasta)

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
}