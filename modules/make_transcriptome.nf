


process MAKE_TRANSCRIPTOME {
    label 'process_low'

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hd03093a_1' :
        'quay.io/biocontainers/gffread:0.12.7--hd03093a_1' }"

    input:
    path(gtf)
    path(fasta)
    val organism
    val ensembl_version

    output:
    path("*.fa"), emit: transcripts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    gffread \\
        -w ${prefix}.transcripts.fa \\
        -g $fasta \\
        $gtf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    touch ${prefix}.transcripts.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: 0.12.7
    END_VERSIONS
    """
}