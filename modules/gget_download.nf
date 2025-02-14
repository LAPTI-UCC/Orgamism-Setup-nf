

process DOWNLOAD_GENOME_GGET {
    tag "${params.organism}_${params.ensembl_version}"

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'
    
    conda "bioconda::gget=0.29.0 conda-forge::openmm=7.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gget:0.27.0--pyhdfd78af_0' :
        'quay.io/biocontainers/gget:0.27.0--pyhdfd78af_0' }"

    input:
    val(species)
    val(ensembl_release)
    val(which)

    output:
    path "*.fa", emit: genome_fasta, optional: true
    path "*.gtf", emit: genome_gtf, optional: true
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def which_arg = which ? "--which ${which.join(',')}" : ''
    def release_arg = ensembl_release ? "--release ${ensembl_release}" : ''
    """
    gget ref ${species} \\
        ${which_arg} \\
        ${release_arg} \\
        ${args} \\
        --download 

    gzip -d *

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gget: \$(gget --version | sed 's/gget //')
    END_VERSIONS
    """
}