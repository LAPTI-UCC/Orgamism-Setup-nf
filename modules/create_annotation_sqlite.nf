

process BUILD_BOWTIE_INDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::bowtie=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hed8969a_1' :
        'quay.io/biocontainers/bowtie:1.3.0--py38hed8969a_1' }"

    input:
    path fasta

    output:
    path "bowtie_index"    , emit: index
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie_index
    bowtie-build $args --threads $task.cpus $fasta bowtie_index/${fasta.baseName}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}