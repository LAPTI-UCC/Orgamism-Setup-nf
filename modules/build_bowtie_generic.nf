


process BUILD_BOWTIE_INDEX_GENERIC {
    tag "$fasta"
    label 'process_high'

    publishDir "${params.base}/${organism}/${version}", mode: 'copy'

    conda "bioconda::bowtie=1.3.1 conda-forge::python_abi=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie:1.3.0--py38hed8969a_1' :
        'quay.io/biocontainers/bowtie:1.3.0--py38hed8969a_1' }"

    input:
    path fasta
    val index_name
    val organism
    val version

    output:
    path "${index_name}_index", emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir ${index_name}_index
    bowtie-build $args --threads $task.cpus $fasta ${index_name}_index/${index_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${index_name}_index
    touch ${index_name}_index/${index_name}.1.ebwt
    touch ${index_name}_index/${index_name}.2.ebwt
    touch ${index_name}_index/${index_name}.3.ebwt
    touch ${index_name}_index/${index_name}.4.ebwt
    touch ${index_name}_index/${index_name}.rev.1.ebwt
    touch ${index_name}_index/${index_name}.rev.2.ebwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: 1.3.1
    END_VERSIONS
    """
}
