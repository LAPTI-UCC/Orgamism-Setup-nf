process PREPARE_STAR_GENOME_INDEX {
    tag "$genome"
    label 'process_high'

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'


    conda "bioconda::star=2.7.10b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.10b--h9ee0642_0' :
        'biocontainers/star:2.7.10b--h9ee0642_0' }"

    input:
    path genome

    output:
    path 'star_index', emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir -p star_index

    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome} \\
        --runThreadN ${task.cpus} \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}