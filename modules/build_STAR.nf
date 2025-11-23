process PREPARE_STAR_GENOME_INDEX {
    tag "$genome"
    label 'process_high'

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'


    conda "bioconda::star=2.7.11b"

    input:
    path genome
    val organism
    val ensembl_version

    output:
    path 'star_index', emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Only set memory limit if explicitly provided and reasonable (>= 8GB)
    // Otherwise let STAR use what it needs
    def memory = (task.memory && task.memory.toGiga() >= 8) ?
        "--limitGenomeGenerateRAM ${task.memory.toBytes()}" : ''
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

    stub:
    """
    mkdir -p star_index
    touch star_index/SA
    touch star_index/SAindex
    touch star_index/Genome
    touch star_index/chrName.txt
    touch star_index/chrLength.txt
    touch star_index/chrStart.txt
    touch star_index/chrNameLength.txt
    touch star_index/genomeParameters.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: 2.7.10b
    END_VERSIONS
    """
}