

process DOWNLOAD_GENOME_URL {
    tag "${organism}"

    conda "conda-forge::curl=8.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' }"

    publishDir "${params.base}/${organism}/${version}", mode: 'copy'

    input:
    val(genome_fasta_url)
    val(genome_gtf_url)
    val(organism)
    val(version)

    output:
    path "*.fa", emit: genome_fasta
    path "*.gtf", emit: genome_gtf
    path "versions.yml", emit: versions

    script:
    def fasta_filename = genome_fasta_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    def gtf_filename = genome_gtf_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    """
    # Download genome FASTA
    curl --silent -o genome.fa.gz ${genome_fasta_url}

    # Download GTF annotation
    curl --silent -o genome.gtf.gz ${genome_gtf_url}

    # Decompress files
    gzip -d genome.fa.gz
    gzip -d genome.gtf.gz

    # Rename to clean names
    mv genome.fa ${fasta_filename}
    mv genome.gtf ${gtf_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n1 | sed 's/curl //')
    END_VERSIONS
    """

    stub:
    def fasta_filename = genome_fasta_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    def gtf_filename = genome_gtf_url.tokenize('/')[-1].replaceAll(/\.gz$/, '')
    """
    touch ${fasta_filename}
    touch ${gtf_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.21.3
    END_VERSIONS
    """
}
