

process DOWNLOAD_RRNA_SILVA {
    tag "SILVA_${params.silva_version}"
    label 'process_low'

    conda "conda-forge::curl=8.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' :
        'biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:1196f19ebf5dec139b02f2b5077300fdd8d0aa94-0' }"

    publishDir "${params.base}/${organism}/${version}", mode: 'copy'

    input:
    val(silva_url)
    val(organism)
    val(version)

    output:
    path "rRNA.fa", emit: rrna_fasta
    path "versions.yml", emit: versions

    script:
    """
    # Download SILVA rRNA database
    # TODO: Update silva_url in params to match your desired SILVA version
    # Common URLs:
    # - SSU (16S/18S): https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
    # - LSU (23S/28S): https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz

    curl --silent -o silva_rrna.fasta.gz ${silva_url}

    # Decompress
    gzip -d silva_rrna.fasta.gz

    # Rename to standard output name
    mv silva_rrna.fasta rRNA.fa

    echo "Downloaded SILVA rRNA sequences: \$(grep -c '^>' rRNA.fa) sequences"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | head -n1 | sed 's/curl //')
    END_VERSIONS
    """

    stub:
    """
    touch rRNA.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.21.3
    END_VERSIONS
    """
}
