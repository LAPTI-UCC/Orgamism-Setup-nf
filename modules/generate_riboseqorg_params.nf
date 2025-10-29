

process GENERATE_RIBOSEQORG_PARAMS {
    tag "${organism}_${ensembl_version}"
    label 'process_low'

    publishDir "${params.base}/${organism}/${ensembl_version}", mode: 'copy'

    input:
    path star_index
    path bowtie_index
    path rrna_index
    path gtf
    path genome_fasta
    path chrom_sizes
    path ribometric_tsv
    // path annotation_sqlite
    val organism
    val ensembl_version

    output:
    path "riboseqorg_params.config", emit: config
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_file = "riboseqorg_params.config"
    def base_path = "${params.base}/${organism}/${ensembl_version}"
    // Get absolute path for the base directory
    def absolute_base = new File(base_path).absolutePath
    """
    generate_riboseqorg_config.py \\
        --star_index ${absolute_base}/${star_index} \\
        --bowtie_index ${absolute_base}/bowtie_index \\
        --rrna_index ${absolute_base}/${rrna_index} \\
        --gtf ${absolute_base}/${gtf} \\
        --genome_fasta ${absolute_base}/${genome_fasta} \\
        --chrom_sizes ${absolute_base}/${chrom_sizes} \\
        --ribometric_tsv ${absolute_base}/${ribometric_tsv} \\
        # --annotation_sqlite ${absolute_base}/${annotation_sqlite} \\
        --organism "${organism}" \\
        --ensembl_version "${ensembl_version}" \\
        --output ${output_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def base_path = "${params.base}/${organism}/${ensembl_version}"
    def absolute_base = new File(base_path).absolutePath
    """
    cat > riboseqorg_params.config <<EOF
// Auto-generated params.config for riboseqorg-nf pipeline (STUB)
// Organism: ${organism}
// Ensembl Version: ${ensembl_version}

params {
    // Reference Genome Paths
    star_index = "${absolute_base}/${star_index}"
    bowtie_index = "${absolute_base}/bowtie_index"
    rRNA_index = "${absolute_base}/${rrna_index}"
    gtf = "${absolute_base}/${gtf}"
    genome_fasta = "${absolute_base}/${genome_fasta}"
    chrom_sizes_file = "${absolute_base}/${chrom_sizes}"
    ribometric_annotation = "${absolute_base}/${ribometric_tsv}"
    #annotation_sqlite = "${absolute_base}/${annotation_sqlite}"

    // Alignment Parameters
    mismatches = 3
    alignment_type = "EndToEnd"
    allow_introns = true
    max_multimappers = 10

    // Sample Sheet
    sample_sheet = "\\\${projectDir}/sample_sheet.csv"
}
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9.0
    END_VERSIONS
    """
}
