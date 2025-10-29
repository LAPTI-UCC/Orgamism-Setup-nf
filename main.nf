#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { PREPARE_STAR_GENOME_INDEX } from './modules/build_STAR'
include { BUILD_BOWTIE_INDEX } from './modules/build_bowtie'
include { BUILD_BOWTIE_INDEX_GENERIC } from './modules/build_bowtie_generic'
include { EXTRACT_RRNA_FROM_GTF } from './modules/extract_rRNA_from_gtf'
include { DOWNLOAD_RRNA_SILVA } from './modules/download_rRNA_silva'
include { GENERATE_CHROM_SIZES } from './modules/generate_chrom_sizes'
include { CREATE_ANNOTATION_SQLITE } from './modules/create_annotation_sqlite'
// include { DOWNLOAD_GENOME_ENSEMBL } from './modules/ensembl_download'
include { DOWNLOAD_GENOME_GGET } from './modules/gget_download'
include { DOWNLOAD_GENOME_URL } from './modules/url_download'
include { MAKE_TRANSCRIPTOME } from './modules/make_transcriptome'
include { RIBOMETRIC_PREPARE } from './modules/RiboMetric_prepare'
include { GENERATE_RIBOSEQORG_PARAMS } from './modules/generate_riboseqorg_params'



workflow {
    main:
        // Normalize organism name: lowercase and replace spaces with underscores
        def organism_normalized = params.organism.toLowerCase().replace(' ', '_')

        // Set version based on download method
        def version = params.download_method == 'gget' ? params.ensembl_version : (params.genome_version ?: 'custom')

        if (params.download_method == 'gget') {
            DOWNLOAD_GENOME_GGET(
                organism_normalized,
                params.ensembl_version,
                params.gget_which
            )
            genome_fasta = DOWNLOAD_GENOME_GGET.out.genome_fasta
            genome_gtf = DOWNLOAD_GENOME_GGET.out.genome_gtf
        } else if (params.download_method == 'url') {
            DOWNLOAD_GENOME_URL(
                params.genome_fasta_url,
                params.genome_gtf_url,
                organism_normalized,
                version
            )
            genome_fasta = DOWNLOAD_GENOME_URL.out.genome_fasta
            genome_gtf = DOWNLOAD_GENOME_URL.out.genome_gtf
        } else {
            error "Invalid download method. Choose 'gget' or 'url'."
        }

        // Build indices and generate reference files
        PREPARE_STAR_GENOME_INDEX(genome_fasta, organism_normalized, version)
        MAKE_TRANSCRIPTOME(genome_gtf, genome_fasta, organism_normalized, version)
        BUILD_BOWTIE_INDEX(MAKE_TRANSCRIPTOME.out.transcripts)
        GENERATE_CHROM_SIZES(genome_fasta, organism_normalized, version)
        CREATE_ANNOTATION_SQLITE(genome_gtf, MAKE_TRANSCRIPTOME.out.transcripts, organism_normalized, version)
        RIBOMETRIC_PREPARE(genome_gtf, genome_fasta, organism_normalized, version)

        // rRNA handling - choose method based on params.rrna_source
        if (params.rrna_source == 'gtf') {
            EXTRACT_RRNA_FROM_GTF(genome_gtf, genome_fasta, organism_normalized, version)
            rrna_fasta = EXTRACT_RRNA_FROM_GTF.out.rrna_fasta
        } else if (params.rrna_source == 'silva') {
            DOWNLOAD_RRNA_SILVA(params.silva_url, organism_normalized, version)
            rrna_fasta = DOWNLOAD_RRNA_SILVA.out.rrna_fasta
        } else {
            error "Invalid rRNA source. Choose 'gtf' or 'silva'."
        }

        // Build rRNA bowtie index
        BUILD_BOWTIE_INDEX_GENERIC(rrna_fasta, 'rRNA', organism_normalized, version)

        // Generate riboseqorg-nf params.config
        GENERATE_RIBOSEQORG_PARAMS(
            PREPARE_STAR_GENOME_INDEX.out.index,
            BUILD_BOWTIE_INDEX.out.index,
            BUILD_BOWTIE_INDEX_GENERIC.out.index,
            genome_gtf,
            genome_fasta,
            GENERATE_CHROM_SIZES.out.chrom_sizes,
            RIBOMETRIC_PREPARE.out.ribometric_tsv,
            // CREATE_ANNOTATION_SQLITE.out.sqlite_db,
            organism_normalized,
            version
        )
}

// Function to print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --organism "Homo sapiens" --ensembl_version 108 --download_method gget --rrna_source gtf

    Required arguments:
      --download_method   Method to download genome: 'gget' or 'url'
      --rrna_source       Method to obtain rRNA sequences: 'gtf' or 'silva'

    For 'gget' download method:
      --organism          Organism name (e.g., "Homo sapiens")
      --ensembl_version   Ensembl version number
      --gget_which        List of file types to download (default: ['dna', 'gtf'])

    For 'url' download method:
      --organism          Organism name for output directory (e.g., "Homo_sapiens_T2T")
      --genome_fasta_url  URL to genome FASTA file (can be .gz)
      --genome_gtf_url    URL to GTF annotation file (can be .gz)

    For 'gtf' rRNA source:
      Uses GTF annotations to extract rRNA sequences from genome

    For 'silva' rRNA source:
      --silva_url         URL to SILVA rRNA database FASTA (e.g., SSU or LSU)

    Optional arguments:
      --base              Output base directory (default: 'results')
    """.stripIndent()
}

// Show help message if --help flag is used
if (params.help) {
    helpMessage()
    exit 0
}