#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import processes
include { PREPARE_STAR_GENOME_INDEX } from './modules/build_STAR'
include { BUILD_BOWTIE_INDEX } from './modules/build_bowtie'
// include { DOWNLOAD_GENOME_ENSEMBL } from './modules/ensembl_download'
include { DOWNLOAD_GENOME_GGET } from './modules/gget_download'
include { MAKE_TRANSCRIPTOME } from './modules/make_transcriptome'
include { RIBOMETRIC_PREPARE } from './modules/RiboMetric_prepare'



workflow {
    main:
        if (params.download_method == 'gget') {
            DOWNLOAD_GENOME_GGET(
                params.organism.toLowerCase().replace(' ', '_'),
                params.ensembl_version,
                params.gget_which
            )
            genome_fasta = DOWNLOAD_GENOME_GGET.out.genome_fasta
            genome_gtf = DOWNLOAD_GENOME_GGET.out.genome_gtf
        } else {
            error "Invalid download method. Choose 'ensembl' or 'gget'."
        }

        PREPARE_STAR_GENOME_INDEX(genome_fasta)
        MAKE_TRANSCRIPTOME(genome_gtf, genome_fasta)
        BUILD_BOWTIE_INDEX(MAKE_TRANSCRIPTOME.out.transcripts)
}

// Function to print help message
def helpMessage() {
    log.info"""
    Usage:ll
    nextflow run main.nf --organism "Homo sapiens" --ensembl_version 108 --download_method ensembl

    Required arguments:
      --organism          Organism name (e.g., "Homo sapiens")
      --ensembl_version   Ensembl version number
      --download_method   Method to download genome: 'ensembl' or 'gget'

    Optional arguments:
      --outdir            Output directory (default: 'results')
      --gget_which        List of file types to download when using gget (default: ['dna', 'gtf'])
    """.stripIndent()
}

// Show help message if --help flag is used
if (params.help) {
    helpMessage()
    exit 0
}