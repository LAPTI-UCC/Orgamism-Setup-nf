#!/usr/bin/env python3

"""
Generate params.config for riboseqorg-nf pipeline
"""

import os
import sys
import argparse
from pathlib import Path
from datetime import datetime


def get_absolute_path(base_path, filename):
    """Construct absolute path from base path and filename."""
    return str(Path(base_path) / filename)


def generate_riboseqorg_config(
    star_index,
    bowtie_index,
    rrna_index,
    gtf,
    genome_fasta,
    chrom_sizes,
    ribometric_tsv,
    annotation_sqlite,
    organism,
    ensembl_version,
    output_file
):
    """
    Generate a params.config file for the riboseqorg-nf pipeline.

    Args:
        star_index: Path to STAR index directory
        bowtie_index: Path to Bowtie transcriptome index directory
        rrna_index: Path to rRNA Bowtie index directory
        gtf: Path to GTF annotation file
        genome_fasta: Path to genome FASTA file
        chrom_sizes: Path to chromosome sizes file
        ribometric_tsv: Path to RiboMetric annotation TSV
        annotation_sqlite: Path to annotation SQLite database
        organism: Organism name
        ensembl_version: Ensembl version
        output_file: Output path for generated config
    """

    # Get absolute paths
    star_index_abs = str(Path(star_index).resolve())
    bowtie_index_abs = str(Path(bowtie_index).resolve())
    rrna_index_abs = str(Path(rrna_index).resolve())
    gtf_abs = str(Path(gtf).resolve())
    genome_fasta_abs = str(Path(genome_fasta).resolve())
    chrom_sizes_abs = str(Path(chrom_sizes).resolve())
    ribometric_tsv_abs = str(Path(ribometric_tsv).resolve())
    annotation_sqlite_abs = str(Path(annotation_sqlite).resolve())

    # Generate config content
    config_content = f"""// Auto-generated params.config for riboseqorg-nf pipeline
// Organism: {organism}
// Ensembl Version: {ensembl_version}
// Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

params {{
    // =====================================
    // Output Directory
    // =====================================
    outdir = "./results"

    // =====================================
    // Reference Genome Paths
    // =====================================
    // STAR index for genome alignment
    star_index = "{star_index_abs}"

    // Bowtie index for transcriptome alignment
    bowtie_index = "{bowtie_index_abs}"

    // Bowtie index for rRNA filtering
    rRNA_index = "{rrna_index_abs}"

    // GTF annotation file
    gtf = "{gtf_abs}"

    // Genome FASTA file
    genome_fasta = "{genome_fasta_abs}"

    // Chromosome sizes file
    chrom_sizes_file = "{chrom_sizes_abs}"

    // RiboMetric annotation file
    ribometric_annotation = "{ribometric_tsv_abs}"

    // Annotation SQLite database
    annotation_sqlite = "{annotation_sqlite_abs}"

    // =====================================
    // Alignment Parameters
    // =====================================
    // Maximum number of mismatches allowed in alignment
    mismatches = 3

    // Alignment type: "EndToEnd" or "Local"
    alignment_type = "EndToEnd"

    // Allow alignments spanning introns
    allow_introns = true

    // Maximum number of multimapper alignments to keep
    max_multimappers = 10

    // =====================================
    // Sample Sheet
    // =====================================
    // Path to sample sheet CSV file
    // Required columns: sample_name, fastq (or fastq_1, fastq_2 for paired-end)
    sample_sheet = "${{projectDir}}/sample_sheet.csv"

    // =====================================
    // Adapter Finding
    // =====================================
    // Automatically detect adapter sequences using architecture detection
    use_architecture_detection = false

    // =====================================
    // Data Fetching
    // =====================================
    // Fetch data from external sources (e.g., SRA)
    fetch = true

    // Force re-fetching of data even if it already exists
    force_fetch = false

    // Path to collapsed read files (optional)
    // collapsed_read_path = "/path/to/collapsed/reads"

    // Header pattern for collapsed reads (optional)
    // collapsed_read_header_pattern = ".*"

    // =====================================
    // Output Options
    // =====================================
    // Save BAM files from STAR transcriptome alignment
    save_star_transcriptome_bam = true

    // =====================================
    // Trimming
    // =====================================
    // Number of bases to trim from the front of reads
    trim_front = 0
}}
"""

    # Write to output file
    with open(output_file, 'w') as f:
        f.write(config_content)

    print(f"✓ Generated riboseqorg-nf params.config: {output_file}")
    print(f"  Organism: {organism}")
    print(f"  Version: {ensembl_version}")
    print(f"  STAR index: {star_index_abs}")
    print(f"  Bowtie index: {bowtie_index_abs}")
    print(f"  rRNA index: {rrna_index_abs}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate params.config for riboseqorg-nf pipeline'
    )

    parser.add_argument('--star_index', required=True,
                        help='Path to STAR index directory')
    parser.add_argument('--bowtie_index', required=True,
                        help='Path to Bowtie transcriptome index directory')
    parser.add_argument('--rrna_index', required=True,
                        help='Path to rRNA Bowtie index directory')
    parser.add_argument('--gtf', required=True,
                        help='Path to GTF annotation file')
    parser.add_argument('--genome_fasta', required=True,
                        help='Path to genome FASTA file')
    parser.add_argument('--chrom_sizes', required=True,
                        help='Path to chromosome sizes file')
    parser.add_argument('--ribometric_tsv', required=True,
                        help='Path to RiboMetric annotation TSV')
    parser.add_argument('--annotation_sqlite', required=True,
                        help='Path to annotation SQLite database')
    parser.add_argument('--organism', required=True,
                        help='Organism name')
    parser.add_argument('--ensembl_version', required=True,
                        help='Ensembl version')
    parser.add_argument('--output', required=True,
                        help='Output path for generated config')

    args = parser.parse_args()

    generate_riboseqorg_config(
        star_index=args.star_index,
        bowtie_index=args.bowtie_index,
        rrna_index=args.rrna_index,
        gtf=args.gtf,
        genome_fasta=args.genome_fasta,
        chrom_sizes=args.chrom_sizes,
        ribometric_tsv=args.ribometric_tsv,
        annotation_sqlite=args.annotation_sqlite,
        organism=args.organism,
        ensembl_version=args.ensembl_version,
        output_file=args.output
    )


if __name__ == '__main__':
    main()
