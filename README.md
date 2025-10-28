# Organism Reference Genome Setup Pipeline

A Nextflow pipeline for downloading, processing, and indexing reference genomes for [LAPTI-UCC/riboseqorg-nf](https://github.com/LAPTI-UCC/riboseqorg-nf).

## Features

## Quick Start

### 1. Basic Usage

Run the pipeline with default settings (Yeast genome):

```bash
nextflow run main.nf -c params.config
```

### 2. Custom Organism

Download and process a different organism:

```bash
nextflow run main.nf -c params.config \
    --organism "Mus musculus" \
    --ensembl_version 112
```

### 3. Stub Testing

Test the workflow structure without running bioinformatics tools:

```bash
nextflow run main.nf -stub -c params.config
```


## Configuration

### params.config

The main configuration file. Key parameters:

```groovy
params {
    // Output directory
    base = "results"

    // Download method: 'gget' or 'url'
    download_method = "gget"

    // Organism (for gget method)
    organism = "Saccharomyces cerevisiae"
    ensembl_version = '112'

    // rRNA source: 'gtf' or 'silva'
    rrna_source = "gtf"
}
```


## Automatic riboseqorg-nf Config Generation

The pipeline automatically generates a `riboseqorg_params.config` file compatible with the [riboseqorg-nf](https://github.com/LAPTI-UCC/riboseqorg-nf) pipeline. This config is created at the end of every run and contains all necessary reference paths.

### Output Location

```
results/<organism_normalized>/<ensembl_version>/riboseqorg_params.config
```

### Using with riboseqorg-nf

```bash
# Clone the riboseqorg-nf pipeline
git clone https://github.com/LAPTI-UCC/riboseqorg-nf.git
cd riboseqorg-nf

# Copy the generated config
cp /path/to/results/saccharomyces_cerevisiae/112/riboseqorg_params.config ./params.config

# Run the pipeline
nextflow run main.nf -c params.config --sample_sheet samples.csv
```

## Advanced Usage

### Custom URLs (T2T-CHM13 Example)

```groovy
// In params.config
params {
    download_method = "url"
    organism = "Homo_sapiens_T2T"
    genome_fasta_url = "https://example.com/chm13v2.0.fa.gz"
    genome_gtf_url = "https://example.com/chm13v2.0.gtf.gz"
}
```

### SILVA rRNA Database

```groovy
// In params.config
params {
    rrna_source = "silva"
    silva_url = "https://www.arb-silva.de/.../SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
    silva_version = "138.1"
}
```

## Module Documentation

### Available Processes

- `DOWNLOAD_GENOME_GGET` - Download from Ensembl via gget
- `DOWNLOAD_GENOME_URL` - Download from custom URLs
- `PREPARE_STAR_GENOME_INDEX` - Build STAR index
- `MAKE_TRANSCRIPTOME` - Extract transcriptome from genome
- `BUILD_BOWTIE_INDEX` - Build Bowtie index
- `BUILD_BOWTIE_INDEX_GENERIC` - Build Bowtie index with custom name
- `GENERATE_CHROM_SIZES` - Generate chromosome sizes file
- `CREATE_ANNOTATION_SQLITE` - Create SQLite database from GTF
- `RIBOMETRIC_PREPARE` - Generate RiboMetric annotations
- `EXTRACT_RRNA_FROM_GTF` - Extract rRNA sequences from GTF
- `DOWNLOAD_RRNA_SILVA` - Download rRNA from SILVA database

### All Modules Have Stub Support

Every module includes stub tests for rapid workflow validation:

```bash
nextflow run main.nf -stub -c params.config
```


## License

MIT License - See LICENSE file for details
