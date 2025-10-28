# Examples

This directory contains example configurations and outputs from the pipeline.

## Example Input Configuration Files

These demonstrate different ways to configure the genome setup pipeline for various organisms and data sources.

### Default: params.config (in project root)
- **Organism**: Saccharomyces cerevisiae (yeast)
- **Download method**: gget (Ensembl)
- **rRNA source**: GTF extraction
- **Use case**: Quick testing with small genome (~12 Mb)

### params_mouse_ensembl.config
- **Organism**: Mus musculus (mouse)
- **Download method**: gget (Ensembl release 112)
- **rRNA source**: GTF extraction
- **Use case**: Standard mammalian genome from Ensembl

Usage:
```bash
nextflow run main.nf -c examples/params_mouse_ensembl.config
```

### params_t2t_silva.config
- **Organism**: Homo sapiens T2T-CHM13 v2.0
- **Download method**: URL (custom download from S3/NCBI)
- **rRNA source**: SILVA SSU database (16S/18S)
- **Use case**: Complete telomere-to-telomere human genome with external rRNA database

Usage:
```bash
nextflow run main.nf -c examples/params_t2t_silva.config
```

### params_zebrafish_url_silva.config
- **Organism**: Danio rerio (zebrafish)
- **Download method**: URL (NCBI GRCz11)
- **rRNA source**: SILVA LSU database (23S/28S)
- **Use case**: NCBI genome with SILVA LSU rRNA for ribosome profiling studies

Usage:
```bash
nextflow run main.nf -c examples/params_zebrafish_url_silva.config
```

## Generated Output Configuration

### riboseqorg_params.config.example

An example of the automatically generated params.config file for the [riboseqorg-nf](https://github.com/LAPTI-UCC/riboseqorg-nf) pipeline.

This file is generated automatically at the end of each pipeline run and contains:
- Absolute paths to all reference files (STAR index, Bowtie indices, GTF, FASTA, etc.)
- Default alignment parameters
- Template configuration for sample sheets and other pipeline options

### How It Works

1. **Run the genome setup pipeline** with your desired organism:
   ```bash
   nextflow run main.nf -c params.config \
       --organism "Saccharomyces cerevisiae" \
       --ensembl_version 112
   ```

2. **Pipeline generates** `riboseqorg_params.config` automatically:
   ```
   results/saccharomyces_cerevisiae/112/riboseqorg_params.config
   ```

3. **Use with riboseqorg-nf**:
   ```bash
   # Clone riboseqorg-nf
   git clone https://github.com/LAPTI-UCC/riboseqorg-nf.git
   cd riboseqorg-nf

   # Copy generated config
   cp /path/to/results/saccharomyces_cerevisiae/112/riboseqorg_params.config ./params.config

   # Run analysis
   nextflow run main.nf -c params.config --sample_sheet samples.csv
   ```

## Generated Config Structure

The generated config includes all parameters needed for riboseqorg-nf:

```groovy
params {
    // Reference paths (automatically populated)
    star_index = "/absolute/path/to/star_index"
    bowtie_index = "/absolute/path/to/bowtie_index"
    rRNA_index = "/absolute/path/to/rRNA_index"
    gtf = "/absolute/path/to/annotation.gtf"
    genome_fasta = "/absolute/path/to/genome.fa"
    chrom_sizes_file = "/absolute/path/to/chrom.sizes"
    ribometric_annotation = "/absolute/path/to/ribometric.tsv"
    annotation_sqlite = "/absolute/path/to/annotation.sqlite"

    // Default parameters (customize as needed)
    mismatches = 3
    alignment_type = "EndToEnd"
    allow_introns = true
    max_multimappers = 10

    // User-provided (update these)
    sample_sheet = "${projectDir}/sample_sheet.csv"
}
```

## Python Script

The config generation is handled by `bin/generate_riboseqorg_config.py`, which:
- Takes paths to all reference files as input
- Generates a properly formatted params.config
- Uses absolute paths for portability
- Includes all required parameters for riboseqorg-nf

## Module

The Nextflow module `modules/generate_riboseqorg_params.nf`:
- Runs at the end of the main workflow
- Collects outputs from all reference generation processes
- Calls the Python script to generate the config
- Publishes the config alongside other reference files

## Customization

You can modify the generated config to:
- Change alignment parameters (mismatches, alignment_type, etc.)
- Enable/disable features (use_architecture_detection, fetch, etc.)
- Set trimming options (trim_front)
- Configure output options (save_star_transcriptome_bam)

All reference file paths will already be set correctly and don't need to be changed.
