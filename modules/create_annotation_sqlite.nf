

process CREATE_ANNOTATION_SQLITE {
    tag "$gtf"
    label 'process_medium'

    publishDir "${params.base}/${organism}/${version}", mode: 'copy'

    conda "bioconda::gffutils=0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffutils:0.12--pyh7cba7a3_0' :
        'quay.io/biocontainers/gffutils:0.12--pyh7cba7a3_0' }"

    input:
    path gtf
    val organism
    val version

    output:
    path "*.sqlite", emit: sqlite_db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = gtf.baseName
    """
    # Create SQLite database from GTF using gffutils
    # TODO: Adjust merge_strategy and other parameters based on your needs
    # Common merge strategies: 'create_unique', 'merge', 'error', 'warning', 'replace'

    python << EOF
import gffutils
import sys

try:
    # Create database from GTF file
    # disable_infer_genes=True if your GTF already has gene features
    # disable_infer_transcripts=True if your GTF already has transcript features
    db = gffutils.create_db(
        data='$gtf',
        dbfn='${prefix}.sqlite',
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True,
        disable_infer_genes=False,
        disable_infer_transcripts=False
    )

    # Print summary statistics
    print(f"Database created successfully: ${prefix}.sqlite")
    print(f"Total features: {db.count_features_of_type()}")

except Exception as e:
    print(f"Error creating database: {e}", file=sys.stderr)
    sys.exit(1)
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffutils: \$(python -c "import gffutils; print(gffutils.__version__)")
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = gtf.baseName
    """
    touch ${prefix}.sqlite

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffutils: 0.12
        python: 3.9.0
    END_VERSIONS
    """
}