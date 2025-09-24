nextflow.enable.dsl=2

process BUILD_KRAKEN2_DB {
    tag "kraken2_db"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "resources/kraken2/${file(fn).name}" } )

    input:
    path ref_fasta
    path taxonomy_dir
    val taxid

    output:
    path("kraken2_db"), emit: db_dir

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p kraken2_db/library kraken2_db/taxonomy
    cp -r ${taxonomy_dir}/* kraken2_db/taxonomy/
    awk -v tx=${taxid} '/^>/{sub(">","",\$1); print ">" \$1 "|kraken:taxid|" tx; next} {print}' ${ref_fasta} > kraken2_db/reference.tax.fasta
    kraken2-build --db kraken2_db --add-to-library kraken2_db/reference.tax.fasta --threads ${threads}
    kraken2-build --db kraken2_db --build --threads ${threads}
    """
}


