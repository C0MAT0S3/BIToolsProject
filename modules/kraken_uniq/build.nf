nextflow.enable.dsl=2

process BUILD_KRAKENUNIQ_DB {
    tag "krakenuniq_db"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "resources/krakenuniq/${file(fn).name}" } )

    input:
    path ref_fasta
    path taxonomy_dir
    val taxid

    output:
    path("krakenuniq_db"), emit: db_dir

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p krakenuniq_db/library krakenuniq_db/taxonomy
    cp -r ${taxonomy_dir}/* krakenuniq_db/taxonomy/
    awk -v tx=${taxid} '/^>/{sub(">","",\$1); print ">" \$1 "|kraken:taxid|" tx; next} {print}' ${ref_fasta} > krakenuniq_db/reference.tax.fasta
    cp krakenuniq_db/reference.tax.fasta krakenuniq_db/library/
    # Locate Jellyfish binary from the active environment
    JELLY="\$(command -v jellyfish || true)"
    if [ -z "\${JELLY}" ]; then
        echo "ERROR: jellyfish not found in PATH. Ensure it is installed in the bitools env." >&2
        exit 1
    fi
    krakenuniq-build --db krakenuniq_db --kmer-len 31 \
      --taxids-for-genomes --taxids-for-sequences \
      --build --threads ${threads} --jellyfish-bin "\${JELLY}"
    """
}


