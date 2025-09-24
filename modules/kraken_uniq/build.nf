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
    cp ${ref_fasta} krakenuniq_db/library/reference.fasta
    awk -v tx=${taxid} '/^>/{hdr=\$1; sub(/^>/,"",hdr); print hdr "\t" tx}' krakenuniq_db/library/reference.fasta > krakenuniq_db/library/reference.map
    krakenuniq-build --db krakenuniq_db --kmer-len 31 --threads ${threads} --taxids-for-genomes --taxids-for-sequences --jellyfish-bin jellyfish
    """
}


