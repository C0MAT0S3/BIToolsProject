nextflow.enable.dsl=2

process BUILD_KMCP_DB {
    tag "kmcp_db"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "resources/kmcp/${file(fn).name}" } )

    input:
    path ref_fasta

    output:
    path("genomes.kmcp"), emit: db_dir

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p genomes_src
    cp ${ref_fasta} genomes_src/reference.fasta
    kmcp compute -k 21 --split-number 10 --split-overlap 150 --in-dir genomes_src/ --out-dir genomes-k21-n10 --threads ${threads}
    kmcp index --false-positive-rate 0.1 --num-hash 1 --in-dir genomes-k21-n10/ --out-dir genomes.kmcp
    """
}


