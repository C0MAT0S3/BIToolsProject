nextflow.enable.dsl=2

process BUILD_BWA_INDEX {
    tag "bwa_index"
    cpus 4
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "resources/bwa/${file(fn).name}" } )

    input:
    path ref_fasta

    output:
    path("bwa_index"), emit: index_dir

    script:
    """
    set -euo pipefail
    mkdir -p bwa_index
    cp ${ref_fasta} bwa_index/reference.fasta
    bwa index bwa_index/reference.fasta
    """
}


