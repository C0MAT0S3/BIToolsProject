nextflow.enable.dsl=2

process BUILD_BOWTIE2_INDEX {
    tag "bowtie2_index"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "resources/bowtie2/${file(fn).name}" } )

    input:
    path ref_fasta

    output:
    path("bowtie2_index"), emit: index_dir

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p bowtie2_index
    cp ${ref_fasta} bowtie2_index/reference.fasta
    bowtie2-build -f bowtie2_index/reference.fasta bowtie2_index/index -p ${threads}
    """
}


