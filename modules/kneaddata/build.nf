nextflow.enable.dsl=2

process BUILD_KNEADDATA_DB {
    tag "kneaddata_db"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "resources/kneaddata/${file(fn).name}" } )

    input:
    path ref_fasta

    output:
    path("kneaddata_db"), emit: db_dir

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p kneaddata_db
    cp ${ref_fasta} kneaddata_db/reference.fasta
    bowtie2-build -f kneaddata_db/reference.fasta kneaddata_db/index -p ${threads}
    """
}


