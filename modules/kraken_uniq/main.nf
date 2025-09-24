nextflow.enable.dsl=2

process RUN_KRAKENUNIQ_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/reads/krakenuniq/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path krakenuniq_db
    val host_taxid

    output:
    tuple val(sample_id), path("krakenuniq/${sample_id}.clean_R1.fastq.gz"), path("krakenuniq/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("krakenuniq/${sample_id}.host_R1.fastq.gz"),  path("krakenuniq/${sample_id}.host_R2.fastq.gz"),  emit: host
    
    script:
    def report = "krakenuniq/${sample_id}.report"
    def output = "krakenuniq/${sample_id}.output"
    def remain1 = "krakenuniq/${sample_id}.clean_R1.fastq"
    def remain2 = "krakenuniq/${sample_id}.clean_R2.fastq"
    def contam1 = "krakenuniq/${sample_id}.host_R1.fastq"
    def contam2 = "krakenuniq/${sample_id}.host_R2.fastq"
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p krakenuniq
    krakenuniq --db ${krakenuniq_db} --threads ${threads} --report-file ${report} --output ${output} --paired ${reads1} ${reads2}
    # Use krakenuniq-extract-reads (per-mate) with taxonomy children
    krakenuniq-extract-reads -t ${krakenuniq_db}/taxDB ${host_taxid} ${output} ${reads1} > ${contam1}
    krakenuniq-extract-reads -t ${krakenuniq_db}/taxDB ${host_taxid} ${output} ${reads2} > ${contam2}
    krakenuniq-extract-reads -i -t ${krakenuniq_db}/taxDB ${host_taxid} ${output} ${reads1} > ${remain1}
    krakenuniq-extract-reads -i -t ${krakenuniq_db}/taxDB ${host_taxid} ${output} ${reads2} > ${remain2}
    pigz -p ${threads} ${remain1} ${remain2} ${contam1} ${contam2}

    
    """
}
