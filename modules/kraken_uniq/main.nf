nextflow.enable.dsl=2

process RUN_KRAKENUNIQ_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/krakenuniq/${file(fn).name}" } )

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
    def remain = "krakenuniq/${sample_id}.clean.fastq"
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p krakenuniq
    krakenuniq --db ${krakenuniq_db} --threads ${threads} --report-file ${report} --output ${output} --paired ${reads1} ${reads2} --hll-precision 12
    # Extract reads not classified as host (taxid ${host_taxid})
    krakenuniq-extract-reads -p ${host_taxid} ${output} ${reads1} ${reads2} > ${remain}
    # Split paired reads back and gzip
    awk -v out1="krakenuniq/${sample_id}.clean_R1.fastq" -v out2="krakenuniq/${sample_id}.clean_R2.fastq" '{ print > ((NR%8>0 && NR%8<=4) ? out1 : out2) }' ${remain}
    pigz -p ${threads} krakenuniq/${sample_id}.clean_R1.fastq krakenuniq/${sample_id}.clean_R2.fastq || true
    # Note: host reads extraction as a separate step isn't trivially supported; emit empty gz files for now
    : > krakenuniq/${sample_id}.host_R1.fastq
    : > krakenuniq/${sample_id}.host_R2.fastq
    pigz -p ${threads} krakenuniq/${sample_id}.host_R1.fastq krakenuniq/${sample_id}.host_R2.fastq
    """
}


