nextflow.enable.dsl=2

process RUN_KRAKEN2_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/kraken2/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path kraken2_db
    val host_taxid

    output:
    tuple val(sample_id), path("kraken2/${sample_id}.clean_R1.fastq.gz"), path("kraken2/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("kraken2/${sample_id}.host_R1.fastq.gz"),  path("kraken2/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def report = "kraken2/${sample_id}.report"
    def output = "kraken2/${sample_id}.output"
    def remain1 = "kraken2/${sample_id}.clean_R1.fastq"
    def remain2 = "kraken2/${sample_id}.clean_R2.fastq"
    def contam1 = "kraken2/${sample_id}.host_R1.fastq"
    def contam2 = "kraken2/${sample_id}.host_R2.fastq"
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p kraken2
    kraken2 --db ${kraken2_db} -1 ${reads1} -2 ${reads2} --threads ${threads} --use-names --report-zero-counts --report ${report} --output ${output}
    extract_kraken_reads.py -k ${output} -r ${report} -1 ${reads1} -2 ${reads2} -t ${host_taxid} --include-children --exclude --fastq-output -o ${remain1} -o2 ${remain2}
    extract_kraken_reads.py -k ${output} -r ${report} -1 ${reads1} -2 ${reads2} -t ${host_taxid} --include-children --fastq-output -o ${contam1} -o2 ${contam2}
    pigz -p ${threads} ${remain1} ${remain2} ${contam1} ${contam2}
    """
}


