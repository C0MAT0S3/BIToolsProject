nextflow.enable.dsl=2

process RUN_KMCP_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> file(fn).name.endsWith('.fastq.gz') ? "${sample_id}/reads/kmcp/${file(fn).name}" : null } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path kmcp_db_dir

    output:
    tuple val(sample_id), path("kmcp/${sample_id}.clean_R1.fastq.gz"), path("kmcp/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("kmcp/${sample_id}.host_R1.fastq.gz"),  path("kmcp/${sample_id}.host_R2.fastq.gz"),  emit: host
    tuple val(sample_id), path("kmcp/${sample_id}.host_ids.txt"), emit: host_ids

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p kmcp
    kmcp search --db-dir ${kmcp_db_dir} -1 ${reads1} -2 ${reads2} --out-file kmcp/${sample_id}.tsv.gz --load-whole-db --threads ${threads}
    gunzip -c kmcp/${sample_id}.tsv.gz | grep -v '^#' | cut -f1 | sort -u > kmcp/${sample_id}.host_ids.txt
    seqkit grep -v -f kmcp/${sample_id}.host_ids.txt ${reads1} | pigz -p ${threads} > kmcp/${sample_id}.clean_R1.fastq.gz
    seqkit grep -v -f kmcp/${sample_id}.host_ids.txt ${reads2} | pigz -p ${threads} > kmcp/${sample_id}.clean_R2.fastq.gz
    seqkit grep    -f kmcp/${sample_id}.host_ids.txt ${reads1} | pigz -p ${threads} > kmcp/${sample_id}.host_R1.fastq.gz
    seqkit grep    -f kmcp/${sample_id}.host_ids.txt ${reads2} | pigz -p ${threads} > kmcp/${sample_id}.host_R2.fastq.gz
    """
}


