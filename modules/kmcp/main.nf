nextflow.enable.dsl=2

process RUN_KMCP_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/kmcp/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path kmcp_db_dir

    output:
    tuple val(sample_id), path("kmcp/${sample_id}.clean_R1.fastq.gz"), path("kmcp/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("kmcp/${sample_id}.host_ids.txt"), emit: host_ids

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p kmcp
    kmcp search --db-dir ${kmcp_db_dir} -1 ${reads1} -2 ${reads2} --out-file kmcp/${sample_id}.tsv.gz --load-whole-db --threads ${threads}
    gunzip -c kmcp/${sample_id}.tsv.gz | grep -v '^#' | cut -f1 | sort -u > kmcp/${sample_id}.host_ids.txt
    # Filter reads not in host_ids to produce cleaned fastqs (IDs assumed to match first column without leading @)
    seqkit grep -v -f kmcp/${sample_id}.host_ids.txt ${reads1} | pigz -p ${threads} > kmcp/${sample_id}.clean_R1.fastq.gz
    seqkit grep -v -f kmcp/${sample_id}.host_ids.txt ${reads2} | pigz -p ${threads} > kmcp/${sample_id}.clean_R2.fastq.gz
    """
}


