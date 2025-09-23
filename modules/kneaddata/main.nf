nextflow.enable.dsl=2

process RUN_KNEADDATA_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/kneaddata/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path kneaddata_db

    output:
    tuple val(sample_id), path("kneaddata/${sample_id}.clean_R1.fastq.gz"), path("kneaddata/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("kneaddata/${sample_id}.host_R1.fastq.gz"),  path("kneaddata/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def threads = task.cpus
    def paired1 = "kneaddata/${sample_id}_kneaddata_paired_1.fastq"
    def paired2 = "kneaddata/${sample_id}_kneaddata_paired_2.fastq"
    def contam1 = "kneaddata/${sample_id}_kneaddata_${kneaddata_db.baseName}_bowtie2_paired_contam_1.fastq"
    def contam2 = "kneaddata/${sample_id}_kneaddata_${kneaddata_db.baseName}_bowtie2_paired_contam_2.fastq"
    """
    set -euo pipefail
    mkdir -p kneaddata
    kneaddata -i1 ${reads1} -i2 ${reads2} -o kneaddata -db ${kneaddata_db} --bypass-trim --bypass-trf -t ${threads}
    # Standardize output names
    pigz -p ${threads} -c ${paired1} > kneaddata/${sample_id}.clean_R1.fastq.gz
    pigz -p ${threads} -c ${paired2} > kneaddata/${sample_id}.clean_R2.fastq.gz
    if [ -f "${contam1}" ] && [ -f "${contam2}" ]; then
        pigz -p ${threads} -c ${contam1} > kneaddata/${sample_id}.host_R1.fastq.gz
        pigz -p ${threads} -c ${contam2} > kneaddata/${sample_id}.host_R2.fastq.gz
    else
        : > kneaddata/${sample_id}.host_R1.fastq && : > kneaddata/${sample_id}.host_R2.fastq
        pigz -p ${threads} kneaddata/${sample_id}.host_R1.fastq kneaddata/${sample_id}.host_R2.fastq
    fi
    """
}


