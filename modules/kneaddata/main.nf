nextflow.enable.dsl=2

process RUN_KNEADDATA_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/reads/kneaddata/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path kneaddata_db

    output:
    tuple val(sample_id), path("kneaddata/${sample_id}.clean_R1.fastq.gz"), path("kneaddata/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("kneaddata/${sample_id}.host_R1.fastq.gz"),  path("kneaddata/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p kneaddata
    kneaddata -i1 ${reads1} -i2 ${reads2} -o kneaddata -db ${kneaddata_db} --bypass-trim --bypass-trf -t ${threads} || true
    p1=\$(ls -1 kneaddata/*_kneaddata_paired_1.fastq 2>/dev/null | head -n1 || true)
    p2=\$(ls -1 kneaddata/*_kneaddata_paired_2.fastq 2>/dev/null | head -n1 || true)
    if [ -n "\${p1}" ]; then pigz -p ${threads} -c "\${p1}" > kneaddata/${sample_id}.clean_R1.fastq.gz; else : > kneaddata/${sample_id}.clean_R1.fastq && pigz -p ${threads} kneaddata/${sample_id}.clean_R1.fastq; fi
    if [ -n "\${p2}" ]; then pigz -p ${threads} -c "\${p2}" > kneaddata/${sample_id}.clean_R2.fastq.gz; else : > kneaddata/${sample_id}.clean_R2.fastq && pigz -p ${threads} kneaddata/${sample_id}.clean_R2.fastq; fi
    c1=\$(ls -1 kneaddata/*_bowtie2_paired_contam_1.fastq 2>/dev/null | head -n1 || true)
    c2=\$(ls -1 kneaddata/*_bowtie2_paired_contam_2.fastq 2>/dev/null | head -n1 || true)
    if [ -n "\${c1}" ] && [ -f "\${c1}" ] && [ -n "\${c2}" ] && [ -f "\${c2}" ]; then
        pigz -p ${threads} -c "\${c1}" > kneaddata/${sample_id}.host_R1.fastq.gz
        pigz -p ${threads} -c "\${c2}" > kneaddata/${sample_id}.host_R2.fastq.gz
    else
        : > kneaddata/${sample_id}.host_R1.fastq && : > kneaddata/${sample_id}.host_R2.fastq
        pigz -p ${threads} kneaddata/${sample_id}.host_R1.fastq kneaddata/${sample_id}.host_R2.fastq
    fi
    """
}


