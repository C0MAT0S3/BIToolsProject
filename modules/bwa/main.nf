nextflow.enable.dsl=2

process RUN_BWA_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/bwa/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path ref_fasta

    output:
    tuple val(sample_id), path("bwa/${sample_id}.clean_R1.fastq.gz"), path("bwa/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("bwa/${sample_id}.host_R1.fastq.gz"),  path("bwa/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def sorted_bam = "bwa/${sample_id}.sorted.bam"
    def aligned_bam = "bwa/${sample_id}.aligned.bam"
    def unaligned_bam = "bwa/${sample_id}.unaligned.bam"
    def aligned_fq1 = "bwa/${sample_id}.host_R1.fastq"
    def aligned_fq2 = "bwa/${sample_id}.host_R2.fastq"
    def unaligned_fq1 = "bwa/${sample_id}.clean_R1.fastq"
    def unaligned_fq2 = "bwa/${sample_id}.clean_R2.fastq"
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p bwa
    bwa mem -t ${threads} ${ref_fasta} ${reads1} ${reads2} | samtools sort -@ ${threads} -o ${sorted_bam}
    samtools index ${sorted_bam}
    samtools view -@ ${threads} -b -F 4 ${sorted_bam} > ${aligned_bam}
    samtools view -@ ${threads} -b -f 4 ${sorted_bam} > ${unaligned_bam}
    samtools fastq -@ ${threads} -1 ${aligned_fq1} -2 ${aligned_fq2} -0 /dev/null -s /dev/null -n ${aligned_bam}
    samtools fastq -@ ${threads} -1 ${unaligned_fq1} -2 ${unaligned_fq2} -0 /dev/null -s /dev/null -n ${unaligned_bam}
    pigz -p ${threads} ${aligned_fq1} ${aligned_fq2} ${unaligned_fq1} ${unaligned_fq2}
    """
}


