nextflow.enable.dsl=2

process RUN_BOWTIE2_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/bowtie2/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path index_dir

    output:
    tuple val(sample_id), path("bowtie2/${sample_id}.clean_R1.fastq.gz"), path("bowtie2/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("bowtie2/${sample_id}.host_R1.fastq.gz"),  path("bowtie2/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def sam_out = "bowtie2/${sample_id}.sam"
    def sorted_bam = "bowtie2/${sample_id}.sorted.bam"
    def aligned_bam = "bowtie2/${sample_id}.aligned.bam"
    def unaligned_bam = "bowtie2/${sample_id}.unaligned.bam"
    def aligned_fq1 = "bowtie2/${sample_id}.host_R1.fastq"
    def aligned_fq2 = "bowtie2/${sample_id}.host_R2.fastq"
    def unaligned_fq1 = "bowtie2/${sample_id}.clean_R1.fastq"
    def unaligned_fq2 = "bowtie2/${sample_id}.clean_R2.fastq"
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p bowtie2
    bowtie2 -x ${index_dir}/index -1 ${reads1} -2 ${reads2} -S ${sam_out} -p ${threads}
    samtools sort -@ ${threads} -o ${sorted_bam} ${sam_out}
    samtools index ${sorted_bam}
    samtools view -@ ${threads} -b -F 4 ${sorted_bam} > ${aligned_bam}
    samtools view -@ ${threads} -b -f 4 ${sorted_bam} > ${unaligned_bam}
    samtools fastq -@ ${threads} -1 ${aligned_fq1} -2 ${aligned_fq2} -0 /dev/null -s /dev/null -n ${aligned_bam}
    samtools fastq -@ ${threads} -1 ${unaligned_fq1} -2 ${unaligned_fq2} -0 /dev/null -s /dev/null -n ${unaligned_bam}
    pigz -p ${threads} ${aligned_fq1} ${aligned_fq2} ${unaligned_fq1} ${unaligned_fq2}
    """
}


