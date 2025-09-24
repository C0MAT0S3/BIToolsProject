nextflow.enable.dsl=2

process RUN_BOWTIE2_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/reads/bowtie2/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path index_dir

    output:
    tuple val(sample_id), path("bowtie2/${sample_id}.clean_R1.fastq.gz"), path("bowtie2/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("bowtie2/${sample_id}.host_R1.fastq.gz"),  path("bowtie2/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def sorted_bam = "bowtie2/${sample_id}.sorted.bam"
    def aligned_bam = "bowtie2/${sample_id}.aligned.bam"
    def aligned_fq1 = "bowtie2/${sample_id}.host_R1.fastq"
    def aligned_fq2 = "bowtie2/${sample_id}.host_R2.fastq"
    def host_ids = "bowtie2/${sample_id}.host_ids.txt"
    def threads = task.cpus
    """
    set -euo pipefail
    mkdir -p bowtie2
    bowtie2 -x ${index_dir}/index -1 ${reads1} -2 ${reads2} -p ${threads} | samtools sort -@ ${threads} -o ${sorted_bam}
    samtools view -@ ${threads} -b -F 4 ${sorted_bam} > ${aligned_bam}
    samtools fastq -@ ${threads} -1 ${aligned_fq1} -2 ${aligned_fq2} -0 /dev/null -s /dev/null -n ${aligned_bam}
    pigz -p ${threads} ${aligned_fq1} ${aligned_fq2}
    pigz -dc ${aligned_fq1}.gz | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > ${host_ids}
    pigz -dc ${aligned_fq2}.gz | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' >> ${host_ids}
    sort -u -o ${host_ids} ${host_ids}
    seqkit grep -v -f ${host_ids} ${reads1} | pigz -p ${threads} > bowtie2/${sample_id}.clean_R1.fastq.gz
    seqkit grep -v -f ${host_ids} ${reads2} | pigz -p ${threads} > bowtie2/${sample_id}.clean_R2.fastq.gz
    """
}


