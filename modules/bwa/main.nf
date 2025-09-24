nextflow.enable.dsl=2

process RUN_BWA_DECONTAM {
    tag "${sample_id}"
    cpus 8
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/reads/bwa/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path bwa_index_dir

    output:
    tuple val(sample_id), path("bwa/${sample_id}.clean_R1.fastq.gz"), path("bwa/${sample_id}.clean_R2.fastq.gz"), emit: cleaned
    tuple val(sample_id), path("bwa/${sample_id}.host_R1.fastq.gz"),  path("bwa/${sample_id}.host_R2.fastq.gz"),  emit: host

    script:
    def sorted_bam = "bwa/${sample_id}.sorted.bam"
    def aligned_bam = "bwa/${sample_id}.aligned.bam"
    def aligned_fq1 = "bwa/${sample_id}.host_R1.fastq"
    def aligned_fq2 = "bwa/${sample_id}.host_R2.fastq"
    def host_ids = "bwa/${sample_id}.host_ids.txt"
    def threads = task.cpus
    def ref_fasta = "${bwa_index_dir}/reference.fasta"
    """
    set -euo pipefail
    mkdir -p bwa
    bwa mem -t ${threads} ${ref_fasta} ${reads1} ${reads2} | samtools sort -@ ${threads} -o ${sorted_bam}
    samtools view -@ ${threads} -b -F 4 ${sorted_bam} > ${aligned_bam}
    samtools fastq -@ ${threads} -1 ${aligned_fq1} -2 ${aligned_fq2} -0 /dev/null -s /dev/null -n ${aligned_bam}
    pigz -p ${threads} ${aligned_fq1} ${aligned_fq2}
    pigz -dc ${aligned_fq1}.gz | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > ${host_ids}
    pigz -dc ${aligned_fq2}.gz | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' >> ${host_ids}
    sort -u -o ${host_ids} ${host_ids}
    seqkit grep -v -f ${host_ids} ${reads1} | pigz -p ${threads} > bwa/${sample_id}.clean_R1.fastq.gz
    seqkit grep -v -f ${host_ids} ${reads2} | pigz -p ${threads} > bwa/${sample_id}.clean_R2.fastq.gz
    """
}


