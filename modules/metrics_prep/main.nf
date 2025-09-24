nextflow.enable.dsl=2

process METRICS_PREP {
    tag "${sample_id}_${tool}"
    cpus 1
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/metrics/${file(fn).name}" } )

    input:
    tuple val(sample_id), val(tool), path(truth_dir), path(host_ids_r1), path(host_ids_r2), path(clean_ids_r1), path(clean_ids_r2)

    output:
    tuple val(sample_id), val(tool), path("${tool}"), emit: f1dir
    path("${tool}/*.txt")

    script:
    def outd = "${tool}"
    """
    set -euo pipefail
    mkdir -p ${outd}
    cp ${truth_dir}/HOST_R1.ids ${outd}/HOST_R1.txt
    cp ${truth_dir}/HOST_R2.ids ${outd}/HOST_R2.txt
    cp ${truth_dir}/MICROBIOME_R1.ids ${outd}/MICROBIOME_R1.txt
    cp ${truth_dir}/MICROBIOME_R2.ids ${outd}/MICROBIOME_R2.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${host_ids_r1} ${outd}/HOST_R1.txt > ${outd}/TP_R1.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${host_ids_r2} ${outd}/HOST_R2.txt > ${outd}/TP_R2.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${host_ids_r1} ${outd}/MICROBIOME_R1.txt > ${outd}/FP_R1.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${host_ids_r2} ${outd}/MICROBIOME_R2.txt > ${outd}/FP_R2.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${clean_ids_r1} ${outd}/HOST_R1.txt > ${outd}/FN_R1.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${clean_ids_r2} ${outd}/HOST_R2.txt > ${outd}/FN_R2.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${clean_ids_r1} ${outd}/MICROBIOME_R1.txt > ${outd}/TN_R1.txt
    awk 'FNR==NR{a[\$1]; next} \$1 in a' ${clean_ids_r2} ${outd}/MICROBIOME_R2.txt > ${outd}/TN_R2.txt
    """
}


