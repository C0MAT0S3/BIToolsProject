nextflow.enable.dsl=2

process PREPARE_TRUTH {
    tag "${sample_id}"
    cpus 2
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/truth/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(raw_r1), path(raw_r2), path(neg_r1), path(neg_r2)

    output:
    tuple val(sample_id), path("truth/${sample_id}"), emit: truth_dir

    script:
    def outd = "truth/${sample_id}"
    """
    set -euo pipefail
    mkdir -p ${outd}

    # extract read IDs (without @) from FASTQ; handle gz automatically
    if [[ "${raw_r1}" == *.gz ]]; then gzip -dc "${raw_r1}" | awk 'NR%4==1{gsub(/^@/,""); print \$1}' > ${outd}/RAW_R1.ids; else awk 'NR%4==1{gsub(/^@/,""); print \$1}' "${raw_r1}" > ${outd}/RAW_R1.ids; fi
    if [[ "${raw_r2}" == *.gz ]]; then gzip -dc "${raw_r2}" | awk 'NR%4==1{gsub(/^@/,""); print \$1}' > ${outd}/RAW_R2.ids; else awk 'NR%4==1{gsub(/^@/,""); print \$1}' "${raw_r2}" > ${outd}/RAW_R2.ids; fi
    if [[ "${neg_r1}" == *.gz ]]; then gzip -dc "${neg_r1}" | awk 'NR%4==1{gsub(/^@/,""); print \$1}' > ${outd}/MICROBIOME_R1.ids; else awk 'NR%4==1{gsub(/^@/,""); print \$1}' "${neg_r1}" > ${outd}/MICROBIOME_R1.ids; fi
    if [[ "${neg_r2}" == *.gz ]]; then gzip -dc "${neg_r2}" | awk 'NR%4==1{gsub(/^@/,""); print \$1}' > ${outd}/MICROBIOME_R2.ids; else awk 'NR%4==1{gsub(/^@/,""); print \$1}' "${neg_r2}" > ${outd}/MICROBIOME_R2.ids; fi

    # host = raw - microbiome
    awk 'FNR==NR{a[\$1]; next} !(\$1 in a)' ${outd}/MICROBIOME_R1.ids ${outd}/RAW_R1.ids > ${outd}/HOST_R1.ids
    awk 'FNR==NR{a[\$1]; next} !(\$1 in a)' ${outd}/MICROBIOME_R2.ids ${outd}/RAW_R2.ids > ${outd}/HOST_R2.ids
    """
}


