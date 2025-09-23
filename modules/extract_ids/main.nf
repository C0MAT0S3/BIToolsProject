nextflow.enable.dsl=2

process EXTRACT_IDS {
    tag "${sample_id}_${tool}_${label}"
    cpus 1

    input:
    tuple val(sample_id), val(tool), val(label), path(src)

    output:
    tuple val(sample_id), val(tool), val(label), path("ids/${tool}_${sample_id}_${label}.txt"), emit: ids

    script:
    def out = "ids/${tool}_${sample_id}_${label}.txt"
    """
    set -euo pipefail
    mkdir -p ids
    name="${src}"
    if [[ "\$name" == *.fastq || "\$name" == *.fq || "\$name" == *.fastq.gz || "\$name" == *.fq.gz ]]; then
        if [[ "\$name" == *.gz ]]; then
            gzip -dc "${src}" | awk 'NR%4==1{gsub(/^@/,""); print \$1}' > "${out}"
        else
            awk 'NR%4==1{gsub(/^@/,""); print \$1}' "${src}" > "${out}"
        fi
    else
        # Assume already an ID list; normalize by copying
        cp "${src}" "${out}"
    fi
    """
}


