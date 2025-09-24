nextflow.enable.dsl=2

process EXTRACT_IDS {
    tag "${sample_id}_${tool}_${label}"
    cpus 1
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/ids/${tool}/${file(fn).name}" } )

    input:
    tuple val(sample_id), val(tool), val(label), path(src)

    output:
    tuple val(sample_id), val(tool), val(label), path("ids/${tool}/${sample_id}.${label}.txt"), emit: ids

    script:
    def out = "ids/${tool}/${sample_id}.${label}.txt"
    """
    set -euo pipefail
    mkdir -p "\$(dirname "${out}")"
    name="${src}"
    if [[ "\$name" == *.fastq || "\$name" == *.fq || "\$name" == *.fastq.gz || "\$name" == *.fq.gz ]]; then
        if [[ "\$name" == *.gz ]]; then
            gzip -dc "${src}" | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); h=gensub("/[12]\$","",1,h); n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > "${out}"
        else
            awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); h=gensub("/[12]\$","",1,h); n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' "${src}" > "${out}"
        fi
    else
        cp "${src}" "${out}"
    fi
    """
}


