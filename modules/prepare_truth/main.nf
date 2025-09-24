nextflow.enable.dsl=2

process PREPARE_TRUTH {
    tag "${sample_id}"
    cpus 2

    input:
    tuple val(sample_id), path(raw_r1), path(raw_r2), path(neg_r1), path(neg_r2)

    output:
    tuple val(sample_id), path("ids"), emit: truth_dir

    script:
    def outd = "ids"
    """
    set -euo pipefail
    mkdir -p ${outd}
    if [[ "${raw_r1}" == *.gz ]]; then gzip -dc "${raw_r1}" | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > ${outd}/RAW_R1.ids; else awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' "${raw_r1}" > ${outd}/RAW_R1.ids; fi
    if [[ "${raw_r2}" == *.gz ]]; then gzip -dc "${raw_r2}" | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > ${outd}/RAW_R2.ids; else awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' "${raw_r2}" > ${outd}/RAW_R2.ids; fi
    if [[ "${neg_r1}" == *.gz ]]; then gzip -dc "${neg_r1}" | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > ${outd}/MICROBIOME_R1.ids; else awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' "${neg_r1}" > ${outd}/MICROBIOME_R1.ids; fi
    if [[ "${neg_r2}" == *.gz ]]; then gzip -dc "${neg_r2}" | awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' > ${outd}/MICROBIOME_R2.ids; else awk 'NR%4==1{h=\$1; sub(/^@/,"",h); sub(/ .*/,"",h); sub(/#.*/,"",h); if(length(h)>=2){t=substr(h,length(h)-1); if(t=="/1"||t=="/2") h=substr(h,1,length(h)-2)}; n=split(h,a,"."); if(n>=2) print a[1] "." a[2]; else print h}' "${neg_r2}" > ${outd}/MICROBIOME_R2.ids; fi
    awk 'FNR==NR{a[\$1]; next} !(\$1 in a)' ${outd}/MICROBIOME_R1.ids ${outd}/RAW_R1.ids > ${outd}/HOST_R1.ids
    awk 'FNR==NR{a[\$1]; next} !(\$1 in a)' ${outd}/MICROBIOME_R2.ids ${outd}/RAW_R2.ids > ${outd}/HOST_R2.ids
    """
}


