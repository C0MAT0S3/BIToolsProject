nextflow.enable.dsl=2

process CALC_CLASSIFY_METRICS {
    tag "${sample_id}_${tool}"

    input:
    tuple val(sample_id), val(tool), path(metrics_dir)

    output:
    tuple val(sample_id), val(tool), path("metrics/${tool}.tsv"), emit: metrics

    script:
    def out = "metrics/${tool}.tsv"
    """
    set -euo pipefail
    mkdir -p metrics
    tp_r1=\$([ -f "${metrics_dir}/TP_R1.txt" ] && wc -l < "${metrics_dir}/TP_R1.txt" || echo 0)
    tp_r2=\$([ -f "${metrics_dir}/TP_R2.txt" ] && wc -l < "${metrics_dir}/TP_R2.txt" || echo 0)
    fp_r1=\$([ -f "${metrics_dir}/FP_R1.txt" ] && wc -l < "${metrics_dir}/FP_R1.txt" || echo 0)
    fp_r2=\$([ -f "${metrics_dir}/FP_R2.txt" ] && wc -l < "${metrics_dir}/FP_R2.txt" || echo 0)
    fn_r1=\$([ -f "${metrics_dir}/FN_R1.txt" ] && wc -l < "${metrics_dir}/FN_R1.txt" || echo 0)
    fn_r2=\$([ -f "${metrics_dir}/FN_R2.txt" ] && wc -l < "${metrics_dir}/FN_R2.txt" || echo 0)
    tn_r1=\$([ -f "${metrics_dir}/TN_R1.txt" ] && wc -l < "${metrics_dir}/TN_R1.txt" || echo 0)
    tn_r2=\$([ -f "${metrics_dir}/TN_R2.txt" ] && wc -l < "${metrics_dir}/TN_R2.txt" || echo 0)
    total_tp=\$((tp_r1 + tp_r2))
    total_fp=\$((fp_r1 + fp_r2))
    total_fn=\$((fn_r1 + fn_r2))
    total_tn=\$((tn_r1 + tn_r2))
    printf "sample_id\ttool\tTP\tFP\tFN\tTN\tprecision\trecall\taccuracy\tf1\n" > "${out}"
    awk -v sid='${sample_id}' -v tool='${tool}' \
        -v TP="\${total_tp}" -v FP="\${total_fp}" -v FN="\${total_fn}" -v TN="\${total_tn}" 'BEGIN {
            tp = TP + 0; fp = FP + 0; fn = FN + 0; tn = TN + 0;
            pr_d = (tp + fp); re_d = (tp + fn); acc_d = (tp + fp + fn + tn);
            precision = (pr_d>0) ? tp / pr_d : 0;
            recall    = (re_d>0) ? tp / re_d : 0;
            f1_d = (2*tp + fp + fn);
            f1 = (f1_d>0) ? (2.0 * tp) / f1_d : 0;
            accuracy  = (acc_d>0) ? (tp + tn) / acc_d : 0;
            OFS = "\t";
            print sid, tool, tp, fp, fn, tn, precision, recall, accuracy, f1;
        }' >> "${out}"
    """
}


process MERGE_SAMPLE_METRICS {
    tag "${sample_id}"
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/metrics/${file(fn).name}" } )

    input:
    tuple val(sample_id), path(tool_tsvs)

    output:
    tuple val(sample_id), path("metrics/summary.tsv"), emit: summaries

    script:
    def out = "metrics/summary.tsv"
    """
    set -euo pipefail
    mkdir -p metrics
    printf "sample_id\ttool\tTP\tFP\tFN\tTN\tprecision\trecall\taccuracy\tf1\\n" > "${out}"
    TMP_FILE=\$(mktemp)
    for f in ${tool_tsvs}; do
        tail -n +2 "\$f" >> "\$TMP_FILE"
    done
    LC_ALL=C sort -t\$'\t' -k2,2 -k1,1 "\$TMP_FILE" >> "${out}"
    """
}


