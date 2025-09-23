nextflow.enable.dsl=2

process CALC_CLASSIFY_METRICS {
    tag "${sample_id}_${tool}"
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${sample_id}/metrics/${file(fn).name}" } )

    input:
    tuple val(sample_id), val(tool), path(metrics_dir)

    output:
    tuple val(sample_id), val(tool), path("metrics/${tool}_${sample_id}.tsv"), emit: metrics

    script:
    def out = "metrics/${tool}_${sample_id}.tsv"
    """
    set -euo pipefail
    mkdir -p metrics

    # count lines or 0 if file missing
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

    # compute metrics with zero-division guards
    awk -v sid='${sample_id}' -v tool='${tool}' \
        -v TP=\"\${total_tp}\" -v FP=\"\${total_fp}\" -v FN=\"\${total_fn}\" -v TN=\"\${total_tn}\" 'BEGIN {
            tp = TP + 0; fp = FP + 0; fn = FN + 0; tn = TN + 0;
            pr_d = (tp + fp); re_d = (tp + fn); acc_d = (tp + fp + fn + tn);
            precision = (pr_d>0) ? tp / pr_d : 0;
            recall    = (re_d>0) ? tp / re_d : 0;
            f1_d = (2*tp + fp + fn);
            f1 = (f1_d>0) ? (2.0 * tp) / f1_d : 0;
            accuracy  = (acc_d>0) ? (tp + tn) / acc_d : 0;
            printf("sample_id\ttool\tTP\tFP\tFN\tTN\tprecision\trecall\taccuracy\tf1\n") > "${out}";
            printf("%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n", sid, tool, tp, fp, fn, tn, precision, recall, accuracy, f1) >> "${out}";
        }'
    """
}


