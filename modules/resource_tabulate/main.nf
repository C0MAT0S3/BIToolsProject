nextflow.enable.dsl=2

process RESOURCE_TABULATE {
    tag "resource_tabulate"
    cpus 1

    input:
    val trace_file

    output:
    path("resources/tools_resources.tsv"), emit: table

    script:
    """
    set -euo pipefail
    mkdir -p resources

    # Header: process\ttool\ttag\trealtime_s\tpeak_rss_mb
    echo -e "process\ttool\ttag\trealtime_s\tpeak_rss_mb" > resources/tools_resources.tsv

    awk -F '\t' 'NR>1 {
        proc=\$1; tag=\$2; rt=\$3; mem=\$4;
        tool="unknown";
        if (proc ~ /RUN_BOWTIE2_DECONTAM/) tool="bowtie2";
        else if (proc ~ /RUN_BWA_DECONTAM/) tool="bwa";
        else if (proc ~ /RUN_KRAKEN2_DECONTAM/) tool="kraken2";
        else if (proc ~ /RUN_KRAKENUNIQ_DECONTAM/) tool="krakenuniq";
        else if (proc ~ /RUN_KNEADDATA_DECONTAM/) tool="kneaddata";
        else if (proc ~ /RUN_KMCP_DECONTAM/) tool="kmcp";
        # Convert realtime to seconds: try HH:MM:SS(.ms) or numeric already
        secs=rt;
        if (rt ~ /:/) {
            n=split(rt, a, ":");
            if (n==3) {
                split(a[3], b, ".");
                secs = a[1]*3600 + a[2]*60 + b[1];
            }
        }
        # Convert peak_rss (bytes) to MB if numeric
        mb=mem;
        if (mem ~ /^[-]?[0-9]+(\.[0-9]+)?$/) { mb = mem/1024/1024; }
        printf("%s\t%s\t%s\t%.2f\t%.2f\n", proc, tool, tag, secs, mb) >> "resources/tools_resources.tsv";
    }' "${trace_file}"
    """
}


