nextflow.enable.dsl=2

process RESOURCE_SUMMARY {
    tag "resource_summary"
    cpus 1

    input:
    val trace_file
    val metrics_gate

    output:
    path("resources/tools_resources.tsv"), emit: table

    script:
    """
    set -euo pipefail
    mkdir -p resources
    { echo -e "process\ttag\trealtime\tpeak_rss"; [ -s "${trace_file}" ] && tail -n +2 "${trace_file}" | cut -f1-4 || true; } > resources/tools_resources.tsv
    """
}