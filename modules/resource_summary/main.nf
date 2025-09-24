nextflow.enable.dsl=2

process RESOURCE_SUMMARY {
    tag "resource_summary"
    cpus 1
    publishDir( params.output_dir ?: "results", mode: 'copy', overwrite: true, saveAs: { fn -> "${file(fn).name}" } )

    input:
    path trace_file
    val metrics_gate

    output:
    path("resources_by_sample.tsv"), emit: by_sample

    script:
    """
    set -euo pipefail
    :
    {
      echo -e "sample_id\ttool\ttotal_realtime_hms\tmax_peak_rss_GB";
      if [ -s "${trace_file}" ]; then
        awk 'BEGIN{FS=OFS="\t"}
          NR==1{for(i=1;i<=NF;i++) idx[\$i]=i; next}
          {
            proc = (idx["process"] ? \$(idx["process"]) : "");
            tag  = (idx["tag"] ? \$(idx["tag"]) : "");
            s_ms = (idx["realtime"] ? \$(idx["realtime"]) + 0 : 0);
            s    = s_ms / 1000.0;
            b    = (idx["peak_rss"] ? \$(idx["peak_rss"]) + 0 : 0);
            cpu  = (idx["%cpu"] ? \$(idx["%cpu"]) + 0 : 0);
            # Include both RUN_* and BUILD_* processes
            if (proc ~ /^RUN_ALL_TOOLS:/) {
              step=substr(proc, index(proc, ":")+1);
              tool=step;
              sub(/^BUILD_|^RUN_/ , "", tool);
              sub(/_.*/, "", tool);
              tool_l=tolower(tool);
              key=tag "\t" tool_l;
              cnt[key]++;
              sum_s[key]+=s;
              sum_cpu[key]+=cpu;
              if (b>maxb[key]) maxb[key]=b;
              if (pvm>maxvm[key]) maxvm[key]=pvm;
            }
          }
          END{
            for (k in cnt){
              s=sum_s[k]+0; n=cnt[k]+0; h=int(s/3600); m=int((s%3600)/60); sec=int(s%60);
              hms=sprintf("%02d:%02d:%02d",h,m,sec);
              gb=maxb[k]/1024/1024/1024;
              printf "%s\t%s\t%.3f\\n", k, hms, gb;
            }
          }' "${trace_file}" | LC_ALL=C sort -t\$'\t' -k2,2 -k1,1
      fi
      } > resources_by_sample.tsv
    """
}