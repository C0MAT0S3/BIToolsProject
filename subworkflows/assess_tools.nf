nextflow.enable.dsl=2

include { PREPARE_TRUTH         } from '../modules/prepare_truth/main.nf'
include { METRICS_PREP          } from '../modules/metrics_prep/main.nf'
include { CALC_CLASSIFY_METRICS } from '../modules/metrics/main.nf'
include { MERGE_SAMPLE_METRICS  } from '../modules/metrics/main.nf'
include { EXTRACT_IDS           } from '../modules/extract_ids/main.nf'

workflow ASSESS_TOOLS {
    take:
    samples_ch
    tool_ids_ch

    main:
    def truth_ch = PREPARE_TRUTH(samples_ch)

    def keyed_truth = truth_ch.truth_dir.map { sid, tdir -> tuple(sid, tdir) }
    def keyed_tool  = tool_ids_ch.map { sid, tool, h1, h2, c1, c2 -> tuple(sid, tool, h1, h2, c1, c2) }

    def to_labels = keyed_tool.flatMap { sid, tool, h1, h2, c1, c2 ->
        [
            tuple(sid, tool, 'host_R1', h1),
            tuple(sid, tool, 'host_R2', h2),
            tuple(sid, tool, 'clean_R1', c1),
            tuple(sid, tool, 'clean_R2', c2)
        ]
    }

    def extracted = EXTRACT_IDS(to_labels)
        .ids
        .map { sid, tool, label, idsfile -> tuple(sid, tool, label, idsfile) }

    def _host_r1 = extracted.filter { _sid, _tool, label, _f -> label == 'host_R1' }.map { sid, tool, _label, f -> tuple(sid, tool, f) }
    def _host_r2 = extracted.filter { _sid, _tool, label, _f -> label == 'host_R2' }.map { sid, tool, _label, f -> tuple(sid, tool, f) }
    def _clean_r1 = extracted.filter { _sid, _tool, label, _f -> label == 'clean_R1' }.map { sid, tool, _label, f -> tuple(sid, tool, f) }
    def _clean_r2 = extracted.filter { _sid, _tool, label, _f -> label == 'clean_R2' }.map { sid, tool, _label, f -> tuple(sid, tool, f) }

    def merged_ids = extracted
        .groupTuple(by: [0,1])
        .map { sid, tool, labels, files ->
            def mapv = [:]
            labels.eachWithIndex { l, i -> mapv[l] = files[i] }
            def nullf = file('/dev/null')
            tuple(sid, tool,
                  mapv['host_R1'] ?: nullf,
                  mapv['host_R2'] ?: nullf,
                  mapv['clean_R1'] ?: nullf,
                  mapv['clean_R2'] ?: nullf)
        }

    def joined = merged_ids.combine(keyed_truth, by:0)
        .map { sid, tool, h1, h2, c1, c2, tdir -> tuple(sid, tool, tdir, h1, h2, c1, c2) }

    def f1_dirs = METRICS_PREP(joined)
        .f1dir
        .map { sid, tool, dir -> tuple(sid, tool, dir) }

    def metrics = CALC_CLASSIFY_METRICS(
            f1_dirs.map { sid, tool, dir -> tuple(sid, tool, dir) }
        )

    def per_sample = metrics.metrics
        .groupTuple(by: 0)
        .map { sid, _tools, files -> tuple(sid, files) }

    def summaries = MERGE_SAMPLE_METRICS(per_sample)

    emit:
    metrics_tsv = summaries.summaries
}


