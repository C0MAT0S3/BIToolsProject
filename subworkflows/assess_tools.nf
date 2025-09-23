nextflow.enable.dsl=2

include { PREPARE_TRUTH         } from '../modules/prepare_truth/main.nf'
include { METRICS_PREP          } from '../modules/metrics_prep/main.nf'
include { CALC_CLASSIFY_METRICS } from '../modules/metrics/main.nf'
include { EXTRACT_IDS           } from '../modules/extract_ids/main.nf'

workflow ASSESS_TOOLS {
    take:
    samples_ch      // tuple(sample_id, raw_r1, raw_r2, neg_r1, neg_r2)
    tool_ids_ch     // tuples per tool with host/clean IDs per read: see below

    /*
      tool_ids_ch emits per sample and tool:
        tuple(sample_id, tool, host_r1_fq_or_ids, host_r2_fq_or_ids, clean_r1_fq_or_ids, clean_r2_fq_or_ids)
      This flow will extract IDs from fastq(.gz) if necessary.
    */

    main:
    def truth_ch = PREPARE_TRUTH(samples_ch)

    // join truth with tool ids on sample_id
    def keyed_truth = truth_ch.truth_dir.map { sid, tdir -> tuple(sid, tdir) }
    def keyed_tool  = tool_ids_ch.map { sid, tool, h1, h2, c1, c2 -> tuple(sid, tool, h1, h2, c1, c2) }

    // Normalize to ID lists using a single module that handles FastQ or ID lists
    def to_labels = keyed_tool.flatMap { sid, tool, h1, h2, c1, c2 ->
        Channel.of(
            tuple(sid, tool, 'host_R1', h1),
            tuple(sid, tool, 'host_R2', h2),
            tuple(sid, tool, 'clean_R1', c1),
            tuple(sid, tool, 'clean_R2', c2)
        )
    }

    def extracted = EXTRACT_IDS(to_labels)
        .ids
        .map { sid, tool, label, idsfile -> tuple(sid, tool, label, idsfile) }

    def merged_ids = extracted
        .groupTuple(by: [0,1])
        .map { key, rows ->
            def sid = key[0]; def tool = key[1]
            def mapv = rows.collectEntries { [ (it[2]) : it[3] ] }
            tuple(sid, tool, mapv['host_R1'], mapv['host_R2'], mapv['clean_R1'], mapv['clean_R2'])
        }

    def joined = merged_ids.join(keyed_truth, by:0).map { sid, tool, h1, h2, c1, c2, tdir ->
        tuple(sid, tool, tdir, h1, h2, c1, c2)
    }

    def f1_dirs = METRICS_PREP(joined)
        .f1dir
        .map { sid, tool, dir -> tuple(sid, tool, dir) }

    def metrics = CALC_CLASSIFY_METRICS(
            f1_dirs.map { sid, tool, dir -> tuple(sid, tool, dir) }
        )

    emit:
    metrics_tsv = metrics.metrics
}


