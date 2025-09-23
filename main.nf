nextflow.enable.dsl=2

include { RUN_ALL_TOOLS    } from './subworkflows/run_all_tools.nf'
include { ASSESS_TOOLS     } from './subworkflows/assess_tools.nf'
include { RESOURCE_SUMMARY } from './modules/resource_summary/main.nf'

workflow {
    // samplename,R1,R2,neg_R1,neg_R2
    if (!params.input) {
        error "Missing required parameter: --input samplesheet.csv (columns: samplename,R1,R2)"
    }

    // Required resources for builds and runs
    def host_ref_fasta     = params.bwa_ref ?: params.host_ref ?: null
    def kraken_tax_dir     = params.kraken_taxonomy_dir ?: null
    def host_taxid         = params.host_taxid ?: '9606' // Human

    log.info """\
    HOST-DECONTAMINATION COMPARISON PIPELINE
    =======================================
    Input samplesheet               : ${params.input}
    Output directory                : ${params.output_dir ?: 'results'}
    Host reference FASTA            : ${host_ref_fasta}
    Kraken taxonomy dir             : ${kraken_tax_dir}
    Host taxid                      : ${host_taxid}
    """.stripIndent()

    reads_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true)
        .map { row ->
            def sid = row.samplename?.toString()?.trim()
            def r1  = file(row.R1)
            def r2  = file(row.R2)
            def nr1 = file(row.neg_R1)
            def nr2 = file(row.neg_R2)
            if (!sid || !r1 || !r2 || !nr1 || !nr2 || !host_ref_fasta || !kraken_tax_dir) {
                error "Invalid samplesheet row (need samplename,R1,R2,neg_R1,neg_R2): ${row}"
            }
            tuple(sid, r1, r2)
        }

    def all_tools = RUN_ALL_TOOLS(
        reads_ch,
        file(host_ref_fasta),
        host_taxid,
        file(kraken_tax_dir)
    )

    // Build samples channel for assessment with negative controls
    samples_for_assess = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true)
        .map { row ->
            tuple(row.samplename.toString().trim(), file(row.R1), file(row.R2), file(row.neg_R1), file(row.neg_R2))
        }

    def metrics_out = ASSESS_TOOLS(
        samples_for_assess,
        all_tools.tool_ids_for_assess
    )

    // Post-run summary
    RESOURCE_SUMMARY("${params.output_dir ?: 'results'}/trace.tsv", metrics_out.metrics_tsv.map { true }.first())
}