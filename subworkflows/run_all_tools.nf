nextflow.enable.dsl=2

include { RUN_BOWTIE2_DECONTAM    } from '../modules/bowtie2/main.nf'
include { BUILD_BOWTIE2_INDEX     } from '../modules/bowtie2/build.nf'
include { RUN_BWA_DECONTAM        } from '../modules/bwa/main.nf'
include { BUILD_BWA_INDEX         } from '../modules/bwa/build.nf'
include { RUN_KRAKEN2_DECONTAM    } from '../modules/kraken2/main.nf'
include { BUILD_KRAKEN2_DB        } from '../modules/kraken2/build.nf'
include { RUN_KRAKENUNIQ_DECONTAM } from '../modules/kraken_uniq/main.nf'
include { BUILD_KRAKENUNIQ_DB     } from '../modules/kraken_uniq/build.nf'
include { RUN_KNEADDATA_DECONTAM  } from '../modules/kneaddata/main.nf'
include { BUILD_KNEADDATA_DB      } from '../modules/kneaddata/build.nf'
include { RUN_KMCP_DECONTAM       } from '../modules/kmcp/main.nf'
include { BUILD_KMCP_DB           } from '../modules/kmcp/build.nf'

workflow RUN_ALL_TOOLS {
    take:
    reads_ch            // tuple(sample_id, R1, R2)
    host_ref_fasta      // path fasta
    host_taxid          // val taxid string or number
    kraken_taxonomy_dir // path taxonomy directory

    main:
    // Always build indices/DBs
    def bowtie2_idx_dir_ch   = params.run_bowtie2    ? BUILD_BOWTIE2_INDEX(file(host_ref_fasta)).index_dir : Channel.empty()
    def _bwa_idx_dir         = params.run_bwa        ? BUILD_BWA_INDEX(file(host_ref_fasta)).index_dir : Channel.empty()
    def bwa_ref_fasta_ch     = params.run_bwa        ? _bwa_idx_dir.map { dir -> file("${dir}/reference.fasta") } : Channel.empty()
    def k2_db_dir            = params.run_kraken2    ? BUILD_KRAKEN2_DB(file(host_ref_fasta), file(kraken_taxonomy_dir), host_taxid).db_dir : Channel.empty()
    def ku_db_dir            = params.run_krakenuniq ? BUILD_KRAKENUNIQ_DB(file(host_ref_fasta), file(kraken_taxonomy_dir), host_taxid).db_dir : Channel.empty()
    def kd_db_dir            = params.run_kneaddata  ? BUILD_KNEADDATA_DB(file(host_ref_fasta)).db_dir : Channel.empty()
    def kmcp_db              = params.run_kmcp       ? BUILD_KMCP_DB(file(host_ref_fasta)).db_dir : Channel.empty()

    // Fan-out to each tool using built resources
    def bowtie2_out = params.run_bowtie2    ? RUN_BOWTIE2_DECONTAM(reads_ch, bowtie2_idx_dir_ch) : [ cleaned: Channel.empty(), host: Channel.empty() ]
    def bwa_out     = params.run_bwa        ? RUN_BWA_DECONTAM(reads_ch, bwa_ref_fasta_ch) : [ cleaned: Channel.empty(), host: Channel.empty() ]
    def k2_out      = params.run_kraken2    ? RUN_KRAKEN2_DECONTAM(reads_ch, k2_db_dir, host_taxid) : [ cleaned: Channel.empty(), host: Channel.empty() ]
    def ku_out      = params.run_krakenuniq ? RUN_KRAKENUNIQ_DECONTAM(reads_ch, ku_db_dir, host_taxid) : [ cleaned: Channel.empty(), host: Channel.empty() ]
    def kd_out      = params.run_kneaddata  ? RUN_KNEADDATA_DECONTAM(reads_ch, kd_db_dir) : [ cleaned: Channel.empty(), host: Channel.empty() ]
    def kmcp_out    = params.run_kmcp       ? RUN_KMCP_DECONTAM(reads_ch, kmcp_db) : [ cleaned: Channel.empty(), host_ids: Channel.empty() ]

    // Also emit IDs for assessment subworkflow (host/clean per tool)
    // For FastQ outputs, emit the fastq.gz paths; assessment will extract IDs downstream
    def bowtie2_ids = bowtie2_out.host.combine(bowtie2_out.cleaned).map { sid, h1, h2, sid2, c1, c2 ->
        assert sid==sid2; tuple(sid, 'bowtie2', h1, h2, c1, c2)
    }
    def bwa_ids = bwa_out.host.combine(bwa_out.cleaned).map { sid, h1, h2, sid2, c1, c2 ->
        assert sid==sid2; tuple(sid, 'bwa', h1, h2, c1, c2)
    }
    def kraken2_ids = k2_out.host.combine(k2_out.cleaned).map { sid, h1, h2, sid2, c1, c2 ->
        assert sid==sid2; tuple(sid, 'kraken2', h1, h2, c1, c2)
    }
    def krakenuniq_ids = ku_out.host.combine(ku_out.cleaned).map { sid, h1, h2, sid2, c1, c2 ->
        assert sid==sid2; tuple(sid, 'krakenuniq', h1, h2, c1, c2)
    }
    def kneaddata_ids = kd_out.host.combine(kd_out.cleaned).map { sid, h1, h2, sid2, c1, c2 ->
        assert sid==sid2; tuple(sid, 'kneaddata', h1, h2, c1, c2)
    }
    def kmcp_ids = kmcp_out.host_ids.combine(kmcp_out.cleaned).map { sid, ids, sid2, c1, c2 ->
        assert sid==sid2; tuple(sid, 'kmcp', ids, ids, c1, c2)
    }

    def tool_ids = bowtie2_ids
        .mix(bwa_ids)
        .mix(kraken2_ids)
        .mix(krakenuniq_ids)
        .mix(kneaddata_ids)
        .mix(kmcp_ids)

    emit:
    bowtie2_cleaned = bowtie2_out.cleaned
    bwa_cleaned     = bwa_out.cleaned
    kraken2_cleaned = k2_out.cleaned
    krakenuniq_cleaned = ku_out.cleaned
    kneaddata_cleaned  = kd_out.cleaned
    kmcp_cleaned       = kmcp_out.cleaned
    tool_ids_for_assess = tool_ids
}


