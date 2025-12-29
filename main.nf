#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CastRNA2 - Consensus Projection Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Derives RNA secondary structures by projecting Rfam consensus structures
    onto query sequences via covariance model alignment.
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
    log.info """
    =========================================
     CastRNA2 v${workflow.manifest.version}
    =========================================

    Usage:
      nextflow run main.nf --input <input_file> [options]

    Mandatory arguments:
      --input               Path to input file. Can be:
                            - FASTA file: all sequences will be scanned against Rfam
                            - CSV samplesheet: columns 'seq_id', 'fasta', 'rfam_id'

    Optional arguments:
      --mapping             Path to mapping CSV (seq_id,rfam_id) for use with multi-FASTA
      --outdir              Output directory (default: 'results')
      --rfam_cm             Path to Rfam.cm or Rfam.cm.gz (optional, downloads if not provided)
      --rfam_clanin         Path to Rfam.clanin (optional, not required for pipeline)
      --rfam_seed           Path to Rfam.seed or Rfam.seed.gz (for seed alignment method)
      --rfam_release        Rfam release version, e.g., '14.9' (default: CURRENT)
                            Note: older Rfam versions may not include the clanin file

    Alignment method:
      --alignment_method    Method for aligning sequences to Rfam families:
                            'cmalign' (default): Use Infernal cmalign to compute alignments
                            'seed': Use pre-computed Rfam seed alignments (like bpRNA)
                            
                            Note: 'seed' method only works for sequences that are already
                            in the Rfam seed alignment. Use this to reproduce bpRNA results.
                            The 'cmalign' method works for any sequence but may produce
                            slightly different alignments.

    Thresholds:
      --cmscan_evalue       E-value threshold for cmscan (default: 1e-5)
      --cmscan_score        Bit score threshold for cmscan (default: 20.0)

    Projection options:
      --noncanonical_pairs  How to handle non-canonical base pairs:
                            'keep' (default), 'remove', or 'flag'

    MSA output:
      --output_msa          Generate multi-sequence alignments grouped by Rfam family
                            Outputs Stockholm (.sto) and aligned FASTA (.afa) per family

    Examples:
      # Inference mode (all sequences scanned against Rfam)
      nextflow run main.nf --input sequences.fasta

      # With mapping file (known Rfam IDs)
      nextflow run main.nf --input sequences.fasta --mapping mapping.csv

      # Use specific Rfam version
      nextflow run main.nf --input sequences.fasta --rfam_release 14.9

      # Reproduce bpRNA results using seed alignments
      nextflow run main.nf --input sequences.fasta --alignment_method seed --rfam_release 12.2

      # Generate MSA outputs grouped by Rfam family
      nextflow run main.nf --input sequences.fasta --mapping mapping.csv --output_msa

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.input) {
    error "ERROR: Please provide an input file with --input"
}

log.info """
=========================================
 CastRNA2 v${workflow.manifest.version}
=========================================
Input       : ${params.input}
Mapping     : ${params.mapping ?: 'none (inference mode)'}
Output dir  : ${params.outdir}
Rfam release: ${params.rfam_release ?: 'CURRENT'}
Alignment   : ${params.alignment_method ?: 'cmalign'}
Output MSA  : ${params.output_msa}
cmscan E    : ${params.cmscan_evalue}
cmscan T    : ${params.cmscan_score}
-----------------------------------------
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_RFAM         } from './subworkflows/local/prepare_rfam'
include { INFERNAL_CMSCAN      } from './modules/local/infernal_cmscan/main'
include { INFERNAL_CMALIGN     } from './modules/local/infernal_cmalign/main'
include { SEED_ALIGNMENT       } from './modules/local/seed_alignment/main'
include { CONSENSUS_PROJECTION } from './modules/local/consensus_projection/main'
include { GENERATE_SUMMARY     } from './modules/local/generate_summary/main'
include { CMALIGN_GROUP        } from './modules/local/cmalign_group/main'
include { GENERATE_MSA         } from './modules/local/generate_msa/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL PROCESSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process WRITE_FASTA {
    // Write a sequence string to a FASTA file
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), val(sequence)

    output:
    tuple val(meta), path("*.fa")

    script:
    """
    echo ">${meta.id}" > ${meta.id}.fa
    echo "${sequence}" >> ${meta.id}.fa
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    ch_versions = Channel.empty()

    //
    // STEP 1: Prepare Rfam database
    //
    PREPARE_RFAM()
    ch_versions = ch_versions.mix(PREPARE_RFAM.out.versions)

    //
    // STEP 2: Parse input files and create channels
    //
    def input_path = file(params.input)
    def is_csv = input_path.name.endsWith('.csv')

    if (is_csv) {
        // CSV samplesheet mode: seq_id, fasta, rfam_id
        Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                def meta = [
                    id: row.seq_id,
                    rfam_id: row.rfam_id?.trim() ?: null
                ]
                def fasta = file(row.fasta, checkIfExists: true)
                [meta, fasta]
            }
            .set { ch_input }

    } else if (params.mapping) {
        // FASTA + mapping CSV mode
        // Load mapping into a map
        Channel
            .fromPath(params.mapping)
            .splitCsv(header: true)
            .map { row -> [row.seq_id, row.rfam_id?.trim() ?: null] }
            .set { ch_mapping }

        // Split FASTA and join with mapping
        Channel
            .fromPath(params.input)
            .splitFasta(record: [id: true, seqString: true])
            .map { record -> [record.id, record.seqString] }
            .join(ch_mapping, by: 0, remainder: true)
            .map { seq_id, sequence, rfam_id ->
                def meta = [id: seq_id, rfam_id: rfam_id]
                [meta, sequence]
            }
            .set { ch_sequences_with_mapping }

        // Write sequences to FASTA files
        WRITE_FASTA(ch_sequences_with_mapping)
        ch_input = WRITE_FASTA.out

    } else {
        // Pure FASTA mode - all sequences go through inference
        Channel
            .fromPath(params.input)
            .splitFasta(record: [id: true, seqString: true])
            .map { record ->
                def meta = [id: record.id, rfam_id: null]
                [meta, record.seqString]
            }
            .set { ch_sequences }

        WRITE_FASTA(ch_sequences)
        ch_input = WRITE_FASTA.out
    }

    //
    // STEP 3: Branch into known vs unknown Rfam IDs
    //
    ch_input
        .branch {
            meta, fasta ->
                known: meta.rfam_id != null
                unknown: true
        }
        .set { ch_branched }

    //
    // STEP 4: Run cmscan on unknown sequences
    //
    INFERNAL_CMSCAN(
        ch_branched.unknown,
        PREPARE_RFAM.out.cm_index.collect()
    )
    ch_versions = ch_versions.mix(INFERNAL_CMSCAN.out.versions.first().ifEmpty([]))

    //
    // STEP 5: Parse cmscan results and select best hits
    //
    INFERNAL_CMSCAN.out.tblout
        .join(ch_branched.unknown)
        .map { meta, tblout, fasta ->
            // Parse tblout to find best hit
            def hits = []
            tblout.eachLine { line ->
                if (!line.startsWith('#') && line.trim()) {
                    def parts = line.split(/\s+/)
                    if (parts.size() >= 16) {
                        try {
                            def score = parts[14] as float
                            def evalue = parts[15] as float
                            if (score >= params.cmscan_score && evalue <= params.cmscan_evalue) {
                                hits << [rfam_id: parts[1], score: score, evalue: evalue]
                            }
                        } catch (Exception e) {
                            // Skip malformed lines
                        }
                    }
                }
            }

            if (hits) {
                def best = hits.max { it.score }
                def new_meta = meta + [
                    rfam_id: best.rfam_id,
                    cmscan_score: best.score,
                    cmscan_evalue: best.evalue
                ]
                [new_meta, fasta, 'PASS']
            } else {
                [meta, fasta, 'NO_HIT']
            }
        }
        .branch {
            meta, fasta, status ->
                pass: status == 'PASS'
                    return [meta, fasta]
                fail: true
                    return [meta, fasta]
        }
        .set { ch_cmscan_parsed }

    // Log excluded sequences
    ch_cmscan_parsed.fail
        .subscribe { meta, fasta ->
            log.warn "No Rfam hit found for ${meta.id} - excluding from output"
        }

    //
    // STEP 6: Combine known and inferred sequences for alignment
    //
    ch_branched.known
        .mix(ch_cmscan_parsed.pass)
        .set { ch_for_alignment }

    //
    // STEP 6b: Generate MSA outputs (if --output_msa is enabled)
    //
    if (params.output_msa) {
        // Group sequences by Rfam family
        ch_for_alignment
            .map { meta, fasta -> [meta.rfam_id, fasta] }
            .groupTuple()
            .set { ch_grouped_by_family }

        // Run cmalign on grouped sequences to get proper MSA
        CMALIGN_GROUP(
            ch_grouped_by_family,
            PREPARE_RFAM.out.cm_index.collect()
        )
        ch_versions = ch_versions.mix(CMALIGN_GROUP.out.versions.first().ifEmpty([]))

        // Generate MSA outputs (Stockholm with SS, aligned FASTA, TSV)
        GENERATE_MSA(
            CMALIGN_GROUP.out.stockholm
        )
        ch_versions = ch_versions.mix(GENERATE_MSA.out.versions.first().ifEmpty([]))
    }

    //
    // STEP 7: Run alignment (cmalign or seed-based depending on method)
    //
    if (params.alignment_method == 'seed') {
        // Use pre-computed Rfam seed alignments (like bpRNA)
        SEED_ALIGNMENT(
            ch_for_alignment,
            PREPARE_RFAM.out.seed.collect()
        )
        ch_stockholm = SEED_ALIGNMENT.out.stockholm
        ch_versions = ch_versions.mix(SEED_ALIGNMENT.out.versions.first().ifEmpty([]))
    } else {
        // Default: Use cmalign to compute alignments
        INFERNAL_CMALIGN(
            ch_for_alignment,
            PREPARE_RFAM.out.cm_index.collect()
        )
        ch_stockholm = INFERNAL_CMALIGN.out.stockholm
        ch_versions = ch_versions.mix(INFERNAL_CMALIGN.out.versions.first().ifEmpty([]))
    }

    //
    // STEP 8: Project consensus structure
    //
    CONSENSUS_PROJECTION(
        ch_stockholm
    )
    ch_versions = ch_versions.mix(CONSENSUS_PROJECTION.out.versions.first().ifEmpty([]))

    //
    // STEP 9: Generate summary
    //
    ch_reports = CONSENSUS_PROJECTION.out.report.map { meta, report -> report }

    GENERATE_SUMMARY(
        ch_reports.collect().ifEmpty([]),
        Channel.empty().collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix(GENERATE_SUMMARY.out.versions.ifEmpty([]))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    log.info """
    =========================================
     Pipeline completed!
    =========================================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Output dir: ${params.outdir}
    -----------------------------------------
    """.stripIndent()
}
