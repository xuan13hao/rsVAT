/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { MERGE_FEATURECOUNTS                } from '../../modules/local/merge_featurecounts'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { checkSamplesAfterGrouping      } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { methodsDescriptionText         } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SUBREAD_FEATURECOUNTS      } from '../../modules/nf-core/subread/featurecounts'
include { MULTIQC                    } from '../../modules/nf-core/multiqc'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_ALIGN_VAT                  } from '../../subworkflows/local/fastq_align_vat'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    ch_samplesheet       // channel: path(sample_sheet.csv)
    ch_versions          // channel: [ path(versions.yml) ]
    ch_fasta             // channel: path(genome.fasta)
    ch_gtf               // channel: path(genome.gtf)
    ch_fai               // channel: path(genome.fai)
    ch_chrom_sizes       // channel: path(genome.sizes)
    ch_gene_bed          // channel: path(gene.bed)
    ch_transcript_fasta  // channel: path(transcript.fasta)
    ch_star_index        // channel: path(star/index/)
    ch_rsem_index        // channel: path(rsem/index/)
    ch_hisat2_index      // channel: path(hisat2/index/)
    ch_vat_index         // channel: path(vat/index/)
    ch_salmon_index      // channel: path(salmon/index/)
    ch_kallisto_index    // channel: [ meta, path(kallisto/index/) ]
    ch_bbsplit_index     // channel: path(bbsplit/index/)
    ch_ribo_db           // channel: path(sortmerna_fasta_list)
    ch_sortmerna_index   // channel: path(sortmerna/index/)
    ch_splicesites       // channel: path(genome.splicesites.txt)

    main:

    ch_multiqc_files = Channel.empty()
    ch_trim_status   = Channel.empty()
    ch_map_status    = Channel.empty()    // unused in the VAT route, kept for the emit signature
    ch_strand_status = Channel.empty()    // unused in the VAT route, kept for the emit signature

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2, genome_bam, transcriptome_bam ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { id, metas, reads ->
            checkSamplesAfterGrouping([ id, metas, reads ])
        }
        .map { meta, reads -> [ meta.findAll { it.key != 'percent_mapped' }, reads ] }
        .set { ch_fastq }

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        ch_sortmerna_index,
        ch_bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit || ! params.fasta,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming,
        params.skip_umi_extract,
        !salmon_index_available,
        false,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        params.remove_ribo_rna,
        params.with_umi,
        params.umi_discard_read,
        params.stranded_threshold,
        params.unstranded_threshold,
        params.skip_linting,
        false
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }

    //
    // SUBWORKFLOW: Alignment with VAT
    //
    FASTQ_ALIGN_VAT (
        ch_strand_inferred_filtered_fastq,
        ch_vat_index.map { [ [:], it ] },
        ch_fasta.map { [ [:], it ] }
    )
    ch_genome_bam       = FASTQ_ALIGN_VAT.out.bam
    ch_genome_bam_index = params.bam_csi_index ? FASTQ_ALIGN_VAT.out.csi : FASTQ_ALIGN_VAT.out.bai
    ch_multiqc_files    = ch_multiqc_files
        .mix(FASTQ_ALIGN_VAT.out.summary.collect{it[1]})
        .mix(FASTQ_ALIGN_VAT.out.stats.collect{it[1]})
        .mix(FASTQ_ALIGN_VAT.out.flagstat.collect{it[1]})
        .mix(FASTQ_ALIGN_VAT.out.idxstats.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_ALIGN_VAT.out.versions)

    //
    // The VAT aligner concatenates R1/R2 into a single read stream and emits a
    // single-end-style BAM without paired-end SAM FLAGs. When vat_single_end_bam
    // is enabled, force downstream featureCounts to treat the BAM as single-end
    // (drop -p/--countReadPairs).
    //
    def vat_force_single_end = params.aligner == 'vat' && params.vat_single_end_bam
    def asSingleEndBam = { meta, bam -> [ vat_force_single_end ? meta + [ single_end: true ] : meta, bam ] }

    //
    // MODULE: Gene-level quantification with featureCounts
    //
    ch_genome_bam
        .map(asSingleEndBam)
        .combine(ch_gtf)
        .set { ch_featurecounts }

    SUBREAD_FEATURECOUNTS (
        ch_featurecounts
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // MODULE: Merge per-sample featureCounts tables into a single gene-level
    // count matrix (Geneid + one column per sample) for DESeq2/edgeR/limma.
    //
    MERGE_FEATURECOUNTS (
        SUBREAD_FEATURECOUNTS.out.counts.map { meta, counts -> counts }.collect()
    )
    ch_versions = ch_versions.mix(MERGE_FEATURECOUNTS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()

    if (!params.skip_multiqc) {

        // Load MultiQC configuration files
        ch_multiqc_config        = Channel.fromPath("$projectDir/workflows/rnaseq/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()

        // Prepare the workflow summary
        ch_workflow_summary = Channel.value(
            paramsSummaryMultiqc(
                paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            )
        ).collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        ch_methods_description = Channel.value(
            methodsDescriptionText(
                params.multiqc_methods_description
                    ? file(params.multiqc_methods_description)
                    : file("$projectDir/workflows/rnaseq/assets/multiqc/methods_description_template.yml", checkIfExists: true)
            )
        ).collectFile(name: 'methods_description_mqc.yaml')

        // Add summary, versions, and methods to the MultiQC input file list
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
            .mix(ch_methods_description)

        // Provide MultiQC with rename patterns to ensure it uses sample names
        // for single-techrep samples not processed by CAT_FASTQ.
        //
        // We only add mappings when the FASTQ simpleName differs from the sample ID.
        // This prevents duplicate/conflicting mappings when multiple samples share
        // the same FASTQ filename in different directories (see #1657).
        //
        // Note: _raw/_trimmed suffixes are handled via extra_fn_clean_exts in multiqc_config.yml
        ch_name_replacements = ch_fastq
            .map{ meta, reads ->
                def paired = reads[0][1] as boolean
                def suffixes = paired ? ['_1', '_2'] : ['']
                def mappings = []

                def fastq1_simplename = file(reads[0][0]).simpleName
                if (fastq1_simplename != meta.id) {
                    mappings << [fastq1_simplename, "${meta.id}${suffixes[0]}"]
                    if (paired) {
                        mappings << [file(reads[0][1]).simpleName, "${meta.id}${suffixes[1]}"]
                    }
                }

                return mappings.collect { it.join('\t') }
            }
            .flatten()
            .collectFile(name: 'name_replacement.txt', newLine: true)
            .ifEmpty([])

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    trim_status    = ch_trim_status    // channel: [id, boolean]
    map_status     = ch_map_status     // channel: [id, boolean]
    strand_status  = ch_strand_status  // channel: [id, boolean]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
