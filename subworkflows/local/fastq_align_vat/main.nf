include { VAT_ALIGN             } from '../../../modules/local/vat/align'
include { SAMTOOLS_VIEW         } from '../../../modules/local/samtools/view'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_VAT {

    take:
    reads    // channel: [ val(meta), [ reads ] ]
    index    // channel: [ val(meta), path(vat/index) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = Channel.empty()

    VAT_ALIGN(reads, index)
    ch_versions = ch_versions.mix(VAT_ALIGN.out.versions.first())

    SAMTOOLS_VIEW(VAT_ALIGN.out.sam)
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    BAM_SORT_STATS_SAMTOOLS(SAMTOOLS_VIEW.out.bam, ch_fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    sam      = VAT_ALIGN.out.sam                    // channel: [ val(meta), sam ]
    orig_bam = SAMTOOLS_VIEW.out.bam                // channel: [ val(meta), bam ]
    summary  = VAT_ALIGN.out.summary                // channel: [ val(meta), log ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}
