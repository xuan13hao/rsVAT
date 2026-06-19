# rsVAT — VAT RNA-seq quantification pipeline
nf-core/rnaseq is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet with FASTQ files or pre-aligned BAM files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.
This repository is a trimmed-down fork of [nf-core/rnaseq](https://nf-co.re/rnaseq).


## Pipeline overview

```
FASTQ (R1[, R2])
   │  trim / QC
   ▼
VAT_ALIGN 
   ▼
SAMTOOLS_VIEW ─── SAM → BAM
   ▼
BAM_SORT_STATS_SAMTOOLS ─── coordinate-sort + index + stats
   ▼
SUBREAD_FEATURECOUNTS 
   ▼
MERGE_FEATURECOUNTS ─── per-sample tables → gene_counts_matrix.tsv
   ▼
MULTIQC ─── aggregates FastQC/trimming, alignment stats, featureCounts summary
```

---

## Quick start

```bash
nextflow run /path/to/rsVAT \
  --input    samplesheet.csv \
  --outdir   results \
  --fasta    /path/to/genome.fa \
  --gtf      /path/to/genes.gtf \
  --aligner vat\
  -profile   singularity \
```

`--aligner vat` is the default, so it does not need to be passed. VAT runs inside
the container referenced by `--vat_container`, so use `-profile singularity` (or
`apptainer`).

A ready-made launch script is provided in [`run_bash2.sh`](run_bash2.sh).

### Samplesheet

Standard nf-core/rnaseq format. Paired-end and single-end rows may be mixed.

```csv
sample,fastq_1,fastq_2,strandedness
AT54-L,/data/AT54-L_R1.fastq.gz,/data/AT54-L_R2.fastq.gz,auto
AT54-R,/data/AT54-R_R1.fastq.gz,/data/AT54-R_R2.fastq.gz,auto
```

For paired-end rows, R1 and R2 are aligned together as one merged single-end
stream (see below). For single-end rows the single FASTQ is used as-is.

---

## VAT-specific parameters

| Parameter             | Default                       | Description |
|-----------------------|-------------------------------|-------------|
| `--aligner`           | `vat`                         | Alignment route. `vat` is the only accepted value. |
| `--vat_index`         | `null`                        | Path to a pre-built VAT index directory. If omitted, an index is built from `--fasta` with `VAT makevatdb`. |
| `--vat_container`     | `bin/VAT_latest.sif`          | Singularity/Apptainer image that provides the `VAT` binary. |
| `--vat_single_end_bam`| `true`                        | Treat the VAT BAM as single-end in downstream featureCounts, even when the input FASTQ was paired-end. |

---

## Output

The merged gene-level count matrix is published to:

```
<outdir>/vat/featurecounts/gene_counts_matrix.tsv
```

with **genes as rows** and **sample IDs as columns**:

```
Geneid	AT54-L	AT54-R
TrnP	3	1
TrnT	0	0
CYTB	101	114
...
```

The same directory also contains the per-sample `*.featureCounts.tsv` tables and
their `.summary` files. The sorted BAMs, FastQC/trimming reports, and the
aggregated MultiQC report are published under `<outdir>` as usual.

---

## Citations

If you use the VAT alignment route, please also cite:

> Xuan, Hao, et al. "A general and extensible algorithmic framework to biological sequence alignment across scales and applications." _bioRxiv_ (2026): 2026-01.

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat
