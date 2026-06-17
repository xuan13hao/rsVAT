#!/bin/bash
nextflow run /home/h392x566/rsVAT \
  --input /home/h392x566/JZ/sample1.csv \
  --outdir /home/h392x566/JZ/results_mouse_rnaseq4 \
  --fasta /home/h392x566/JZ/ref/mm39.fa \
  --gtf   /home/h392x566/JZ/ref/mm39.ncbiRefSeq.gtf \
  --aligner vat \
  --max_cpus 32 --max_memory '120.GB'\
  --multiqc_title "Mouse_TotalRNA_PE1011_VAT" \
  --vat_single_end_bam true \
  -profile singularity \
  -resume

# The run you asked for:
# bash run_bash2.sh          # uses --aligner vat; vat_single_end_bam defaults to true

# To disable the behavior (e.g. compare/confirm the old failure):
# nextflow run /home/h392x566/rsVAT --input sample1.csv --aligner vat \
#   --vat_single_end_bam false -profile singularity -resume
