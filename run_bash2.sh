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
  -profile singularity 
