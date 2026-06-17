process MERGE_FEATURECOUNTS {
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path ('featurecounts/*')

    output:
    path "gene_counts_matrix.tsv", emit: matrix
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p tmp

    # featureCounts output layout (per sample, one BAM):
    #   line 1  -> '# Program:...' comment
    #   line 2  -> header (Geneid Chr Start End Strand Length <bam>)
    #   3..N    -> data, gene id in col 1, count in the last (7th) col
    #
    # All per-sample tables come from the same GTF/featureCounts run, so the
    # gene order is identical across samples and we can paste columns together.

    first=\$(ls ./featurecounts/* | head -n 1)
    echo "Geneid" > geneids.txt
    grep -v '^#' "\$first" | tail -n +2 | cut -f1 >> geneids.txt

    for fileid in \$(ls ./featurecounts/* | sort); do
        samplename=\$(basename "\$fileid" | sed 's/\\.featureCounts\\.tsv\$//g')
        echo "\$samplename" > tmp/\${samplename}.counts.txt
        grep -v '^#' "\$fileid" | tail -n +2 | cut -f7 >> tmp/\${samplename}.counts.txt
    done

    paste geneids.txt tmp/*.counts.txt > gene_counts_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch gene_counts_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(echo \$(sed --version 2>&1) | sed 's/^.*GNU sed) //; s/ .*\$//')
    END_VERSIONS
    """
}
