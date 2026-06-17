process VAT_ALIGN {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? params.vat_container : '' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val(meta), path("*.vat.log"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_inputs = [reads].flatten().join(' ')
    """
    REF=\$(find -L ${index} -maxdepth 1 -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) | head -n 1)
    if [ -z "\$REF" ]; then
        echo "No FASTA reference found in VAT index directory: ${index}" >&2
        exit 1
    fi

    for read in ${read_inputs}; do
        if [[ "\$read" == *.gz ]]; then
            gzip -cd "\$read"
        else
            cat "\$read"
        fi
    done > ${prefix}.vat_input.fastq

    VAT nucl short RNAseq \\
        -d "\$REF" \\
        -q ${prefix}.vat_input.fastq \\
        -o ${prefix}.sam \\
        -f sam \\
        -p ${task.cpus < 16 ? 16 : task.cpus} \\
        ${args} \\
        > ${prefix}.vat.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vat: \$(VAT --version 2>&1 | grep -oE '[0-9]+[0-9a-z._-]*' | head -n 1 || VAT -version 2>&1 | grep -oE '[0-9]+[0-9a-z._-]*' | head -n 1 || echo "Unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sam
    touch ${prefix}.vat.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vat: stub
    END_VERSIONS
    """
}
