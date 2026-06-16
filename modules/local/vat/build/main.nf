process VAT_BUILD {
    tag "$fasta"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? params.vat_container : '' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("vat"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir vat
    cp ${fasta} vat/${fasta.name}

    cd vat
    VAT makevatdb \\
        --in ${fasta.name} \\
        --dbtype nucl \\
        -p ${task.cpus} \\
        ${args}
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vat: \$(VAT --version 2>&1 | head -n 1 | sed 's/^VAT[ _-]*//; s/^version[ _-]*//' || true)
    END_VERSIONS
    """

    stub:
    """
    mkdir vat
    touch vat/${fasta.name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vat: stub
    END_VERSIONS
    """
}
