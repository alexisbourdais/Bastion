process gffread {

    label 'process_medium'

    publishDir "${params.resultsDir}/Annotation/", mode: 'copy'

    input:
    path(assembly)
    path(annotation)

    output:
    path("${assembly.simpleName}_proteins.fasta")

    script:
    """
    gffread -g ${assembly} -y "${assembly.simpleName}_proteins.fasta" ${annotation}
    """
}