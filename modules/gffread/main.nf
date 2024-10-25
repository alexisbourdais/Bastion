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
    filename=\$(basename -- "${assembly}")
    filename="\${filename%%.*}"

    gffread -g ${assembly} -y "\${filename}_proteins.fasta" ${annotation}
    """
}