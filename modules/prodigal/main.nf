process prodigal {

    label 'process_medium'

    publishDir "${params.resultsDir}/Annotation/", mode: 'copy', pattern: "*.gff3"

    input:
    path(assembly)

    output:
    path(assembly)
    path("${assembly.simpleName}.gff3")

    script:
    """
    prodigal -i ${assembly} -o ${assembly.simpleName}.gff3 -f gff
    """
}