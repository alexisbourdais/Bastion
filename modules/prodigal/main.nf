process prodigal {

    label 'process_medium'

    publishDir "${params.resultsDir}/Annotation/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.simpleName}.gff3")

    script:
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%%.*}"

    prodigal -i ${assembly} -o \${filename}.gff3 -f gff
    """
}