process quast {

    label 'process_medium'

    //publishDir ""${baseDir}/${params.resultsDir}/Quast/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.baseName}/transposed_report.tsv")

    script:
    """
    quast -o ${assembly.baseName} -t ${task.cpus} ${assembly}
    """
}