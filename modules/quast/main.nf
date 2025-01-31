process quast {

    label 'process_medium'

    //publishDir "${params.resultsDir}/Quast/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.baseName}/transposed_report.tsv")

    script:

    if (params.mode_quast=="bacteria" || params.mode_quast =="") {
    """
    quast -o ${assembly.baseName} -t ${task.cpus} ${assembly}
    """
    }

    else {
    """
    quast -o ${assembly.baseName} -t ${task.cpus} --${params.mode_quast} ${assembly}
    """
    }
}