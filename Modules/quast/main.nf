process quast {

    label 'process_medium'

    publishDir "${results}/Quast/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.baseName}/transposed_report.tsv")

    script:

    if (params.mode_quast=="bacteria" || params.mode_quast =="") {
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    quast -o \${filename} -t ${task.cpus} ${assembly}
    """
    }

    else {
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    quast -o \${filename} -t ${task.cpus} --${params.mode_quast} ${assembly}
    """
    }
}