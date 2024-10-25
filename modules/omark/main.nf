process setup_Omark {

    publishDir "${params.database}/", mode: 'move'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    output:
    path("omark_db_LUCA.h5")

    script:
    """
    wget https://omabrowser.org/All/LUCA.h5
    mv LUCA.h5 omark_db_LUCA.h5
    """
}

process omark {

    label 'process_medium'

    publishDir "${params.resultsDir}/Omark/", mode: 'copy'

    input:
    path(proteins)

    output:
    path("${proteins.simpleName}_detailed_summary.txt"), emit: report
    path("${proteins.simpleName}.sum"), emit: plot
    path("${proteins.simpleName}.png")

    script:
    """
    filename=\$(basename -- "${proteins}")
    filename="\${filename%%.*}"

    omamer search --db ${params.db_omark} --query ${proteins} --out \${filename}.tsv --nthreads ${task.cpus}
    omark -f \${filename}.tsv -d ${params.db_omark} -o \${filename}
    """
}