process setup_Gunc {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("gunc_db_progenomes2.1.dmnd")

    script:
    """
    gunc download_db ./
    """
}

process gunc {

    label 'process_high'

    publishDir "${baseDir}/${params.resultsDir}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Gunc/GUNC.progenomes_2.1.maxCSS_level.tsv"), emit: report
    path("Gunc/diamond_output/*.html")

    script:
    """
    dir_res="Gunc/"
    mkdir "\${dir_res}"

    gunc run \
    --db_file ${params.db_gunc} \
    --input_dir ${assemblyDir} \
    --file_suffix ${params.format} \
    --threads ${task.cpus} \
    --out_dir "\${dir_res}"


    diamond_results="\${dir_res}diamond_output/"
    cd \${diamond_results}
    for file in *.out
    do
        gunc plot --diamond_file \${file}
    done
    """
}