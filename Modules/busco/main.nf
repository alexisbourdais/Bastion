process setup_Busco {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${database}/", mode: 'move'

    output:
    path("busco_db/")

    script:
    """
    busco --download prokaryota virus
    mv busco_downloads/ busco_db/
    """
}

process busco {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Busco/batch_summary.txt"), emit: report
    path("Busco/busco_figure.png")

    script:
    """
    busco --download_path ${params.db_busco} \
    -i ${assemblyDir} \
    -m ${params.mode_busco} \
    -l ${params.lineage_busco} \
    -c ${task.cpus} \
    -o Busco

    mv Busco/*/short_summary*.txt Busco/

    generate_plot.py -wd Busco
    """
}