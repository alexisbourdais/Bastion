process setup_Gtdbtk {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("gtdbtk_db/")

    script:
    """
    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
    tar -xzf gtdbtk_data.tar.gz
    mv release*/ gtdbtk_db/

    #download-db.sh -> not latest db

    conda env config vars set GTDBTK_DATA_PATH="${params.db_gtdbtk}"
    """
}

process gtdbtk {

    label 'process_high'

    publishDir "${params.resultsDir}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Gtdbtk/gtdbtk.bac120.summary.tsv"), emit: report
    path("Gtdbtk/identify")

    script:
    """
    export GTDBTK_DATA_PATH="${params.db_gtdbtk}"

    gtdbtk classify_wf \
    --genome_dir ${assemblyDir} \
    --out_dir "Gtdbtk" \
    --cpus ${task.cpus} \
    --extension ${params.format}
    """
}