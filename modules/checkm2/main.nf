process setup_Checkm2 {

    publishDir "${params.database}/", mode: 'move'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    output:
    path("checkm2_uniref100.KO.1.dmnd")

    script:
    """
    checkm2 database --download --path ./
    mv CheckM2_database/uniref100.KO.1.dmnd checkm2_uniref100.KO.1.dmnd
    """
}

process checkm2 {

    label 'process_high'

    publishDir "${baseDir}/${params.resultsDir}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Checkm2/quality_report.tsv")

    script:
    """
    checkm2 predict \
    --threads ${task.cpus} \
    -x ${params.format} \
    --input ${assemblyDir} \
    --output-directory "Checkm2" \
    --database_path ${params.db_checkm2} \
    --force
    """
}