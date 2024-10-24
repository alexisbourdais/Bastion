process setup_Eukcc {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${database}/", mode: 'move'

    output:
    path("eukcc2_db_ver_1.1/")

    script:
    """
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
    tar -xzf eukcc2_db_ver_1.1.tar.gz
    """
}

process eukcc_folder {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Eukcc/eukcc.csv")

    script:
    """
    eukcc folder \
    --threads ${task.cpus} \
    --suffix ".${params.format}" \
    --out "Eukcc" \
    --db ${params.db_eukcc} \
    ${assemblyDir}
    """
}