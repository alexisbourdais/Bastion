process setup_plasmidfinder {

    label 'process_low'

    publishDir "${params.database}/", mode: 'move'

    output:
    path("plasmidfinder_db/")

    script:
    """
    wget https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/4add282963c7.zip
    mkdir plasmidfinder_db
    unzip *.zip
    rm *.zip
    mv genomicepidemiology-plasmidfinder_db*/* ./plasmidfinder_db/
    cd plasmidfinder_db
    PLASMID_DB=\$(pwd)
    """
}

process plasmidfinder {

    label 'process_low'

    publishDir "${baseDir}/${params.resultsDir}/PlasmidFinder/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.simpleName}.json")

    script:
    """
    results_file="${assembly.simpleName}"
    mkdir \${results_file}
    plasmidfinder.py -i ${assembly} -o \${results_file} -p ${params.db_plasmidfinder}

    mv \${results_file}/data.json ${assembly.simpleName}.json
    """
}