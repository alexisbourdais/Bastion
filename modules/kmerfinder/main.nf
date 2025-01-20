process setup_Kmerfinder {

    label 'process_medium'

    publishDir "${params.database}/", mode: 'move'

    output:
    path("kmerfinder_db/")

    script:
    """
    wget https://cge.food.dtu.dk/services/KmerFinder/etc/kmerfinder_db.tar.gz
    tar -xzvf kmerfinder_db.tar.gz
    """
}

process kmerfinder {

    label 'process_medium'

    publishDir "${params.resultsDir}/Kmerfinder/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.simpleName}/")

    script:
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%%.*}"

    kmerfinder.py -i ${assembly} -o \${filename} -db "${params.db_kmerfinder}${params.taxon_kmerfinder}/${params.taxon_kmerfinder}.ATG" -tax "${params.db_kmerfinder}${params.taxon_kmerfinder}/${params.taxon_kmerfinder}.tax" -x
    """
}