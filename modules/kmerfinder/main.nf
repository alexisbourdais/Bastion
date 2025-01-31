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

    publishDir "${params.resultsDir}/Kmerfinder/", mode: 'copy', pattern: "*_kmerfinder_results.txt"
    
    input:
    path(assembly)

    output:
    path("${assembly.simpleName}_kmerfinder_results.txt")
    path("${assembly.simpleName}_kmerfinder_tophit.txt"), emit: report

    script:
    """
    kmerfinder.py -i ${assembly} -o ${assembly.simpleName} -db "${params.db_kmerfinder}${params.taxon_kmerfinder}/${params.taxon_kmerfinder}.ATG" -tax "${params.db_kmerfinder}${params.taxon_kmerfinder}/${params.taxon_kmerfinder}.tax" -x
    mv ${assembly.simpleName}/results.txt ./${assembly.simpleName}_kmerfinder_results.txt

    sed -n '1,2p' "${assembly.simpleName}_kmerfinder_results.txt" | awk -v prefix="${assembly.simpleName}" '{print prefix "\t" \$0}' > ${assembly.simpleName}_kmerfinder_tophit.txt
    """
}