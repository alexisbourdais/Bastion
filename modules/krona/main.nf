process setup_Krona {

    label 'process_low'

    publishDir "${params.database}/", mode: 'move'

    output:
    path("Krona_db/")

    script:
    """
    mkdir Krona_db
    cd Krona_db
    ktUpdateTaxonomy.sh ./
    """
}

process krona {

    label 'process_low'

    publishDir "${baseDir}/${params.resultsDir}/Krona/", mode: 'copy'

    input:
    path(report)
    val(origin)

    output:
    path("${report.simpleName}_${origin}_krona.html")

    script:
    """
    ktImportTaxonomy -tax "${params.db_krona}" -t 5 -m 3 -o ${report.simpleName}_${origin}_krona.html ${report}
    """
}