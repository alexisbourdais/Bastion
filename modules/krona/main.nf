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

    publishDir "${params.resultsDir}/Krona/", mode: 'copy'

    input:
    path(report)

    output:
    path("${report.simpleName}_krona.html")

    script:
    """
    filename=\$(basename -- "${report}")
    filename="\${filename%%.*}"

    ktImportTaxonomy -tax "${params.db_krona}" -t 5 -m 3 -o \${filename}_krona.html ${report}
    """
}