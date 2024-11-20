process setup_Krona {

    label 'process_low'

    script:
    """
    ktUpdateTaxonomy.sh
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

    ktImportTaxonomy -t 5 -m 3 -o \${filename}_krona.html ${report}
    """
}