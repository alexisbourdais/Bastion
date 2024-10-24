process krona {

    label 'process_low'

    publishDir "${results}/Krona/", mode: 'copy'

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