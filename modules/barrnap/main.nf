process barrnap {

    label 'process_medium'

    publishDir "${baseDir}/${params.resultsDir}/Barrnap/", mode: 'copy'

    input:
    path(assembly)
    val(mode)

    output:
    path("${assembly.simpleName}_rrna.fasta")
    path("${assembly.simpleName}.gff")

    script:
    """
    barrnap --kingdom ${mode} --threads ${task.cpus} --outseq ${assembly.simpleName}_rrna.fasta < ${assembly} > "${assembly.simpleName}.gff"
    """
}