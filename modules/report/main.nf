process final_report {

    label 'contams'

    publishDir "${params.resultsDir}/", mode: 'move'

    input:
    path(busco_report)
    path(quast_report_multi)
    path(checkm2_report)
    path(gunc_report)
    path(checkm1_report)
    path(eukcc_report)
    path(gtdbtk_report)
    //path(omark_report)
    //path(kmerfinder_report)
    //path(kraken2_multi_report)
    //path(physeter_multi_report)

    output:
    path("Bastion_FinalReport.tsv")

    script:
    """
    bastion_companion.py \
    --busco ${busco_report} \
    --quast ${quast_report_multi} \
    --checkm2 ${checkm2_report} \
    --gunc ${gunc_report} \
    --checkm1 ${checkm1_report} \
    --eukcc ${eukcc_report} \
    --gtdbtk ${gtdbtk_report}
    #--omark ${omark_report} \
    #--kmerfinder ${kmerfinder_report} \
    #--physeter ${params} \
    #--kraken ${params}
    """
}