process final_report {

    publishDir "${params.resultsDir}", mode: 'move'

    input:
    path(busco_report_multi)
    path(quast_report_multi)
    path(checkm2_report)
    path(gunc_report)
    path(checkm1_report)
    path(eukcc_report)
    path(gtdbtk_report)
    path(physeter_multi_report)
    path(kraken2_multi_report)
    path(kmerfinder_report)
    path(omark_report_multi)

    output:
    path("Bastion_FinalReport.tsv")
    path(quast_report_multi)

    script:
    """
    report.py \
    --busco ${busco_report_multi} \
    --quast ${quast_report_multi} \
    --checkm2 ${checkm2_report} \
    --gunc ${gunc_report} \
    --checkm1 ${checkm1_report} \
    --eukcc ${eukcc_report} \
    --gtdbtk ${gtdbtk_report} \
    --physeter ${physeter_multi_report} \
    --kraken ${kraken2_multi_report} \
    --kmerfinder ${kmerfinder_report} \
    --omark ${omark_report_multi}
    """
}
