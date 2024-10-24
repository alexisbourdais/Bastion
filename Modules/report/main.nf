process final_report {

    label 'contams'

    publishDir "${results}/", mode: 'move'

    input:
    path(busco_report)
    path(quast_report_multi)
    path(checkm2_report)
    path(gunc_report)
    path(checkm1_report)
    path(kraken2_multi_report)
    path(physeter_multi_report)
    path(eukcc_report)
    path(gtdbtk_report)

    output:
    path("GENERA-contams.table")
    path("positive-list.txt")

    script:
    """
    contams_companion.py \
    ${checkm1_report} \
    --mode=table \
    --ckcompleteness=${params.ckcompleteness} \
    --ckcontamination=${params.ckcontamination} \
    --gunccss=${params.gunccss} \
    --guncrrs=${params.guncrrs} \
    --physetercontamination=${params.physetercontamination} \
    --krakencontamination=${params.krakencontamination} \
    --bucompleteness=${params.bucompleteness} \
    --budups=${params.budups} \
    --numcontigs=${params.numcontigs} \
    --ck2completeness=${params.ck2completeness} \
    --ck2contamination=${params.ck2contamination}
    """
}