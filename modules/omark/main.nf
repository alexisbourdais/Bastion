process setup_Omark {

    publishDir "${params.database}/", mode: 'move'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    output:
    path("omark_db_LUCA.h5")

    script:
    """
    wget https://omabrowser.org/All/LUCA.h5
    mv LUCA.h5 omark_db_LUCA.h5
    """
}

process omark {

    label 'process_medium'

    publishDir "${params.resultsDir}/Omark/", mode: 'copy', pattern: "${proteins.simpleName}/${proteins.simpleName}_detailed_summary.txt"
    publishDir "${params.resultsDir}/Omark/", mode: 'copy', pattern: "${proteins.simpleName}/${proteins.simpleName}.png"

    input:
    path(proteins)

    output:
    path("${proteins.simpleName}/${proteins.simpleName}_detailed_summary.txt")
    path("${proteins.simpleName}/${proteins.simpleName}.png")
    path("${proteins.simpleName}/"), emit: plot
    path("${proteins.simpleName}_omark_besthit.txt"), emit: report

    script:
    """
    omamer search --db ${params.db_omark} --query ${proteins} --out ${proteins.simpleName}.tsv --nthreads ${task.cpus}
    omark -f ${proteins.simpleName}.tsv -d ${params.db_omark} -o ${proteins.simpleName}

    cp ${proteins.simpleName}/${proteins.simpleName}.sum ./
    genome=\$(echo "${proteins.simpleName}" | sed 's/_proteins//')

    if [ -z "\$(sed -n '15p' ${proteins.simpleName}.sum)" ]; then
        sed -n '14p' "${proteins.simpleName}.sum" | cut -f1,4 | while read field1 field4; do echo -e "\${genome};\${field1};\${field4};no_contam"; done >> "${proteins.simpleName}_omark_besthit.txt"
    else
        sed -n '14p' "${proteins.simpleName}.sum" | cut -f1,4 | while read field1 field4; do echo -e "\${genome};\${field1};\${field4};potential_contam"; done >> "${proteins.simpleName}_omark_besthit.txt"
    fi
    """
}

process omark_plot {

    label 'process_low'

    publishDir "${params.resultsDir}/Omark/", mode: 'copy'

    input:
    path(omark_results_multi)

    output:
    path("omark_multiple.png")

    script:
    """
    mkdir omark_results_dir
    mv ${omark_results_multi} ./omark_results_dir/

    omark_plot_all_results.py --input ./omark_results_dir/
    """
}