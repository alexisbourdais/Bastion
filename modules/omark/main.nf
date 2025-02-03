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

    publishDir "${baseDir}/${params.resultsDir}/Omark/", mode: 'copy', pattern: "*summary*"
    publishDir "${baseDir}/${params.resultsDir}/Omark/", mode: 'copy', pattern: "*.png"

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
    taxon=\$(sed -n '14p' "${proteins.simpleName}.sum" | cut -f1)
    score=\$(sed -n '14p' "${proteins.simpleName}.sum" | cut -f4)

    if [ -z "\$(sed -n '15p' ${proteins.simpleName}.sum)" ]; then
        conta="no_contam"
    else
        conta="potential_contam"
    fi

    echo -e "\${genome}\t\${taxon}\t\${score}\t\${conta}" > ${proteins.simpleName}_omark_besthit.txt
    """
}

process omark_plot {

    label 'process_low'

    publishDir "${baseDir}/${params.resultsDir}/Omark/", mode: 'copy'

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