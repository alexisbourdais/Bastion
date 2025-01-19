process setup_Busco {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("busco_db/")

    script:
    """
    busco --download prokaryota virus
    mv busco_downloads/ busco_db/
    """
}

process busco {

    label 'process_high'

    publishDir "${params.resultsDir}/Busco", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)

    output:
    path("${assembly.baseName}_BuscoResults.csv"), emit : report
    path("short_summary.*.txt"), emit : plot

    script:
    if (params.lineage_busco=="auto-lineage" || params.lineage_busco=="auto-lineage*") {
        """
        filename=\$(basename -- "${assembly}")
        filename="\${filename%.*}"

        busco --download_path ${params.db_busco} \
        -i ${assembly} \
        -m ${params.mode_busco} \
        --${params.lineage_busco} \
        -c ${task.cpus} \
        -o \${filename}

        if [ -f "\${filename}"/short_summary.specific.*.json ]; then
            jsonRead.py --input "\${filename}"/short_summary.specific.*.json
            mv "\${filename}"/short_summary.specific.*.txt ./
        else
            jsonRead.py --input "\${filename}"/short_summary.generic.*.json
            mv "\${filename}"/short_summary.generic.*.txt ./
        fi
        """
    }
    
    else {
        """
        filename=\$(basename -- "${assembly}")
        filename="\${filename%.*}"

        busco --download_path ${params.db_busco} \
        -i ${assembly} \
        -m ${params.mode_busco} \
        -l ${params.lineage_busco} \
        -c ${task.cpus} \
        -o \${filename}

        if [ -f "\${filename}"/short_summary.specific.*.json ]; then
            jsonRead.py --input "\${filename}"/short_summary.specific*.json
            mv "\${filename}"/short_summary.specific.*.txt ./
        else
            jsonRead.py --input "\${filename}"/short_summary.generic*.json
            mv "\${filename}"/short_summary.generic.*.txt ./
        fi
        """
    }
}

process busco_plot {

    label 'process_low'

    publishDir "${params.resultsDir}/Busco", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(short_summary_multi)

    output:
    path("busco_figure.png"), emit : report

    script:
        """
        mkdir short_summary_dir
        mv ${short_summary_multi} ./short_summary_dir/

        generate_plot.py -wd ./short_summary_dir/
        mv ./short_summary_dir/busco_figure.png ./
        """
}