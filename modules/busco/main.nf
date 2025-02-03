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

    //publishDir "${baseDir}/${params.resultsDir}/Busco_${mode}", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)
    val(mode)

    output:
    path("${assembly.baseName}_BuscoResults.csv"), optional: true, emit: report
    path("short_summary.*.txt"), emit : plot

    script:
    if (params.lineage_busco=="auto-lineage" || params.lineage_busco=="auto-lineage*") {
        """
        busco --download_path ${params.db_busco} \
        -i ${assembly} \
        -m ${mode} \
        --${params.lineage_busco} \
        -c ${task.cpus} \
        -o ${assembly.baseName}

        if [ -f "${assembly.baseName}"/short_summary.specific.*.json ]; then
            if [ $mode == "genome" ]; then
                jsonRead.py --format ${params.format} --input "${assembly.baseName}"/short_summary.specific.*.json
            fi
            mv "${assembly.baseName}"/short_summary.specific.*.txt ./
        else
            if [ $mode == "genome" ]; then
                jsonRead.py --format ${params.format} --input "${assembly.baseName}"/short_summary.generic.*.json
            fi
            mv "${assembly.baseName}"/short_summary.generic.*.txt ./
        fi
        """
    }
    
    else {
        """
        busco --download_path ${params.db_busco} \
        -i ${assembly} \
        -m ${mode} \
        -l ${params.lineage_busco} \
        -c ${task.cpus} \
        -o ${assembly.baseName}

        if [ -f "${assembly.baseName}"/short_summary.specific.*.json ]; then
            if [ $mode == "genome" ]; then
                jsonRead.py --format ${params.format} --input "${assembly.baseName}"/short_summary.specific*.json
            fi
            mv "${assembly.baseName}"/short_summary.specific.*.txt ./
        else
            if [ $mode == "genome" ]; then
                jsonRead.py --format ${params.format} --input "${assembly.baseName}"/short_summary.generic*.json
            fi
            mv "${assembly.baseName}"/short_summary.generic.*.txt ./
        fi
        """
    }
}

process busco_plot {

    label 'process_low'

    publishDir "${baseDir}/${params.resultsDir}/Busco", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(short_summary_multi)
    val(mode)

    output:
    path("busco_${mode}_figure.png")

    script:
        """
        mkdir short_summary_dir
        mv ${short_summary_multi} ./short_summary_dir/

        generate_plot.py -wd ./short_summary_dir/
        mv ./short_summary_dir/busco_figure.png ./busco_${mode}_figure.png
        """
}