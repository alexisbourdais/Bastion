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
    path("${assembly.baseName}_BuscoResults.csv"), emit: report

    script:
    if (params.lineage_busco=="auto-lineage") {
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
            jsonRead.py --input "\${filename}"/short_summary.specific*.json
        else
            jsonRead.py --input "\${filename}"/short_summary.generic*.json
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
        else
            jsonRead.py --input "\${filename}"/short_summary.generic*.json
        fi
        """
    }
}

process busco_plot {

script:
    """
    mv Busco/*/short_summary*.txt Busco/

    generate_plot.py -wd Busco
    """
}

process busco_batch {

    label 'process_high'

    publishDir "${params.resultsDir}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Busco/batch_summary.txt"), emit: report
    path("Busco/busco_figure.png")

    script:
    if (params.lineage_busco=="auto-lineage") {
        """
        busco --download_path ${params.db_busco} \
        -i ${assemblyDir} \
        -m ${params.mode_busco} \
        --${params.lineage_busco} \
        -c ${task.cpus} \
        -o Busco

        mv Busco/*/short_summary*.txt Busco/

        generate_plot.py -wd Busco
        """
    }
    else {
        """
        busco --download_path ${params.db_busco} \
        -i ${assemblyDir} \
        -m ${params.mode_busco} \
        -l ${params.lineage_busco} \
        -c ${task.cpus} \
        -o Busco

        mv Busco/*/short_summary*.txt Busco/

        generate_plot.py -wd Busco
        """
    }
}
