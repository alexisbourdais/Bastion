process setup_Kraken2 {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("kraken2_db/")

    script:
    """
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20240605.tar.gz
    tar -xzf k2_pluspfp_16gb_20240605.tar.gz
    rm k2_pluspfp_16gb_20240605.tar.gz
    mkdir kraken2_db
    mv * kraken2_db/
    """
}

process kraken2 {

    label 'process_high'

    publishDir "${baseDir}/${params.resultsDir}/Kraken/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    tuple val(name), path(assembly_split)

    output:
    path("${name}-split.report.txt")

    script:
    """
    kraken2 \
    --use-names ${assembly_split} \
    --db ${params.db_kraken2} \
    --threads ${task.cpus} \
    --report "${name}-split.report.txt" \
    --output "${name}.out"
    """
}

process kraken_split {

    label 'process_low'

    //publishDir "${baseDir}/${params.resultsDir}/Kraken_split/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(report)
    path(idl)
    path(taxdir)

    output:
    path("${report.baseName}-parsed.report")

    script:
    """    
    kraken-parser.pl "${report}" \
    --taxdir=${taxdir} \
    --outfile="${report.baseName}-parsed.report" \
    --taxon-list=${idl} \
    --auto-detect=${params.automode}
    """
}