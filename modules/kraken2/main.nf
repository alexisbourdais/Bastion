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

    publishDir "${params.resultsDir}/Kraken/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)

    output:
    path("${assembly.baseName}.report.txt"), emit: krona

    script:
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    kraken2 \
    --use-names ${assembly} \
    --db ${params.db_kraken2} \
    --threads ${task.cpus} \
    --report "\${filename}.report.txt" \
    --output "\${filename}.out"
    """
}

process kraken2_genera {

    label 'contams'
    label 'process_high'

    publishDir "${params.resultsDir}/Kraken/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)
    path(taxdir)

    output:
    path("${assembly.baseName}-parsed.report"), emit: report
    path("${assembly.baseName}.report.txt"), emit: krona

    script:
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    inst-abbr-ids.pl ${assembly} \
    --id-regex=:DEF \
    --id-prefix=\$filename

    inst-split-seqs.pl \$filename-abbr.${params.format} \
    --out=-split

    #labeller
    grep species ${taxdir}/nodes.dmp | cut -f1 > genus.taxid

    create-labeler.pl genus.taxid \
    --taxdir=${taxdir} \
    --level=${params.taxlevel} \
    --kingdoms=Bacteria Archaea Eukaryota \
    > file.idl
    
    kraken2 \
    --use-names \$filename-abbr-split.${params.format} \
    --db ${params.db_kraken2} \
    --threads ${task.cpus} \
    --report "\${filename}.report.txt" \
    --output "\${filename}.out"

    kraken-parser.pl "\${filename}.report.txt" \
    --taxdir=${taxdir} \
    --outfile="\${filename}-parsed.report" \
    --taxon-list=file.idl \
    --auto-detect=${params.automode}
    """
}
