process setup_Physeter {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("physeter_db/")
    path("taxdump/")

    script:
    """
    wget https://figshare.com/ndownloader/files/32405687 -O Cornet-Baurain.tgz
    tar -xzf Cornet-Baurain.tgz
    mkdir physeter_db/
    mv Cornet-2022-GBIO-Figshare/life-tqmd-of73* ./physeter_db/
    mv Cornet-2022-GBIO-Figshare/contam-labels.idl ./physeter_db/

    mkdir taxdump
    setup-taxdir.pl --taxdir=./taxdump/
    """
}

process physeter {

    label 'process_medium'

    publishDir "${params.resultsDir}/Physeter/", mode: 'copy', pattern: "${assembly.baseName}.tsv"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)
    path(taxdir)

    output:
    path("${assembly.baseName}-parsed.report"), emit: report
    path("${assembly.baseName}.tsv"), emit: krona
    tuple val("${assembly.baseName}"), path("${assembly.baseName}-abbr-split.fasta"), emit: kraken
    path("file.idl"), emit: kraken_split

    script:
    """
    #Format
    inst-abbr-ids.pl ${assembly} \
    --id-regex=:DEF \
    --id-prefix=${assembly.baseName}

    inst-split-seqs.pl ${assembly.baseName}-abbr.${params.format} \
    --out=-split

    #labeller
    grep species ${taxdir}/nodes.dmp | cut -f1 > genus.taxid

    create-labeler.pl genus.taxid \
    --taxdir=${taxdir} \
    --level=${params.taxlevel} \
    --kingdoms=Bacteria Archaea Eukaryota \
    > file.idl

    #Run Diamond
    mkdir temp
    diamond blastx \
    -d ${params.db_physeter} \
    -q ${assembly.baseName}-abbr-split.${params.format} \
    -o ${assembly.baseName}-abbr-split.blastx \
    -t temp \
    -k 50 \
    -e 1e-10 \
    -f tab \
    -p ${task.cpus}
    
    rm -rf temp
    
    #Run Physeter
    if [ "${params.format}" != "fasta" ]; then
        mv ${assembly.baseName}-abbr-split.${params.format} ${assembly.baseName}-abbr-split.fasta
    fi

    physeter.pl ${assembly.baseName}-abbr-split.blastx \
    --fasta-dir=./ \
    --outfile=${assembly.baseName}.report \
    --taxdir=${taxdir} \
    --taxon-list=file.idl \
    --auto-detect \
    --kraken

    kraken-parser.pl ${assembly.baseName}-abbr-split-kraken.tsv \
    --taxdir=${taxdir} \
    --outfile=${assembly.baseName}-parsed.report \
    --taxon-list=file.idl \
    --auto-detect=${params.automode}

    mv ${assembly.baseName}-abbr-split-kraken.tsv ${assembly.baseName}.tsv
    """
}
