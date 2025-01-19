process setup_Physeter {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("physeter_db/")

    script:
    """
    wget https://figshare.com/ndownloader/files/32405687 -O Cornet-Baurain.tgz
    tar -xzf Cornet-Baurain.tgz
    mkdir physeter_db/
    mv Cornet-2022-GBIO-Figshare/contam-labels.idl Cornet-2022-GBIO-Figshare/life-tqmd-of73.dmnd Cornet-2022-GBIO-Figshare/life-tqmd-of73.gca Cornet-2022-GBIO-Figshare/taxdump-20211206/* physeter_db/
    
    #setup-taxdir.pl --update-cache --taxdir=physeter_db
    """
}

process physeter {

    label 'process_medium'

    publishDir "${params.resultsDir}/Physeter/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)
    path(taxdir)

    output:
    path("${assembly.baseName}-parsed.report")

    script:
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    #Format
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

    #Run Diamond
    mkdir temp
    diamond blastx \
    -d ${params.db_physeter} \
    -q \$filename-abbr-split.${params.format} \
    -o \$filename.blastx \
    -t temp \
    -k 50 \
    -e 1e-10 \
    -f tab \
    -p ${task.cpus}
    
    rm -rf temp
    
    #Run Physeter
    physeter.pl \$filename.blastx \
    --fasta-dir=./ \
    --outfile=\$filename.report \
    --taxdir=${taxdir} \
    --taxon-list=file.idl \
    --auto-detect \
    --kraken \
    --krona
    
    kraken-parser.pl \$filename-kraken.tsv \
    --taxdir=${taxdir} \
    --outfile=\$filename-parsed.report \
    --taxon-list=file.idl \
    --auto-detect=${params.automode}
    """
}
