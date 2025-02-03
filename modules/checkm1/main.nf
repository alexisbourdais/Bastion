process setup_Checkm1 {

    label 'process_medium'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries 2

    publishDir "${params.database}/", mode: 'move'

    output:
    path("checkm1_db/")

    script:
    """
    mkdir checkm1_db
    cd checkm1_db
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
    tar -xzf checkm_data_2015_01_16.tar.gz
    rm checkm_data_2015_01_16.tar.gz

    conda env config vars set CHECKM_DATA_PATH="${params.db_checkm1}"
    """
}

process checkm1 {

    label 'process_high'

    //publishDir "${baseDir}/${params.resultsDir}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Checkm1_report.txt")

    script:
    """
    export CHECKM_DATA_PATH="${params.db_checkm1}"

    newDir="GenDirCheckm1"
    mkdir "\${newDir}"

    for file in ${assemblyDir}/*.${params.format}
    do
        filename="\$(basename -- "\${file}")"
        filename="\${filename%.*}"
        cp "\${file}" "\${newDir}/Genome_\${filename}.${params.format}"
    done

    checkm lineage_wf -t ${task.cpus} -x ${params.format} \${newDir} checkm1_outfolder > checkm1.results

    rm -rf "\${newDir}"

    echo "#genome,completeness,contamination,str-hetero" > entete.txt
    tr -s " " < checkm1.results | grep "Genome" | cut -f2,14,15,16 -d" " > resultat.txt
    sed -i -e 's/ /,/g' resultat.txt
    cat entete.txt resultat.txt > report_temp.txt
    
    grep "Genome" report_temp.txt | sed "s/Genome_//g" > Checkm1_report.txt
    """
}

/* A TESTER
//Remplacer le checkm1.results par storage/bin_stats_ext.tsv   
grep "NomAsm" bin_stats_ext.tsv | cut -f1 -d',' => 100_FlyeAsm	{'marker lineage': 'g__Pseudomonas'
grep "Flye" bin_stats_ext.tsv | cut -f11 -d','  => 'Completeness': 93.84585410847941
grep "Flye" bin_stats_ext.tsv | cut -f12 -d','  => 'Contamination': 1.6551383399209487
Manque Strain heterogeneity mais ne sert pas pour le rapport ?!
*/