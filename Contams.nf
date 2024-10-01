#!/usr/bin/env nextflow

/*
===============================================================
 This script is largely based on the Genome quality assessment workflow from the Genera toolbox of Luc CORNET from the University of Liege:
 https://github.com/Lcornet/GENERA
 
 ---------------------------------------------------------------
 Contams Analysis Pipeline. Started September 2024.
 #### Homepage / Documentation
 https://github.com/alexisbourdais/
 #### Authors
 Alexis Bourdais
---------------------------------------------------------------
*/

///////////////////////////////////////////////////////////
//////////////////////       HELP       ///////////////////
///////////////////////////////////////////////////////////

def helpMessage() {
    log.info """

    ------------------------------------------------------------------------------------
    >>>>>   Edit conda, singularity and database path in config file before using  <<<<<
    ------------------------------------------------------------------------------------

    Command : nextflow run Contams.nf -profile slurm,conda [option]

    REQUIRED parameter

    -profile [standard]/slurm,      Select profile standard (local) or slurm. Default: standard          
             singularity/conda      Select profile singularity or conda. 
                                    Physeter/Kraken/Report need singularity in both case.
    
    Database directory
    --db_busco          Path to database directory
    --db_checkm2        Path to database directory
    --db_eukcc          Path to database directory
    --db_gunc           Path to database directory
    --db_kraken2        Path to database directory
    --db_gtdbtk         Path to database directory
    --db_checkm1        Path to database directory
    --db_physeter       Path to database directory

    OPTIONAL parameter
    
    Assembly directory
    --assemblyDir            Default: "./Data/"
    --format                 Default: "fasta"

    Results directory
    --resultsDir            Path to results directory, default: "./Results/"

    Quast
    --mode_quast            [bacteria], eukaryote, fungus

    Busco
    --mode_busco            [genome], proteins
    --lineage_busco         [bacteria_odb10], fungi_odb10, auto-lineage ...

    Eukcc2
    --mode_eukcc            [DNA] ...

    Physeter/Kraken
    --taxlevel              [phylum]
    --automode              [label_first]

    Final report            
    --ckcompleteness        Minimum CheckM completeness, default = 95
    --ckcontamination       Maximum CheckM contamination, default = 5
    --gunccss               Maximum GUNC css, default = 0.01
    --guncrrs               Minimum GUNC rrs, default = 0.5
    --physetercontamination Maximum Physeter contamination, default = 100 (unactivated by default)
    --krakencontamination   Maximum Kraken contamination, default = 100 (unactivated by default)
    --bucompleteness        Minimum BUSCO completeness, default = 0 (unactivated by default)
    --budups                Maximum BUSCO duplication, default = 100 (unactivated by default)
    --ck2completeness       Minimum CheckM2 completeness, default = 95
    --ck2contamination      Maximum CheckM2 contamination, default = 5
    --numcontigs            Maximum Number of contigs, default = 1000


    nextflow run Contams.nf --help

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
        helpMessage()
        exit 0
}

///////////////////////////////////////////////////////////
//////////////////////     Process      ///////////////////
///////////////////////////////////////////////////////////

results = file(params.resultsDir)

process busco {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Busco/batch_summary.txt"), emit: report
    path("Busco/busco_figure.png")

    script:
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

process quast {

    label 'process_medium'

    publishDir "${results}/Quast/", mode: 'copy'

    input:
    path(assembly)

    output:
    path("${assembly.baseName}/transposed_report.tsv")

    script:

    if (params.mode_quast=="bacteria" || params.mode_quast =="") {
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    quast -o \${filename} -t ${task.cpus} ${assembly}
    """
    }

    else {
    """
    filename=\$(basename -- "${assembly}")
    filename="\${filename%.*}"

    quast -o \${filename} -t ${task.cpus} --${params.mode_quast} ${assembly}
    """
    }
}

process checkm2 {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Checkm2/quality_report.tsv")

    script:
    """
    checkm2 predict \
    --threads ${task.cpus} \
    -x ${params.format} \
    --input ${assemblyDir} \
    --output-directory "Checkm2" \
    --database_path ${params.db_checkm2} \
    --force
    """
}

process eukcc_folder {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Eukcc/eukcc.csv")

    script:
    """
    eukcc folder \
    --threads ${task.cpus} \
    --suffix ".${params.format}" \
    --out "Eukcc" \
    --db ${params.db_eukcc} \
    ${assemblyDir}
    """
}

process gunc {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Gunc/GUNC.progenomes_2.1.maxCSS_level.tsv"), emit: report
    path("Gunc/diamond_output/*.html")

    script:
    """
    dir_res="Gunc/"
    mkdir "\${dir_res}"

    gunc run \
    --db_file ${params.db_gunc} \
    --input_dir ${assemblyDir} \
    --file_suffix ${params.format} \
    --threads ${task.cpus} \
    --out_dir "\${dir_res}"


    diamond_results="\${dir_res}diamond_output/"
    cd \${diamond_results}
    for file in *.out
    do
        gunc plot --diamond_file \${file}
    done

    """
}

process gtdbtk {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Gtdbtk/gtdbtk.bac120.summary.tsv"), emit: report
    path("Gtdbtk/identify")

    script:
    """
    conda env config vars set GTDBTK_DATA_PATH="${params.db_gtdbtk}"

    gtdbtk classify_wf \
    --genome_dir ${assemblyDir} \
    --out_dir "Gtdbtk" \
    --cpus ${task.cpus} \
    --mash_db "Gtdbtk" \
    --extension ${params.format}
    """
}

process checkm1 {

    label 'process_high'

    publishDir "${results}/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assemblyDir)

    output:
    path("Checkm1_report.txt")

    script:
    """
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

process kraken2 {

    label 'contams'
    label 'process_high'

    //publishDir "${results}/Kraken/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)

    output:
    path("${assembly.baseName}-parsed.report"), emit: report
    path("${assembly.baseName}.report.txt"), emit: krona

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
    grep species ${params.taxdir}/nodes.dmp | cut -f1 > genus.taxid

    create-labeler.pl genus.taxid \
    --taxdir=${params.taxdir} \
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
    --taxdir=${params.taxdir} \
    --outfile="\${filename}-parsed.report" \
    --taxon-list=file.idl \
    --auto-detect=${params.automode}
    """
}

process krona {

    label 'process_low'

    publishDir "${results}/Krona/", mode: 'copy'

    input:
    path(report)

    output:
    path("${report.simpleName}_krona.html")

    script:
    """
    filename=\$(basename -- "${report}")
    filename="\${filename%%.*}"

    ktImportTaxonomy -t 5 -m 3 -o \${filename}_krona.html ${report}
    """
}

process physeter {

    label 'contams'
    label 'process_medium'

    publishDir "${results}/Physeter/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }

    input:
    path(assembly)

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
    grep species ${params.taxdir}/nodes.dmp | cut -f1 > genus.taxid

    create-labeler.pl genus.taxid \
    --taxdir=${params.taxdir} \
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
    --taxdir=${params.taxdir} \
    --taxon-list=file.idl \
    --auto-detect \
    --kraken
    
    kraken-parser.pl \$filename-kraken.tsv \
    --taxdir=${params.taxdir} \
    --outfile=\$filename-parsed.report \
    --taxon-list=file.idl \
    --auto-detect=${params.automode}
    """
}

process final_report {

    label 'contams'

    publishDir "${results}/", mode: 'move'

    input:
    path(busco_report)
    path(quast_report_multi)
    path(checkm2_report)
    path(gunc_report)
    path(checkm1_report)
    path(kraken2_multi_report)
    path(physeter_multi_report)
    path(eukcc_report)
    path(gtdbtk_report)

    output:
    path("GENERA-contams.table")
    path("positive-list.txt")

    script:
    """
    contams_companion.py \
    ${checkm1_report} \
    --mode=table \
    --ckcompleteness=${params.ckcompleteness} \
    --ckcontamination=${params.ckcontamination} \
    --gunccss=${params.gunccss} \
    --guncrrs=${params.guncrrs} \
    --physetercontamination=${params.physetercontamination} \
    --krakencontamination=${params.krakencontamination} \
    --bucompleteness=${params.bucompleteness} \
    --budups=${params.budups} \
    --numcontigs=${params.numcontigs} \
    --ck2completeness=${params.ck2completeness} \
    --ck2contamination=${params.ck2contamination}
    """
}

///////////////////////////////////////////////////////
//////////////////     Workflow     ///////////////////
///////////////////////////////////////////////////////

workflow {

    data_file = Channel.fromPath("${params.assemblyFiles}")
    data_dir = Channel.fromPath("${params.assemblyDir}")
    gunc(data_dir)
    busco(data_dir)
    quast(data_file)
    eukcc_folder(data_dir)
    checkm2(data_dir)
    gtdbtk(data_dir)
    checkm1(data_dir)
    kraken2(data_file)
    krona(kraken2.out.krona)
    physeter(data_file)

    final_report(busco.out.report, quast.out.collectFile(name: 'quast_multi.report'), \
    checkm2.out, gunc.out.report, checkm1.out, kraken2.out.report.collectFile(name: 'kraken2_multi.report'), \
    physeter.out.collectFile(name: 'physeter_multi.report'), eukcc_folder.out, gtdbtk.out.report)
}

workflow.onComplete{println("Workflow execution finished")}