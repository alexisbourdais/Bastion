#!/usr/bin/env nextflow

/*
===============================================================
 Bastion Pipeline. Started September 2024.
 #### Homepage / Documentation
 https://github.com/alexisbourdais/Bastion/
 #### Authors
 Alexis Bourdais
---------------------------------------------------------------
*/

///////////////////////////////////////////////////////////
//////////////////////       HELP       ///////////////////
///////////////////////////////////////////////////////////

def helpMessage() {
    log.info """

    Command : nextflow run Bastion.nf -profile slurm,singularity [option] --workflow

    REQUIRED parameter

    -profile [standard]/slurm,      Select profile standard (local) or slurm. Default: standard          
             singularity/conda      Select profile singularity or conda. 
                                    Physeter need singularity in both case and busco's conda env has mistakes.
    
    --workflow                      Select workflow :  'setup' to download database
                                                       'analysis' to run all analyses

    OPTIONAL parameter

    -resume                         used to resume a workflow from where it was previously stopped or interrupted

    Choose which database you want to set-up when using --workflow setup
    --setAll
    --setEukcc
    --setBusco
    --setGtdbtk
    --setKraken2
    --setKrona
    --setCheckm1
    --setGunc
    --setCheckm2
    --setPhyseter
    --setOmark
    --setKmerfinder
    --setPlasmidfinder

    Path to database directory : automatic if installed with --workflow setup
    --db_busco          Path to database directory
    --db_checkm2        Path to checkm2_uniref100.KO.1.dmnd
    --db_eukcc          Path to database directory
    --db_gunc           Path to gunc_db_progenomes2.1.dmnd
    --db_kraken2        Path to database directory
    --db_krona          Path to database directory
    --db_gtdbtk         Path to database directory
    --db_checkm1        Path to database directory
    --db_omark          Path to database directory
    --db_physeter       Path to life-tqmd-of73.dmnd
    --taxdump           Path to taxdump
    --db_kmerfinder     Path to database directory
    --db_plasmidfinder  Path to database directory

    Assembly directory
    --assemblyDir            Default: "./Data/"
    --format                 Default: "fasta"

    Results directory
    --resultsDir            Path to results directory, default: "Results/"

    Busco
    --lineage_busco         [auto-lineage], auto-lineage-prok, bacteria_odb10

    Physeter
    --taxlevel              [phylum]
    --automode              [label_first]

    nextflow run Bastion.nf --help
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
        helpMessage()
        exit 0
}

//////////////////////////////////////////////
////////////    Modules     //////////////////
//////////////////////////////////////////////

include { setup_Eukcc; eukcc_folder }                                       from './modules/eukcc'
include { setup_Checkm1; checkm1 }                                          from './modules/checkm1'
include { setup_Gunc; gunc }                                                from './modules/gunc'
include { setup_Kraken2; kraken2; kraken_split }                            from './modules/kraken2'
include { setup_Gtdbtk; gtdbtk }                                            from './modules/gtdbtk'
include { setup_Physeter; physeter }                                        from './modules/physeter'
include { setup_Busco; busco; busco_plot }                                  from './modules/busco'
include { setup_Checkm2; checkm2 }                                          from './modules/checkm2'
include { setup_Omark; omark; omark_plot }                                  from './modules/omark'
include { quast }                                                           from './modules/quast'
include { setup_Krona; krona as krona_kraken; krona as krona_physeter }     from './modules/krona'
include { prodigal }                                                        from './modules/prodigal'
include { gffread }                                                         from './modules/gffread'
include { final_report }                                                    from './modules/report'
include { setup_Kmerfinder; kmerfinder }                                    from './modules/kmerfinder'
include { setup_plasmidfinder; plasmidfinder }                              from './modules/plasmidfinder'
include { barrnap }                                                         from './modules/barrnap'

///////////////////////////////////////////////////////////
//////////////////     Sub-Workflow     ///////////////////
///////////////////////////////////////////////////////////

workflow setup_wf {

    if (params.setEukcc) {
        setup_Eukcc()
    }
    if (params.setBusco) {
        setup_Busco()
    }
    if (params.setGtdbtk) {
        setup_Gtdbtk()
    }
    if (params.setKraken2) {
        setup_Kraken2()
    }
    if (params.setKrona) {
        setup_Krona()
    }
    if (params.setCheckm1) {
        setup_Checkm1()
    }
    if (params.setGunc) {
        setup_Gunc()
    }
    if (params.setCheckm2) {
        setup_Checkm2()
    }
    if (params.setPhyseter) {
        setup_Physeter()
    }
    if (params.setOmark) {
        setup_Omark()
    }
    if (params.setKmerfinder) {
        setup_Kmerfinder()
    }
    if (params.setPlasmidfinder) {
        setup_plasmidfinder()
    }
    if (params.setAll) {
        setup_Eukcc()
        setup_Busco()
        setup_Kraken2()
        setup_Krona()
        setup_Checkm1()
        setup_Gunc()
        setup_Checkm2()
        setup_Physeter()
        setup_Omark() 
        setup_Kmerfinder()
        setup_Gtdbtk()
        setup_plasmidfinder()
    }
}

workflow annotation_wf {
    take:
    assembly

    main:
    prodigal(assembly)
    gffread(prodigal.out)
    omark(gffread.out)
    omark_plot(omark.out.plot.collect())
    busco(gffread.out, "proteins")
    busco_plot(busco.out.plot.collect(), "proteins")

    emit:
    omark.out.report
}

workflow analysis_wf {

    data_file = Channel.fromPath("${params.assemblyFiles}")
    data_dir = Channel.fromPath("${params.assemblyDir}")
    gunc(data_dir)
    busco(data_file, "genome")
    busco_plot(busco.out.plot.collect(), "genome")
    quast(data_file)
    eukcc_folder(data_dir)
    checkm2(data_dir)
    checkm1(data_dir)
    annotation_wf(data_file)
    gtdbtk(data_dir)
    physeter(data_file, params.taxdump)
    krona_physeter(physeter.out.krona, "physeter")
    kraken2(physeter.out.kraken)
    kraken_split(kraken2.out, physeter.out.kraken_split, params.taxdump)
    krona_kraken(kraken2.out, "kraken")
    plasmidfinder(data_file)
    kmerfinder(data_file, "bacteria")
    barrnap(data_file, "bac")

    final_report(
        busco.out.report.collectFile(name: 'busco_multi.report'), \
        quast.out.collectFile(name: 'quast_multi.report'), \
        checkm2.out, \
        gunc.out.report, \
        checkm1.out, \
        eukcc_folder.out, \
        gtdbtk.out, \
        physeter.out.report.collectFile(name: 'physeter_multi.report'), \
        kraken_split.out.collectFile(name: 'kraken2_multi.report'), \
        kmerfinder.out.report.collectFile(name: 'kmerfinder_multi.report'), \
        annotation_wf.out.collectFile(name: 'omark_multi.report')
    )
}

///////////////////////////////////////////////////////
//////////////////     Workflow     ///////////////////
///////////////////////////////////////////////////////

workflow {

    if (params.workflow == "setup") {
        setup_wf()
    }

    else if (params.workflow == "analysis") {
        analysis_wf()
    }

    else {
        println """
        
            Error : Any workflow selected !
        
            """
        helpMessage()
        exit 0 
    }
}

workflow.onComplete = {
    println("Workflow execution completed sucessfully !")
}

workflow.onError = {
    println "Error: something when wrong"
}