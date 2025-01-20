#!/usr/bin/env nextflow

/*
===============================================================
 This script is largely based on the Genome quality assessment workflow from the Genera toolbox of Luc CORNET from the University of Liege:
 https://github.com/Lcornet/GENERA
 
 ---------------------------------------------------------------
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

    Command : nextflow run Bastion.nf -profile slurm,singularity [option]

    REQUIRED parameter

    -profile [standard]/slurm,      Select profile standard (local) or slurm. Default: standard          
             singularity/conda      Select profile singularity or conda. 
                                    Physeter/Kraken/Report need singularity in both case.
    
    --workflow           Select workflow :  'setup' to download database
                                            'analysis' to run all analyses

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

    Database directory : automatic if installed with --workflow setup --setAll
    --db_busco          Path to database directory
    --db_checkm2        Path to checkm2_uniref100.KO.1.dmnd
    --db_eukcc          Path to database directory
    --db_gunc           Path to gunc_db_progenomes2.1.dmnd
    --db_kraken2        Path to database directory
    --db_gtdbtk         Path to database directory
    --db_checkm1        Path to database directory
    --db_physeter       Path to life-tqmd-of73.dmnd
    --db_omark          Path to database directory
    --taxdir            Path to taxdump

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
    --lineage_busco         [auto-lineage], bacteria_odb10, fungi_odb10 ...

    Eukcc2
    --mode_eukcc            [DNA] ...

    Physeter/Kraken
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

include { setup_Eukcc; eukcc_folder }       from './modules/eukcc'
include { setup_Checkm1; checkm1 }          from './modules/checkm1'
include { setup_Gunc; gunc }                from './modules/gunc'
include { setup_Kraken2; kraken2 }          from './modules/kraken2'
include { setup_Gtdbtk; gtdbtk }            from './modules/gtdbtk'
include { setup_Physeter; physeter }        from './modules/physeter'
include { setup_Busco; busco; busco_plot }  from './modules/busco'
include { setup_Checkm2; checkm2 }          from './modules/checkm2'
include { setup_Omark; omark; omark_plot }  from './modules/omark'
include { quast }                           from './modules/quast'
include { setup_Krona; krona }              from './modules/krona'
include { prodigal }                        from './modules/prodigal'
include { gffread }                         from './modules/gffread'
include { final_report }                    from './modules/report'
include { setup_Kmerfinder; kmerfinder }    from './modules/kmerfinder'

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

    if (params.setAll) {
        setup_Eukcc()
        setup_Busco()
        setup_Kraken2()
        setup_Krona()
        setup_Checkm1()
        setup_Gunc()
        setup_Checkm2()
        setup_Physeter()
        //setup_Gtdbtk()
        setup_Omark() 
        setup_Kmerfinder()
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

    emit:
    omark.out.report
}

workflow analysis_wf {

    data_file = Channel.fromPath("${params.assemblyFiles}")
    data_dir = Channel.fromPath("${params.assemblyDir}")
    //gunc(data_dir)
    //busco(data_file)
    //busco_plot(busco.out.plot.collect())
    //quast(data_file)
    //eukcc_folder(data_dir)
    //checkm2(data_dir)
    //checkm1(data_dir)
    //kraken2(data_file)
    //krona(kraken2.out.krona)
    physeter(data_file, params.taxdir_physeter)
    //annotation_wf(data_file)
    //gtdbtk(data_dir)
    kmerfinder(data_file)
/*
    final_report(
        busco.out.report.collectFile(name: 'busco_multi.report'), \
        quast.out.collectFile(name: 'quast_multi.report'), \
        checkm2.out, \
        gunc.out.report, \
        checkm1.out, \
        eukcc_folder.out, \
        gtdbtk.out
        //annotation_wf.out.collectFile(name: 'omark_multi.report')
        //kmerfinder.out
        //kraken2.out.report.collectFile(name: 'kraken2_multi.report'), \
        //physeter.out.collectFile(name: 'physeter_multi.report')
    )
*/    
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

workflow.onComplete{
    if (workflow.success) {
        println("Workflow execution completed sucessfully !")
    } 
    else {
        println("Workflow execution completed with errors !")
    }
}
