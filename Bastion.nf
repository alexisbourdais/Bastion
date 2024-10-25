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
    
    --workflow           Select workflow :  'setup' to download database
                                            'analysis' to run all analyses

    Choose which database you want to set-up when using --workflow setup
    --setAll
    --setEukcc
    --setBusco
    --setGtdbtk
    --setKraken2
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
    --db_omark
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

//////////////////////////////////////////////
////////////    Modules     //////////////////
//////////////////////////////////////////////

include { setup_Eukcc; eukcc_folder }   from './modules/eukcc'
include { setup_Checkm1; checkm1 }      from './modules/checkm1'
include { setup_Gunc; gunc }            from './modules/gunc'
include { setup_Kraken2; kraken2 }      from './modules/kraken2'
include { setup_Gtdbtk; gtdbtk }        from './modules/gtdbtk'
include { setup_Physeter; physeter }    from './modules/physeter'
include { setup_Busco; busco }          from './modules/busco'
include { setup_Checkm2; checkm2 }      from './modules/checkm2'
include { setup_Omark; omark }          from './modules/omark'
include { quast }                       from './modules/quast'
include { krona }                       from './modules/krona'
include { prodigal }                    from './modules/prodigal'
include { gffread }                     from './modules/gffread'
include { final_report }                from './modules/report'

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

    if (params.setAll) {
        setup_Eukcc()
        setup_Busco()
        setup_Kraken2()
        setup_Checkm1()
        setup_Gunc()
        setup_Checkm2()
        setup_Physeter()
        //setup_Gtdbtk()
        setup_Omark() 
    }
}

workflow annotation_wf {
    take:
    assembly

    main:
    prodigal(assembly)
    gffread(prodigal.out)
    omark(gffread.out)
}

workflow analysis_wf {

    data_file = Channel.fromPath("${params.assemblyFiles}")
    data_dir = Channel.fromPath("${params.assemblyDir}")
    gunc(data_dir)
    busco(data_dir)
    quast(data_file)
    eukcc_folder(data_dir)
    checkm2(data_dir)
    checkm1(data_dir)
    kraken2(data_file, params.taxdir_kraken)
    krona(kraken2.out.krona)
    physeter(data_file, params.taxdir_physeter)

    if (params.annotation) {
        annotation_wf(data_file)
    }
/*
    if (params.gtdbtk) {
        gtdbtk(data_dir)
    
        final_report(
            busco.out.report, \
            quast.out.collectFile(name: 'quast_multi.report'), \
            checkm2.out, \
            gunc.out.report, \
            checkm1.out, \
            kraken2.out.report.collectFile(name: 'kraken2_multi.report'), \
            physeter.out.collectFile(name: 'physeter_multi.report'), \
            eukcc_folder.out, \
            gtdbtk.out.report
        )
    }
 /*   
    else {
        final_report(
            busco.out.report, \
            quast.out.collectFile(name: 'quast_multi.report'), \
            checkm2.out, \
            gunc.out.report, \
            checkm1.out, \
            kraken2.out.report.collectFile(name: 'kraken2_multi.report'), \
            physeter.out.collectFile(name: 'physeter_multi.report'), \
            eukcc_folder.out
        )
    }
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