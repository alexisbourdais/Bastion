# BASTION

## Overview

**Bacterial ASsembly contaminaTION** is a nextflow pipeline allowing the automatic quality assemment of bacterial assembly. This project is a standalone of the Genome quality assessment workflow from the **Genera** toolbox of Luc CORNET from the University of Liege: 

**https://github.com/Lcornet/GENERA.**


**Changement from original script :**
- Nextflow dsl 2
- Slurm and Conda or singularity profiles
- Includes GTDBTK, Kmerfinder and Omark tools
- Busco plots, Diamond plots from Gunc results and Krona from Kraken/Physeters results
- Annotation with Prodigal

## Quick start

### Set-Up database

- If you don't have any of the necessary databases, just run `nextflow run Bastion.nf -profile slurm,singularity --workflow setup --setAll`
- if you have some database already installed, run `nextflow run Bastion.nf -profile slurm,singularity --workflow setup --setBusco --setKraken2` in order to install Busco and Kraken2 database for exemple

- Busco DB      : 8,4 Go (Prokaryota+virus)
- CheckM1 DB    : 1,4 Go
- CheckM2 DB    : 2.9 Go
- Eukcc2 DB     : 11  Go
- Gunc DB       : 13  Go
- Kraken2 DB    : 16  Go (PlusPF-16 version)
- Krona DB      : ??  Go
- Physeter DB   : 700 Mo
- Omark DB      : 9.4 Go
- GTDBTK DB     : 102 Go
- Kmerfinder DB : ??  Go

Total = ?? Go

### Run analysis

  To completed

## Parameters

nextflow run Bastion.nf --help

    Command : nextflow run Bastion.nf -profile slurm,singularity [option]

    REQUIRED parameter

    -profile [standard]/slurm,      Select profile standard (local) or slurm. Default: standard          
             singularity/conda      Select profile singularity or conda. 
                                    Physeter need singularity in both case.
    
    --workflow                      Select workflow :  'setup' to download database
                                                       'analysis' to run all analyses

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

    Database directory : automatic if installed with --workflow setup
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
    --taxdir            Path to taxdump
    --db_kmerfinder     Path to database directory

    OPTIONAL parameter
    
    Assembly directory
    --assemblyDir            Default: "./Data/"
    --format                 Default: "fasta"

    Results directory
    --resultsDir            Path to results directory, default: "./Results/"

    Quast
    --mode_quast            [bacteria], eukaryote, fungus

    Busco
    --lineage_busco         [auto-lineage], bacteria_odb10, fungi_odb10 ...

    Physeter
    --taxlevel              [phylum]
    --automode              [label_first]

    Kmerfinder
    --taxon_kmerfinder      [bacteria], archae, fungi, protozoa

## Documentation

- GENERA : https://github.com/Lcornet/GENERA

- Busco : https://busco.ezlab.org/

- CheckM1 : https://github.com/Ecogenomics/CheckM/wiki

- CheckM2 : https://github.com/chklovski/CheckM2

- Eukcc2 : https://eukcc.readthedocs.io/en/latest/

- GTDBTK : https://ecogenomics.github.io/GTDBTk/

- Gffread : https://github.com/gpertea/gffread

- Gunc : https://grp-bork.embl-community.io/gunc/

- Kmerfinder : https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/

- Kraken2 : https://github.com/DerrickWood/kraken2/wiki

- Krona : https://github.com/marbl/Krona
  
- OMArk : https://github.com/DessimozLab/OMArk

- Physeter :  https://metacpan.org/dist/Bio-MUST-Apps-Physeter/view/lib/Bio/MUST/Apps/Physeter/Manual.pod

- Prodigal : https://github.com/hyattpd/Prodigal

- Quast : https://quast.sourceforge.net/docs/manual.html


## References

Cornet L, Durieu B, Baert F, D’hooge E, Colignon D, Meunier L, Lupo V, Cleenwerck I, Daniel H-M, Rigouts L, Sirjacobs D, Declerck S, Vandamme P, Wilmotte A, Baurain D, Becker P (2022). The GEN-ERA toolbox: unified and reproducible workflows for research in microbial genomics. Submitted to GIGAscience

Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes.  Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654

Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2014. Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043-1055

Saary, Paul, Alex L. Mitchell, and Robert D. Finn. “Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC.” Genome biology 21.1 (2020): 1-21.

Pierre-Alain Chaumeil, Aaron J Mussig, Philip Hugenholtz, Donovan H Parks, GTDB-Tk v2: memory friendly classification with the genome taxonomy database, Bioinformatics, Volume 38, Issue 23, 1 December 2022, Pages 5315–5316, https://doi.org/10.1093/bioinformatics/btac672

Orakov, A., Fullam, A., Coelho, L.P. et al. GUNC: detection of chimerism and contamination in prokaryotic genomes. Genome Biol 22, 178 (2021). https://doi.org/10.1186/s13059-021-02393-0

Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). https://doi.org/10.1186/s13059-019-1891-0

Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.

OMArk, a tool for gene annotation quality control, reveals erroneous gene inference. Nat Biotechnol 43, 40–41 (2025). https://doi.org/10.1038/s41587-024-02155-w

Alla Mikheenko, Andrey Prjibelski, Vladislav Saveliev, Dmitry Antipov, Alexey Gurevich, Versatile genome assembly evaluation with QUAST-LG, Bioinformatics, Volume 34, Issue 13, July 2018, Pages i142–i150, https://doi.org/10.1093/bioinformatics/bty266

Benchmarking of Methods for Genomic Taxonomy. Larsen MV, Cosentino S, Lukjancenko O, Saputra D, Rasmussen S, Hasman H, Sicheritz-Pontén T, Aarestrup FM, Ussery DW, Lund O. J Clin Microbiol. 2014 Feb 26. [Epub ahead of print]

Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119
