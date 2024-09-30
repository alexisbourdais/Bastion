# Bastion

## Overview

ChloroBras is a nextflow pipeline allowing the automatic quality assemment of bacterial assembly. This project is a standalone of the Genome quality assessment workflow from the Genera toolbox of Luc CORNET from the University of Liege: https://github.com/Lcornet/GENERA

## Quick start



### Manual Set-Up

- Busco DB : https://busco-data.ezlab.org/v5/data/lineages/ `busco --download prokaryota` (Prokaryota : 8,3G)

- CheckM1 DB : https://data.ace.uq.edu.au/public/CheckM_databases/

  or `wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz` (275 MB)

- CheckM2 DB : `checkm2 database --download --path /custom/path/` (2.9 GB)

- Eukcc2 DB : `wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz` (11 GB) -> a verif

- Gunc DB : `gunc download_db /path/to/output/dir/` (13 GB)

- Kraken2 DB : `wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20240605.tar.gz` (16 GB version)

- Physeter DB : 

- GTDBTK DB : https://ecogenomics.github.io/GTDBTk/installing/index.html

  or `wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz` (106 GB)

- Original singularity image : https://dox.uliege.be/index.php/s/fUyVuiLwJV0AkB2 (3.4 GB)

- TaxDir : https://figshare.com/articles/software/taxdump/17124074?file=31658918

  `wget https://figshare.com/ndownloader/files/32405687 -O Cornet-Baurain.tgz`  (65.46 MB)

## Parameters

## Documentation

- GENERA : https://github.com/Lcornet/GENERA

- Busco : https://busco.ezlab.org/

- CheckM1 : https://github.com/Ecogenomics/CheckM/wiki

- CheckM2 : https://github.com/chklovski/CheckM2

- Eukcc2 : https://eukcc.readthedocs.io/en/latest/

- GTDBTK : https://ecogenomics.github.io/GTDBTk/

- Gunc : https://grp-bork.embl-community.io/gunc/

- Kraken2 : https://github.com/DerrickWood/kraken2/wiki

- Physeter :  https://metacpan.org/dist/Bio-MUST-Apps-Physeter/view/lib/Bio/MUST/Apps/Physeter/Manual.pod

- Quast : https://quast.sourceforge.net/docs/manual.html

## References
