# BASTION

## Overview

**Bacterial ASsembly contaminaTION** is a nextflow pipeline allowing the automatic quality assemment of bacterial assembly. This project is a standalone of the Genome quality assessment workflow from the **Genera** toolbox of Luc CORNET from the University of Liege: 

**https://github.com/Lcornet/GENERA.**


**Changement from original script :**
- Nextflow dsl 2
- Includes GTDBTK tool
- Final report includes GTDBTK and Eukcc results
- Busco plots, Diamond plots from Gunc results and Krona from Kraken results
- Slurm, Conda, singularity and original usage profile

## Quick start

### Download contams.sif (3.4 Go) : https://dox.uliege.be/index.php/s/fUyVuiLwJV0AkB2 

### Set-Up database

- If you don't have any of the necessary databases, just run `nextflow run XXX.nf -profile X --workflow setup --setAll`
- if you have some database already installed, run `nextflow run XXX.nf -profile X --workflow setup --setX --setY` where x and y correspond to the desired database. For existing databases, please specify their path with --XXX

- Busco_DB (Prokaryota+virus : 8,4 Go)
- CheckM1_DB (1,4Go)
- CheckM2_DB (2.9 Go)
- Eukcc2_DB : (11 Go)
- Gunc_DB (13 Go)
- Kraken2_DB : (16 Go version)
- Physeter_DB (700 Mo)
- Omark (9.4 Go)
- GTDBTK_DB (102 Go)

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
