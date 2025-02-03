#!  /usr/bin/bash

source /local/env/envnextflow-24.10.sh

nextflow run Bastion.nf \
-profile slurm,singularity \
--singularity "-B /home:/home -B /scratch/abourdais:/scratch/abourdais -B /projects:/projects" \
--workflow analysis \
--db_busco /scratch/abourdais/Database/busco_downloads/ \
--db_kraken2 /scratch/abourdais/Database/kraken2_db/ \
--db_gtdbtk /scratch/abourdais/Database/gtdbtk_db/release220/ \
--db_kmerfinder /scratch/abourdais/Database/kmerfinder_db/

#--format "fna" \
