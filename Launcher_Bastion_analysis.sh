#!  /usr/bin/bash

nextflow run Bastion.nf \
-profile slurm,singularity \
--singularity "-B /home:/home" \
--workflow analysis \
#--format "fna"
