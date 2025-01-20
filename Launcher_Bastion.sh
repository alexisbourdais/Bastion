#!  /usr/bin/bash

nextflow run Bastion.nf \
-profile slurm,singularity \
--workflow analysis \
--singularity "-B /home:/home"
